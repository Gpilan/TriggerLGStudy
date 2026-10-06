#!/usr/bin/env python3
"""Extract pencil-beam timing from an audited input plan and plot length/thickness scans.

Use the plan's frozen analysis modules, never the current simulation build.
ROOT performs all Gaussian fits and renders the delta-t distributions.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def extract_condition(task):
    condition, module_dir, out = task
    sys.path.insert(0, module_dir)
    from verified_timing import event_metrics
    rows, seen = [], set()
    for source in condition['inputs']:
        path = Path(source['canonical'])
        if sha(path) != source['sha256']:
            raise ValueError(f'Canonical input changed: {path}')
        events = json.loads(path.read_text())
        if len(events) != source['events']:
            raise ValueError('Wrong input event count')
        for event in events:
            identity = (source['seed'], event['event'])
            if identity in seen:
                raise ValueError('Duplicate seed/event')
            seen.add(identity)
            b = event['budget']
            if sum(v for k, v in b.items() if k.startswith('fate_')) != b['started']:
                raise ValueError('Photon budget does not close')
            if sum(event['detected']) != b['fate_detected']:
                raise ValueError('Detected count mismatch')
            metric = event_metrics(event, dt=.01, fwhm=3., fraction=.3, half=2)
            rows.append(dict(seed=source['seed'], event=event['event'],
                             detected=event['detected'], generated=b['started'],
                             fates={k[5:]: v for k, v in b.items() if k.startswith('fate_')},
                             times_ns=[t['toa_ns'] for t in metric['triggers']],
                             delta_ns=metric['delta_ns'],
                             quality=[t['quality']['method'] for t in metric['triggers']]))
    if len(rows) != 3000:
        raise ValueError('Expected 3000 incident events per condition')
    result = {k: v for k, v in condition.items() if k != 'inputs'}
    result.update(events=rows, inputs=condition['inputs'])
    Path(out, condition['condition']+'.json').write_text(json.dumps(result, allow_nan=False))
    print('EXTRACTED', condition['condition'], len(rows), flush=True)
    return condition['condition']


def extract(plan_path, out, workers):
    plan = json.loads(plan_path.read_text())
    for p, digest in plan['timing_source_sha256'].items():
        if sha(p) != digest:
            raise ValueError('Frozen timing module changed: '+p)
    data = out/'data'
    data.mkdir(parents=True, exist_ok=True)
    tasks = [(c, plan['timing_module_dir'], str(data)) for c in plan['conditions']]
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        list(pool.map(extract_condition, tasks))
    (data/'extraction-audit.json').write_text(json.dumps(dict(
        conditions=len(tasks), incident_events=plan['total_events'],
        input_plan_sha256=sha(plan_path), CFD_fraction=.3, sample_ns=.01,
        SPE_fwhm_ns=3., fit_half_window=2, LED_cut=False,
        timing_source_sha256=plan['timing_source_sha256'],
        outputs_sha256={c['condition']+'.json':sha(data/(c['condition']+'.json'))
                        for c in plan['conditions']}), indent=2))


def plot(out):
    import numpy as np
    import ROOT
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from collections import Counter
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, 'XYZ')
    ROOT.gStyle.SetTitleFont(42, 'XYZ')
    data, plots = out/'data', out/'plots'
    plots.mkdir(exist_ok=True)
    audit = json.loads((data/'extraction-audit.json').read_text())
    for name, digest in audit['outputs_sha256'].items():
        if sha(data/name) != digest:
            raise ValueError('Extraction changed: '+name)
    sets = [json.loads((data/name).read_text()) for name in audit['outputs_sha256']]
    values = {d['condition']:np.array([e['delta_ns'] for e in d['events']
                                     if e['delta_ns'] is not None]) for d in sets}
    assert len(sets) == 27 and all(len(d['events']) == 3000 for d in sets)
    assert all(len(v) > 3 and np.isfinite(v).all() for v in values.values())
    all_values = np.concatenate(list(values.values()))
    bw = .005
    lo = math.floor((all_values.min()-.02)/bw)*bw
    hi = math.ceil((all_values.max()+.02)/bw)*bw
    nb = round((hi-lo)/bw)
    edges = np.linspace(lo, hi, nb+1)
    archive = ROOT.TFile(str(data/'dt_fits.root'), 'RECREATE')
    keep, objects, summaries = [], {}, []
    def fit(v, tag, save=True):
        h = ROOT.TH1D('h_'+tag, '', nb, lo, hi)
        h.SetDirectory(0)
        for value in v:
            h.Fill(float(value))
        assert h.GetBinContent(0) == h.GetBinContent(nb+1) == 0
        f = ROOT.TF1('gaussian_'+tag, 'gaus', lo, hi)
        f.SetParameters(h.GetMaximum(), float(v.mean()), float(v.std()))
        f.SetParLimits(0, 1e-8, 1e8)
        f.SetParLimits(1, lo, hi)
        f.SetParLimits(2, 1e-4, hi-lo)
        result = h.Fit(f, 'LQRSNI', '', lo, hi)
        if int(result) != 0 or result.CovMatrixStatus() != 3:
            raise ValueError(f'Fit failed: {tag}, {int(result)}, {result.CovMatrixStatus()}')
        expected = np.maximum(np.array([f.Integral(float(a), float(b))/bw
                                       for a,b in zip(edges[:-1],edges[1:])]), 1e-300)
        observed = np.histogram(v, edges)[0]
        term = expected-observed
        mask = observed>0
        term[mask] += observed[mask]*np.log(observed[mask]/expected[mask])
        stats = dict(sigma_ps=f.GetParameter(2)*1000, sigma_error_ps=f.GetParError(2)*1000,
                     mean_ns=f.GetParameter(1), poisson_deviance=float(2*term.sum()),
                     ndf=nb-3, fit_status=int(result), covariance_status=int(result.CovMatrixStatus()))
        if save:
            x = np.linspace(lo, hi, 20001)
            curve = ROOT.TGraph(len(x), x, np.array([f.Eval(float(t)) for t in x]))
            curve.SetName('curve_'+tag)
            curve.SetLineColor(ROOT.kRed+1)
            curve.SetLineWidth(2)
            archive.cd()
            for obj in [h,f,curve]: obj.Write()
            result.Get().Write('fit_'+tag)
            objects[tag] = (h,curve,stats)
            keep.extend([h,f,curve,result])
        return stats
    def mean_sem(v):
        a = np.asarray(v, dtype=float)
        return float(a.mean()), float(a.std(ddof=1)/np.sqrt(len(a)))
    sensitivity = None
    for d in sets:
        tag = d['condition']
        v = values[tag]
        row = {k:d[k] for k in ['condition','mode','width_mm','length_mm','thickness_mm']}
        row.update(fit(v,tag), incident=len(d['events']), valid_pairs=len(v))
        counts = np.array([e['detected'] for e in d['events']])
        generated = np.array([e['generated'] for e in d['events']])
        total = counts.sum(axis=1)
        for label,a in [('detected',total),('detected_T1',counts[:,0]),('detected_T2',counts[:,1]),('generated',generated)]:
            row[label+'_mean'],row[label+'_sem'] = mean_sem(a)
        ratio = float(total.sum()/generated.sum())
        row['collection_fraction'] = ratio
        row['collection_fraction_sem'] = float(np.std(total-ratio*generated,ddof=1)/np.sqrt(len(total))/generated.mean())
        fates = Counter()
        for event in d['events']:fates.update(event['fates'])
        row['no_rindex'] = fates['no_rindex']
        row['fates'] = dict(fates)
        assert sum(fates.values()) == int(generated.sum())
        if row['no_rindex']:
            assert tag=='lg_W10_L40_T5' and row['no_rindex']==1
            bad=[e for e in d['events'] if e['fates'].get('no_rindex',0)]
            assert len(bad)==1 and (bad[0]['seed'],bad[0]['event'])==(106600015,10)
            clean=np.array([e['delta_ns'] for e in d['events'] if e['delta_ns'] is not None and not e['fates'].get('no_rindex',0)])
            other=fit(clean,tag+'_exclude_exception',save=False)
            sensitivity={'condition':tag,'retained':row.copy(),'exclude_one_event':other,
                         'sigma_difference_ps':other['sigma_ps']-row['sigma_ps'],
                         'interpretation':'Sensitivity to excluding the affected event; not a corrected-photon simulation.'}
        summaries.append(row)
        print('FIT',tag,len(v),round(row['sigma_ps'],3),flush=True)
    groups=[('nolg',10),('lg',10),('lg',40)]
    colors=['#0072B2','#009E73','#D55E00'];markers=['o','s','^']
    label=lambda m,w:('noLG' if m=='nolg' else 'LG')+f' | W = {w} mm'
    def scan_rows(axis,mode,width):
        return sorted([r for r in summaries if r['mode']==mode and r['width_mm']==width
                       and (r['thickness_mm']==5 if axis=='length' else r['length_mm']==60)],key=lambda r:r[axis+'_mm'])
    def text(x,y,s,size=.04,color=1):
        t=ROOT.TLatex();t.SetNDC();t.SetTextFont(42);t.SetTextSize(size);t.SetTextColor(color)
        t.DrawLatex(x,y,s);keep.append(t)
    for axis in ['length','thickness']:
        unit='L' if axis=='length' else 'T'
        fixed='T = 5 mm' if axis=='length' else 'L = 60 mm'
        ticks=[30,40,50,60,80] if axis=='length' else [2,3,5,7,10]
        xlabel=f'Scintillator {axis} {unit} [mm]'
        for mode,width in groups:
            c=ROOT.TCanvas(f'dt_{axis}_{mode}_w{width}','',2250,1400);keep.append(c);c.Divide(3,2,.003,.004)
            rr=scan_rows(axis,mode,width)
            for i,r in enumerate(rr):
                pad=c.cd(i+1);pad.SetLeftMargin(.14);pad.SetBottomMargin(.15);pad.SetTopMargin(.1);pad.SetRightMargin(.03);pad.SetTicks(1,1)
                h,curve,stats=objects[r['condition']]
                draw=h.Clone('draw_'+pad.GetName());draw.SetDirectory(0);keep.append(draw)
                draw.SetLineColor(ROOT.kBlack);draw.SetLineWidth(2);draw.SetMinimum(0);draw.SetMaximum(h.GetMaximum()*1.45)
                draw.GetXaxis().SetTitle('#Deltat = T1 - T2 [ns]');draw.GetYaxis().SetTitle('Events / 5 ps')
                for ax in [draw.GetXaxis(),draw.GetYaxis()]:ax.SetLabelSize(.042);ax.SetTitleSize(.05)
                draw.GetYaxis().SetTitleOffset(1.35);draw.Draw('HIST');curve.Draw('L SAME')
                text(.15,.94,f'{label(mode,width)} | {unit} = {r[axis+"_mm"]} mm',.045)
                text(.17,.84,f'#sigma = {r["sigma_ps"]:.2f} #pm {r["sigma_error_ps"]:.2f} ps',.045,ROOT.kRed+1)
                text(.17,.77,f'N = {r["valid_pairs"]}',.036)
                text(.17,.70,f'Poisson D/ndf = {r["poisson_deviance"]:.1f}/{r["ndf"]}',.032,ROOT.kRed+1)
                pad.RedrawAxis()
            c.cd(6);text(.12,.72,label(mode,width),.07);text(.12,.58,fixed,.06)
            text(.12,.43,'Gaussian fit; 5 ps bins',.05);text(.12,.30,'CFD = 0.30; no LED cut',.05)
            c.SaveAs(str(plots/f'dt_{axis}_{mode}_w{width}.png'));archive.cd();c.Write()
        fig,axs=plt.subplots(1,2,figsize=(14,5.5))
        for (m,w),color,marker in zip(groups,colors,markers):
            rr=scan_rows(axis,m,w);xx=[r[axis+'_mm'] for r in rr]
            for ax,key,err in [(axs[0],'sigma_ps','sigma_error_ps'),(axs[1],'detected_mean','detected_sem')]:
                ax.errorbar(xx,[r[key] for r in rr],yerr=[r[err] for r in rr],color=color,marker=marker,lw=2,ms=7,capsize=4,label=label(m,w))
        for ax in axs:ax.set_xlabel(xlabel);ax.set_xticks(ticks);ax.grid(alpha=.25);ax.legend(fontsize=10)
        axs[0].set_ylabel(r'Gaussian fit $\sigma(\Delta t)$ [ps]')
        axs[1].set_ylabel('Mean detected photons / event (T1 + T2)')
        axs[0].set_title('Pair timing');axs[1].set_title('Detected photons')
        fig.suptitle(f'Tile-{axis} dependence | {fixed}',fontsize=17)
        fig.text(.08,.025,'3000 events / condition | timing: fit SE; photon count: SEM',fontsize=10)
        fig.tight_layout(rect=[0,.065,1,.95]);fig.savefig(plots/f'sigma_and_detected_vs_{axis}.png',dpi=180);plt.close(fig)
        fig,axs=plt.subplots(1,2,figsize=(14,5.5))
        for (m,w),color,marker in zip(groups,colors,markers):
            rr=scan_rows(axis,m,w);xx=[r[axis+'_mm'] for r in rr]
            for ax,key,err in [(axs[0],'generated_mean','generated_sem'),(axs[1],'collection_fraction','collection_fraction_sem')]:
                ax.errorbar(xx,[r[key] for r in rr],yerr=[r[err] for r in rr],color=color,marker=marker,lw=2,ms=7,capsize=4,label=label(m,w))
        for ax in axs:ax.set_xlabel(xlabel);ax.set_xticks(ticks);ax.grid(alpha=.25);ax.legend(fontsize=10)
        axs[0].set_ylabel('Mean generated optical photons / event')
        axs[1].set_ylabel('Detected / generated photons (T1 + T2)')
        fig.suptitle(f'Photon yield and collection | {fixed}',fontsize=17)
        fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(plots/f'photon_yield_vs_{axis}.png',dpi=180);plt.close(fig)
        categories=[('detected','Detected','#2ca02c'),('qe_reject','QE reject','#98df8a'),('absorb_scint','Scint. absorption','#ff7f0e'),('absorb_lg','LG bulk','#1f77b4'),('absorb_lg_surface','LG surface','#17becf'),('absorb_glass','Glass','#aec7e8'),('absorb_other_surface','Other surface','#d62728'),('absorb_other','Other bulk','#8c564b'),('no_rindex','NoRINDEX','#000000')]
        fig,axs=plt.subplots(1,3,figsize=(16,5.5),sharey=True)
        for ax,(m,w) in zip(axs,groups):
            rr=scan_rows(axis,m,w);bottom=np.zeros(len(rr));indices=np.arange(len(rr))
            for key,lab,color in categories:
                vv=np.array([r['fates'].get(key,0)/sum(r['fates'].values()) for r in rr])
                ax.bar(indices,vv,bottom=bottom,color=color,label=lab);bottom+=vv
            assert np.allclose(bottom,1,rtol=0,atol=1e-12),'Unplotted fate category'
            ax.set_xticks(indices,[r[axis+'_mm'] for r in rr]);ax.set_xlabel(xlabel);ax.set_title(label(m,w));ax.set_ylim(0,1.02)
        axs[0].set_ylabel('Fraction of generated optical photons')
        handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.055),ncol=5,fontsize=9)
        fig.text(.5,.012,'Other surface includes Al/tape outside the LG-surface category; other bulk is the remaining bulk absorption.',ha='center',fontsize=9)
        fig.suptitle(f'Optical photon fates | {fixed}',fontsize=17);fig.tight_layout(rect=[0,.19,1,.95]);fig.savefig(plots/f'fate_fractions_vs_{axis}.png',dpi=180);plt.close(fig)
    archive.Close()
    check=ROOT.TFile.Open(str(data/'dt_fits.root'))
    for r in summaries:
        tag=r['condition'];assert check.Get('h_'+tag).GetEntries()==r['valid_pairs']
        result=check.Get('fit_'+tag);assert result.Status()==0 and result.CovMatrixStatus()==3
        f,g=check.Get('gaussian_'+tag),check.Get('curve_'+tag)
        assert g.GetN()==20001
        for i in range(g.GetN()):assert math.isclose(g.GetPointY(i),f.Eval(g.GetPointX(i)),rel_tol=1e-12,abs_tol=1e-12)
    check.Close()
    (data/'summary.json').write_text(json.dumps(summaries,indent=2))
    flat=[{k:v for k,v in r.items() if k!='fates'} for r in summaries]
    with (data/'summary.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]),lineterminator='\n');writer.writeheader();writer.writerows(flat)
    (data/'exception-sensitivity.json').write_text(json.dumps(sensitivity,indent=2))
    (data/'plot-audit.json').write_text(json.dumps(dict(status='validated_with_documented_optical_exception',
        conditions=len(summaries),incident_events=sum(r['incident'] for r in summaries),
        valid_pairs=sum(r['valid_pairs'] for r in summaries),no_rindex=sum(r['no_rindex'] for r in summaries),
        ROOT_version=ROOT.gROOT.GetVersion(),fit='ROOT TH1::Fit Gaussian LQRSNI; full common range',
        bin_ns=bw,fit_range_ns=[lo,hi],sqrt2_conversion=False,
        selection='Both CFD times valid; no LED cut; photon statistics use all incident events',
        curve_points=20001,outputs_sha256={str(p.relative_to(out)):sha(p) for p in [*plots.glob('*.png'),data/'summary.csv',data/'dt_fits.root']}
        ),indent=2))
    print('PLOTS_VALIDATED',len(summaries),len(list(plots.glob('*.png'))),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['extract','plot'])
    p.add_argument('--plan',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2)
    a=p.parse_args()
    if a.stage=='extract':
        if a.plan is None:p.error('--plan is required for extraction')
        extract(a.plan,a.out,a.workers)
    else:plot(a.out)


if __name__=='__main__':main()
