#!/usr/bin/env python3
"""Build a local results portal from an explicitly selected JSON catalog.
Never moves source artifacts or runs analysis. Links are relative; no PNG copies.
"""
from pathlib import Path
import argparse, json, os, html, hashlib
from urllib.parse import quote

def build(repo, catalog):
    config=json.loads(catalog.read_text()); records=[]; cards=[]
    root=repo/'analysis/studies'; root.mkdir(parents=True,exist_ok=True)
    def url(target, parent): return quote(os.path.relpath(target,parent),safe='/')
    def link(path,target):
        target=target.resolve(strict=True)
        if path.is_symlink():
            if path.resolve()==target:return
            path.unlink()
        elif path.exists():raise RuntimeError('Refusing to replace regular file: '+str(path))
        path.symlink_to(os.path.relpath(target,path.parent),target_is_directory=target.is_dir())
    style='body{font:16px system-ui;margin:auto;max-width:1250px;padding:30px;background:#f4f6fa;color:#182434}nav{display:flex;gap:20px;flex-wrap:wrap}a{color:#086ea8}article,section{background:white;border:1px solid #dae0e8;border-radius:12px;padding:22px;margin:22px 0}img{max-width:100%;height:auto}h1{font-size:32px}h2{font-size:24px}.meta{color:#526275}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:20px}summary{cursor:pointer;font-weight:600}figure{margin:20px 0}figcaption{overflow-wrap:anywhere}'
    def page(title,body):return '<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+style+'</style><body>'+body+'</body></html>'
    for study in config['studies']:
        slug=study['id']; assert slug.replace('_','').isalnum()
        out=root/slug;out.mkdir(exist_ok=True);(out/'figures').mkdir(exist_ok=True)
        source=repo/study['source']; assert source.is_dir()
        link(out/'records',source)
        if (source/'tools').exists():link(out/'tools',source/'tools')
        if (source/'data').exists():link(out/'data',source/'data')
        previous=json.loads((out/'manifest.json').read_text()) if (out/'manifest.json').exists() else {'figures':[]}
        groups={}; figs=[]
        for selection in study['selections']:
            for f in sorted((repo/selection['folder']).glob(selection['glob'])):
                name=selection.get('prefix','')+f.name
                assert name not in [x['name'] for x in figs],name
                link(out/'figures'/name,f)
                record={'name':name,'source':str(f.relative_to(repo)),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'group':selection['group']}
                figs.append(record);groups.setdefault(record['group'],[]).append(record)
        assert figs and study['hero'] in [x['name'] for x in figs]
        names={x['name'] for x in figs}
        for item in previous['figures']:
            stale=out/'figures'/item['name']
            if item['name'] not in names:
                if not stale.is_symlink():raise RuntimeError('Refusing to remove non-link: '+str(stale))
                stale.unlink()
        title=html.escape(study['title']);meta=html.escape(study['summary']);note=html.escape(study['model_note'])
        intro=f'<nav><a href="../../../RESULTS.html">전체 결과</a><a href="../MODEL_STATUS.md">시뮬레이션 반영 상태</a><a href="records/README.md">상세 기록</a></nav><h1>{title}</h1><p class="meta">{meta}</p><p>{note}</p>'
        body=intro
        for group,items in groups.items():
            body+='<section><details open><summary>'+html.escape(group)+f' · {len(items)}개</summary>'
            for item in items:
                u='figures/'+quote(item['name']);body+=f'<figure><a href="{u}"><img loading="lazy" src="{u}" alt="{html.escape(item["name"])}"></a><figcaption><a href="{u}">{html.escape(item["name"])}</a></figcaption></figure>'
            body+='</details></section>'
        (out/'index.html').write_text(page(study['title'],body))
        (out/'README.md').write_text(f'# {study["title"]}\n\n{study["summary"]}\n\n{study["model_note"]}\n\n[그림 전체](index.html) · [최신 PNG](figures/) · [원본 연구 기록](records/README.md) · [전체 결과](../../../RESULTS.html) · [모델 반영 상태](../MODEL_STATUS.md)\n\n그림은 검증된 원본으로 연결된 상대 심볼릭 링크입니다. 이 폴더의 링크를 통해 PNG를 덮어쓰지 마세요. 새 분석은 새 기록 폴더에서 수행하고 catalog의 선택을 갱신합니다. 원본 ROOT·로그·분석 도구 경로는 유지합니다.\n')
        (out/'manifest.json').write_text(json.dumps({'study':study,'figures':figs},ensure_ascii=False,indent=2)+'\n')
        href=f'analysis/studies/{slug}/index.html';hero=f'analysis/studies/{slug}/figures/{study["hero"]}'
        cards.append(f'<article><h2><a href="{href}">{title}</a></h2><p class="meta">{meta}</p><a href="{href}"><img src="{hero}" alt="{title}"></a><p>{note}</p><a href="{href}">그림 {len(figs)}개 보기 →</a></article>')
        records.extend(figs)
    body='<h1>Trigger · 최신 연구 결과</h1><p class="meta">갱신 '+html.escape(config['updated'])+' · 검증된 결과만 선택 · PNG 원본 연결</p><nav><a href="analysis/studies/MODEL_STATUS.md">현재 시뮬레이션과 결과 버전</a><a href="analysis/tools/README.md">분석 도구 안내</a><a href="README.md">프로젝트</a></nav>'+''.join(cards)+'<footer>개인 로컬 결과 포털입니다. 과거 실행 데이터와 기록은 원래 경로에 보존됩니다.</footer>'
    (repo/'RESULTS.html').write_text(page('Trigger 최신 결과',body))
    (root/'README.md').write_text('# 연구 결과\n\n[최신 결과 대시보드](../../RESULTS.html) · [모델 상태](MODEL_STATUS.md)\n\n'+''.join(f'- [{s["title"]}]({s["id"]}/README.md)\n' for s in config['studies'])+'\n날짜별 원본은 이동하지 않고 최신 검증 결과를 주제별로 연결합니다. 이 로컬 카탈로그는 Git에서 제외합니다.\n')
    return records
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--catalog',type=Path,required=True);args=a.parse_args()
    repo=Path(__file__).resolve().parents[2];r=build(repo,args.catalog.resolve());print(f'Linked {len(r)} figures; wrote RESULTS.html')
