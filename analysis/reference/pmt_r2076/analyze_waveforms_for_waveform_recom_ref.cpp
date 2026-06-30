#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct Parameters {
  double responseRiseNs = 1.0;
  double responseFwhmNs = 3.0;
  double fraction = 0.30;
  double sampleStepNs = 0.20;
  // Gamma-pulse shape and scale calibrated so its continuous 10-90% rise and
  // FWHM equal responseRiseNs and responseFwhmNs for the PDF defaults.
  double gammaShape = 1.915604733026;
  double gammaTauNs = 0.902234354901;
};

struct Pulse {
  std::vector<double> time;
  std::vector<double> voltage;
  double photonCount = 0.0;
  double peak = 0.0;
  double t10 = std::numeric_limits<double>::quiet_NaN();
  double t30 = std::numeric_limits<double>::quiet_NaN();
  double t90 = std::numeric_limits<double>::quiet_NaN();
  double photonT30 = std::numeric_limits<double>::quiet_NaN();
  double arrivalMean = std::numeric_limits<double>::quiet_NaN();
  double arrivalRms = std::numeric_limits<double>::quiet_NaN();
  double arrivalFirst = std::numeric_limits<double>::quiet_NaN();
  double arrivalLast = std::numeric_limits<double>::quiet_NaN();
  double arrivalQ10 = std::numeric_limits<double>::quiet_NaN();
  double arrivalQ50 = std::numeric_limits<double>::quiet_NaN();
  double arrivalQ90 = std::numeric_limits<double>::quiet_NaN();
};

struct EventResult {
  int event = -1;
  double energyDeposit = 0.0;
  double generatedOpticalPhotons = 0.0;
  Pulse trigger[2];
};

struct Dataset {
  std::string label;
  std::string path;
  std::vector<EventResult> events;
  Pulse example[2];
};

struct Summary {
  size_t n = 0;
  double mean = 0.0;
  double sigma = 0.0;
  double median = 0.0;
  double half68 = 0.0;
  double minimum = 0.0;
  double maximum = 0.0;
};

std::pair<double, double> gammaPulseMetrics(double shape) {
  auto response = [shape](double x) {
    if (x <= 0.0) return 0.0;
    return std::pow(x / shape, shape) * std::exp(shape - x);
  };
  auto root = [&](double level, double low, double high, bool rising) {
    for (int i = 0; i < 100; ++i) {
      const double middle = 0.5 * (low + high);
      if ((response(middle) < level) == rising) low = middle;
      else high = middle;
    }
    return 0.5 * (low + high);
  };
  const double t10 = root(0.10, 0.0, shape, true);
  const double t90 = root(0.90, 0.0, shape, true);
  const double halfLeft = root(0.50, 0.0, shape, true);
  const double halfRight = root(0.50, shape, shape + 100.0, false);
  return {t90 - t10, halfRight - halfLeft};
}

void calibrateGammaPulse(Parameters& p) {
  const double targetRatio = p.responseFwhmNs / p.responseRiseNs;
  double low = 0.05;
  double high = 50.0;
  const double maxRatio = gammaPulseMetrics(low).second /
                          gammaPulseMetrics(low).first;
  const double minRatio = gammaPulseMetrics(high).second /
                          gammaPulseMetrics(high).first;
  if (targetRatio < minRatio || targetRatio > maxRatio) {
    throw std::runtime_error("Requested FWHM/rise ratio cannot be represented by gamma pulse");
  }
  for (int i = 0; i < 100; ++i) {
    const double middle = 0.5 * (low + high);
    const auto [rise, width] = gammaPulseMetrics(middle);
    if (width / rise > targetRatio) low = middle;
    else high = middle;
  }
  p.gammaShape = 0.5 * (low + high);
  const auto [dimensionlessRise, dimensionlessWidth] =
      gammaPulseMetrics(p.gammaShape);
  p.gammaTauNs = p.responseRiseNs / dimensionlessRise;
}

double crossing(const std::vector<double>& time,
                const std::vector<double>& voltage,
                size_t peakIndex,
                double level) {
  for (size_t i = 1; i <= peakIndex; ++i) {
    if (voltage[i - 1] < level && voltage[i] >= level) {
      const double dv = voltage[i] - voltage[i - 1];
      if (dv == 0.0) return time[i];
      return time[i - 1] + (level - voltage[i - 1]) *
                               (time[i] - time[i - 1]) / dv;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

Pulse buildPulse(const std::vector<float>& edgeLow,
                 const std::vector<float>& edgeHigh,
                 const std::vector<int>& counts,
                 const Parameters& p,
                 bool retainSamples) {
  if (edgeLow.size() != edgeHigh.size() || edgeLow.size() != counts.size()) {
    throw std::runtime_error("Time-bin edge/count vector sizes do not match");
  }

  Pulse result;
  if (counts.empty()) return result;

  double firstHit = std::numeric_limits<double>::infinity();
  double lastHit = -std::numeric_limits<double>::infinity();
  double weightedTime = 0.0;
  double weightedTime2 = 0.0;
  for (size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] <= 0) continue;
    const double hitTime = 0.5 * (edgeLow[i] + edgeHigh[i]);
    firstHit = std::min(firstHit, hitTime);
    lastHit = std::max(lastHit, hitTime);
    result.photonCount += counts[i];
    weightedTime += counts[i] * hitTime;
    weightedTime2 += counts[i] * hitTime * hitTime;
  }
  if (!std::isfinite(firstHit)) return result;
  result.arrivalFirst = firstHit;
  result.arrivalLast = lastHit;
  result.arrivalMean = weightedTime / result.photonCount;
  result.arrivalRms = std::sqrt(std::max(
      0.0, weightedTime2 / result.photonCount -
               result.arrivalMean * result.arrivalMean));
  auto weightedQuantile = [&](double fraction) {
    const double target = fraction * result.photonCount;
    double cumulative = 0.0;
    for (size_t i = 0; i < counts.size(); ++i) {
      cumulative += std::max(0, counts[i]);
      if (cumulative >= target) return 0.5 * (edgeLow[i] + edgeHigh[i]);
    }
    return lastHit;
  };
  result.arrivalQ10 = weightedQuantile(0.10);
  result.arrivalQ50 = weightedQuantile(0.50);
  result.arrivalQ90 = weightedQuantile(0.90);

  const double start = std::floor((firstHit - 2.0 * p.sampleStepNs) /
                                  p.sampleStepNs) * p.sampleStepNs;
  const double stop = lastHit + 12.0 * p.responseFwhmNs;
  const size_t nSamples = static_cast<size_t>(std::ceil((stop - start) /
                                                        p.sampleStepNs)) + 1;
  std::vector<double> impulse(nSamples, 0.0);
  for (size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] <= 0) continue;
    const double hitTime = 0.5 * (edgeLow[i] + edgeHigh[i]);
    const long index = std::lround((hitTime - start) / p.sampleStepNs);
    if (index >= 0 && static_cast<size_t>(index) < nSamples) {
      impulse[index] += counts[i];
    }
  }

  std::vector<double> time(nSamples);
  for (size_t i = 0; i < nSamples; ++i) {
    time[i] = start + i * p.sampleStepNs;
  }

  size_t photonPeakIndex = 0;
  for (size_t i = 1; i < impulse.size(); ++i) {
    if (impulse[i] > impulse[photonPeakIndex]) photonPeakIndex = i;
  }
  result.photonT30 = crossing(time, impulse, photonPeakIndex,
                              p.fraction * impulse[photonPeakIndex]);

  const size_t kernelSamples = static_cast<size_t>(std::ceil(
      12.0 * p.responseFwhmNs / p.sampleStepNs)) + 1;
  std::vector<double> kernel(kernelSamples, 0.0);
  for (size_t k = 1; k < kernelSamples; ++k) {
    const double x = (k * p.sampleStepNs) / p.gammaTauNs;
    kernel[k] = std::pow(x / p.gammaShape, p.gammaShape) *
                std::exp(p.gammaShape - x);
  }

  std::vector<double> voltage(nSamples, 0.0);
  for (size_t i = 0; i < impulse.size(); ++i) {
    if (impulse[i] == 0.0) continue;
    const size_t available = std::min(kernel.size(), voltage.size() - i);
    for (size_t k = 1; k < available; ++k) {
      voltage[i + k] += impulse[i] * kernel[k];
    }
  }

  size_t peakIndex = 0;
  for (size_t i = 0; i < nSamples; ++i) {
    if (voltage[i] > result.peak) {
      result.peak = voltage[i];
      peakIndex = i;
    }
  }

  result.t10 = crossing(time, voltage, peakIndex, 0.10 * result.peak);
  result.t30 = crossing(time, voltage, peakIndex, p.fraction * result.peak);
  result.t90 = crossing(time, voltage, peakIndex, 0.90 * result.peak);

  if (retainSamples) {
    result.time = std::move(time);
    result.voltage = std::move(voltage);
  }
  return result;
}

Dataset readDataset(const std::string& label,
                    const std::string& path,
                    const Parameters& p) {
  Dataset dataset{label, path, {}, {}};
  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) throw std::runtime_error("Cannot open " + path);

  // Deliberately omit an explicit ROOT key cycle. The final in-place tree header
  // contains all 3000 events, while explicit autosave cycles are incomplete.
  auto* tree = dynamic_cast<TTree*>(file.Get("CBDsim"));
  if (!tree) throw std::runtime_error("CBDsim tree missing in " + path);

  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<float>> low0(reader, "timeMergedEdgeLowTrig0");
  TTreeReaderValue<std::vector<float>> high0(reader, "timeMergedEdgeHighTrig0");
  TTreeReaderValue<std::vector<int>> count0(reader, "timeMergedCountsTrig0");
  TTreeReaderValue<std::vector<float>> low1(reader, "timeMergedEdgeLowTrig1");
  TTreeReaderValue<std::vector<float>> high1(reader, "timeMergedEdgeHighTrig1");
  TTreeReaderValue<std::vector<int>> count1(reader, "timeMergedCountsTrig1");
  TTreeReaderArray<Float_t> energyDeposits(reader, "Edeps.Edep");
  TTreeReaderArray<Int_t> generatedPhotons(
      reader, "opticalPhotons.opticalPhotonNumber");

  dataset.events.reserve(tree->GetEntries());
  int eventNumber = 0;
  while (reader.Next()) {
    EventResult event;
    event.event = eventNumber;
    for (Float_t value : energyDeposits) event.energyDeposit += value;
    for (Int_t value : generatedPhotons) event.generatedOpticalPhotons += value;
    const bool retain = eventNumber == 0;
    event.trigger[0] = buildPulse(*low0, *high0, *count0, p, retain);
    event.trigger[1] = buildPulse(*low1, *high1, *count1, p, retain);
    if (retain) {
      dataset.example[0] = event.trigger[0];
      dataset.example[1] = event.trigger[1];
    }
    dataset.events.push_back(std::move(event));
    ++eventNumber;
  }
  std::cout << label << ": read " << dataset.events.size() << " events\n";
  return dataset;
}

std::vector<double> select(const Dataset& d, int field, int trigger = 0) {
  std::vector<double> values;
  values.reserve(d.events.size());
  for (const auto& event : d.events) {
    double value = std::numeric_limits<double>::quiet_NaN();
    if (field == 0) value = event.trigger[trigger].t30;
    if (field == 1) value = event.trigger[0].t30 - event.trigger[1].t30;
    if (field == 2) value = event.trigger[trigger].photonCount;
    if (field == 3) value = event.trigger[trigger].peak;
    if (field == 4) value = event.trigger[trigger].t90 - event.trigger[trigger].t10;
    if (field == 5) value = event.trigger[trigger].arrivalMean;
    if (field == 6) value = event.trigger[trigger].arrivalRms;
    if (field == 7) value = event.trigger[trigger].arrivalFirst;
    if (field == 8) value = event.trigger[trigger].arrivalQ90 - event.trigger[trigger].arrivalQ10;
    if (field == 9) value = event.energyDeposit;
    if (field == 10) value = event.generatedOpticalPhotons;
    if (field == 11) value = event.trigger[trigger].photonT30;
    if (field == 12) value = event.trigger[0].photonT30 - event.trigger[1].photonT30;
    if (std::isfinite(value)) values.push_back(value);
  }
  return values;
}

double quantile(const std::vector<double>& sorted, double q) {
  if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
  const double position = q * (sorted.size() - 1);
  const size_t lower = static_cast<size_t>(std::floor(position));
  const size_t upper = static_cast<size_t>(std::ceil(position));
  const double weight = position - lower;
  return sorted[lower] * (1.0 - weight) + sorted[upper] * weight;
}

Summary summarize(std::vector<double> values) {
  Summary s;
  s.n = values.size();
  if (values.empty()) return s;
  s.mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
  double sumSquares = 0.0;
  for (double value : values) sumSquares += (value - s.mean) * (value - s.mean);
  s.sigma = values.size() > 1 ? std::sqrt(sumSquares / (values.size() - 1)) : 0.0;
  std::sort(values.begin(), values.end());
  s.minimum = values.front();
  s.maximum = values.back();
  s.median = quantile(values, 0.5);
  s.half68 = 0.5 * (quantile(values, 0.84) - quantile(values, 0.16));
  return s;
}

TH1D* makeHistogram(const std::string& name,
                    const std::string& title,
                    const std::vector<double>& a,
                    const std::vector<double>& b,
                    int color) {
  std::vector<double> all = a;
  all.insert(all.end(), b.begin(), b.end());
  std::sort(all.begin(), all.end());
  double low = quantile(all, 0.005);
  double high = quantile(all, 0.995);
  if (!(high > low)) { low -= 0.5; high += 0.5; }
  const double margin = 0.08 * (high - low);
  auto* histogram = new TH1D(name.c_str(), title.c_str(), 100,
                             low - margin, high + margin);
  histogram->SetDirectory(nullptr);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(2);
  for (double value : a) histogram->Fill(value);
  if (histogram->Integral() > 0.0) histogram->Scale(1.0 / histogram->Integral());
  return histogram;
}

void overlayPlot(const std::string& path,
                 const std::string& title,
                 const std::string& xTitle,
                 const std::vector<double>& lg,
                 const std::vector<double>& noLg) {
  auto* hLg = makeHistogram("h_lg", title, lg, noLg, kBlue + 1);
  auto* hNoLg = makeHistogram("h_nolg", title, noLg, lg, kRed + 1);
  hLg->GetXaxis()->SetTitle(xTitle.c_str());
  hLg->GetYaxis()->SetTitle("Normalized events");
  hLg->SetMaximum(1.15 * std::max(hLg->GetMaximum(), hNoLg->GetMaximum()));
  TCanvas canvas("canvas", "canvas", 900, 650);
  hLg->Draw("HIST");
  hNoLg->Draw("HIST SAME");
  TLegend legend(0.70, 0.76, 0.88, 0.88);
  legend.AddEntry(hLg, "Light guide", "l");
  legend.AddEntry(hNoLg, "No light guide", "l");
  legend.Draw();
  canvas.SaveAs(path.c_str());
  delete hLg;
  delete hNoLg;
}

void waveformPlot(const std::string& path,
                  const Dataset& lg,
                  const Dataset& noLg,
                  int trigger,
                  double fraction) {
  const Pulse& a = lg.example[trigger];
  const Pulse& b = noLg.example[trigger];
  TGraph gA(a.time.size(), a.time.data(), a.voltage.data());
  TGraph gB(b.time.size(), b.time.data(), b.voltage.data());
  gA.SetLineColor(kBlue + 1); gA.SetLineWidth(2);
  gB.SetLineColor(kRed + 1); gB.SetLineWidth(2);
  TCanvas canvas("waveform", "waveform", 900, 650);
  gA.SetTitle(("Event 0 reconstructed PMT waveform, trigger " +
               std::to_string(trigger) + ";Time [ns];Response [a.u.]").c_str());
  gA.Draw("AL");
  gB.Draw("L SAME");
  TLegend legend(0.66, 0.74, 0.88, 0.88);
  legend.AddEntry(&gA, ("LG, t30=" + std::to_string(a.t30) + " ns").c_str(), "l");
  legend.AddEntry(&gB, ("noLG, t30=" + std::to_string(b.t30) + " ns").c_str(), "l");
  legend.Draw();
  canvas.SaveAs(path.c_str());
}

void responseTemplatePlot(const std::string& path, const Parameters& p) {
  const double step = 0.01;
  const size_t n = static_cast<size_t>(std::ceil(12.0 / step)) + 1;
  std::vector<double> time(n), response(n);
  for (size_t i = 0; i < n; ++i) {
    time[i] = i * step;
    const double x = time[i] / p.gammaTauNs;
    response[i] = x > 0.0
        ? std::pow(x / p.gammaShape, p.gammaShape) * std::exp(p.gammaShape - x)
        : 0.0;
  }
  TGraph graph(n, time.data(), response.data());
  graph.SetLineColor(kBlue + 1);
  graph.SetLineWidth(3);
  graph.SetTitle("Calibrated R2076 provisional single-photon response;Time [ns];Normalized response");
  TCanvas canvas("response_template", "response_template", 900, 650);
  graph.Draw("AL");
  canvas.SaveAs(path.c_str());
}

void pdfStyleTimingPlot(const std::string& path,
                        const Dataset& lg,
                        const Dataset& noLg) {
  gStyle->SetOptStat(0);
  TCanvas canvas("pdf_style_timing", "pdf_style_timing", 1800, 1050);
  canvas.Divide(3, 2, 0.012, 0.025);

  struct Panel {
    const Dataset* dataset;
    int field;
    int trigger;
    const char* title;
    const char* xTitle;
  };
  const std::vector<Panel> panels = {
      {&lg, 0, 0, "with Light Guide: T1", "t_{30}(T1) [ns]"},
      {&lg, 0, 1, "with Light Guide: T2", "t_{30}(T2) [ns]"},
      {&lg, 1, 0, "with Light Guide: T1 - T2", "#Deltat_{30} [ns]"},
      {&noLg, 0, 0, "without Light Guide: T1", "t_{30}(T1) [ns]"},
      {&noLg, 0, 1, "without Light Guide: T2", "t_{30}(T2) [ns]"},
      {&noLg, 1, 0, "without Light Guide: T1 - T2", "#Deltat_{30} [ns]"}};

  std::vector<TH1D*> histograms;
  std::vector<TF1*> fits;
  histograms.reserve(panels.size());
  fits.reserve(panels.size());
  for (size_t index = 0; index < panels.size(); ++index) {
    const Panel& panel = panels[index];
    std::vector<double> values = select(*panel.dataset, panel.field, panel.trigger);
    Summary summary = summarize(values);
    double low = summary.mean - 5.0 * summary.sigma;
    double high = summary.mean + 5.0 * summary.sigma;
    auto* histogram = new TH1D(("pdf_hist_" + std::to_string(index)).c_str(),
                               panel.title, 80, low, high);
    histogram->SetDirectory(nullptr);
    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth(2);
    histogram->GetXaxis()->SetTitle(panel.xTitle);
    histogram->GetYaxis()->SetTitle("Events");
    for (double value : values) histogram->Fill(value);

    auto* fit = new TF1(("pdf_fit_" + std::to_string(index)).c_str(), "gaus",
                        summary.mean - 2.5 * summary.sigma,
                        summary.mean + 2.5 * summary.sigma);
    fit->SetParameters(histogram->GetMaximum(), summary.mean, summary.sigma);
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(3);
    histogram->Fit(fit, "RQ0");

    canvas.cd(index + 1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);
    histogram->Draw("HIST");
    fit->Draw("SAME");
    TLatex label;
    label.SetNDC();
    label.SetTextColor(kRed + 1);
    label.SetTextSize(0.050);
    label.DrawLatex(0.16, 0.86,
                    Form("#sigma = %.1f #pm %.1f ps",
                         1000.0 * std::abs(fit->GetParameter(2)),
                         1000.0 * fit->GetParError(2)));
    label.SetTextColor(kBlack);
    label.SetTextSize(0.040);
    label.DrawLatex(0.16, 0.79,
                    Form("#mu = %.4f ns", fit->GetParameter(1)));
    histograms.push_back(histogram);
    fits.push_back(fit);
  }
  canvas.SaveAs(path.c_str());
  for (auto* histogram : histograms) delete histogram;
  for (auto* fit : fits) delete fit;
}

struct GaussianResult {
  double mean = 0.0;
  double sigma = 0.0;
  double sigmaError = 0.0;
};

GaussianResult gaussianResult(const std::vector<double>& values,
                              const std::string& name) {
  Summary summary = summarize(values);
  TH1D histogram((name + "_hist").c_str(), "", 80,
                 summary.mean - 5.0 * summary.sigma,
                 summary.mean + 5.0 * summary.sigma);
  for (double value : values) histogram.Fill(value);
  TF1 fit((name + "_fit").c_str(), "gaus",
          summary.mean - 2.5 * summary.sigma,
          summary.mean + 2.5 * summary.sigma);
  fit.SetParameters(histogram.GetMaximum(), summary.mean, summary.sigma);
  histogram.Fit(&fit, "RQ0");
  return {fit.GetParameter(1), std::abs(fit.GetParameter(2)),
          fit.GetParError(2)};
}

void writePdfStyleSummary(const fs::path& path,
                          const Dataset& lg,
                          const Dataset& noLg) {
  std::ofstream out(path);
  out << std::fixed << std::setprecision(2);
  out << "# PDF-style Gaussian-fit timing result\n\n"
      << "| Configuration | T1 sigma [ps] | T2 sigma [ps] | T1-T2 sigma [ps] |\n"
      << "|---|---:|---:|---:|\n";
  GaussianResult lg0, lg1, lgd, no0, no1, nod;
  int index = 0;
  for (const Dataset* dataset : {&lg, &noLg}) {
    GaussianResult t1 = gaussianResult(select(*dataset, 0, 0),
                                       "summary_" + std::to_string(index++));
    GaussianResult t2 = gaussianResult(select(*dataset, 0, 1),
                                       "summary_" + std::to_string(index++));
    GaussianResult delta = gaussianResult(select(*dataset, 1),
                                          "summary_" + std::to_string(index++));
    out << '|' << dataset->label << '|'
        << 1000.0 * t1.sigma << " +/- " << 1000.0 * t1.sigmaError << '|'
        << 1000.0 * t2.sigma << " +/- " << 1000.0 * t2.sigmaError << '|'
        << 1000.0 * delta.sigma << " +/- " << 1000.0 * delta.sigmaError << "|\n";
    if (dataset == &lg) { lg0 = t1; lg1 = t2; lgd = delta; }
    else { no0 = t1; no1 = t2; nod = delta; }
  }
  out << "\n- Light-guide reduction in Gaussian-fit T1-T2 sigma: "
      << 100.0 * (1.0 - lgd.sigma / nod.sigma) << "%\n"
      << "- Primary presentation figure: `PDF_STYLE_TIMING_RESULT.png`\n"
      << "- This summary uses the Gaussian-fit sigma shown in the figure; `WAVEFORM_COMPARISON.md` also reports the unbinned sample RMS.\n";
}

void writeCsv(const fs::path& path, const Dataset& lg, const Dataset& noLg) {
  std::ofstream out(path);
  out << "configuration,event,t30_trig0_ns,t30_trig1_ns,delta_t_ns,"
         "photons_trig0,photons_trig1,peak_trig0_au,peak_trig1_au,"
         "rise_10_90_trig0_ns,rise_10_90_trig1_ns,"
         "arrival_mean_trig0_ns,arrival_mean_trig1_ns,"
         "arrival_rms_trig0_ns,arrival_rms_trig1_ns,"
         "arrival_first_trig0_ns,arrival_first_trig1_ns,"
         "arrival_last_trig0_ns,arrival_last_trig1_ns,"
         "arrival_q10_trig0_ns,arrival_q10_trig1_ns,"
         "arrival_q50_trig0_ns,arrival_q50_trig1_ns,"
         "arrival_q90_trig0_ns,arrival_q90_trig1_ns,"
         "photon_cfd30_trig0_ns,photon_cfd30_trig1_ns,photon_delta_t_ns,"
         "energy_deposit,generated_optical_photons\n";
  out << std::setprecision(10);
  for (const Dataset* dataset : {&lg, &noLg}) {
    for (const auto& event : dataset->events) {
      out << dataset->label << ',' << event.event << ','
          << event.trigger[0].t30 << ',' << event.trigger[1].t30 << ','
          << event.trigger[0].t30 - event.trigger[1].t30 << ','
          << event.trigger[0].photonCount << ',' << event.trigger[1].photonCount << ','
          << event.trigger[0].peak << ',' << event.trigger[1].peak << ','
          << event.trigger[0].t90 - event.trigger[0].t10 << ','
          << event.trigger[1].t90 - event.trigger[1].t10 << ','
          << event.trigger[0].arrivalMean << ',' << event.trigger[1].arrivalMean << ','
          << event.trigger[0].arrivalRms << ',' << event.trigger[1].arrivalRms << ','
          << event.trigger[0].arrivalFirst << ',' << event.trigger[1].arrivalFirst << ','
          << event.trigger[0].arrivalLast << ',' << event.trigger[1].arrivalLast << ','
          << event.trigger[0].arrivalQ10 << ',' << event.trigger[1].arrivalQ10 << ','
          << event.trigger[0].arrivalQ50 << ',' << event.trigger[1].arrivalQ50 << ','
          << event.trigger[0].arrivalQ90 << ',' << event.trigger[1].arrivalQ90 << ','
          << event.trigger[0].photonT30 << ',' << event.trigger[1].photonT30 << ','
          << event.trigger[0].photonT30 - event.trigger[1].photonT30 << ','
          << event.energyDeposit << ',' << event.generatedOpticalPhotons << '\n';
    }
  }
}

void writeRoot(const fs::path& path,
               const Dataset& lg,
               const Dataset& noLg,
               const Parameters& p) {
  TFile output(path.c_str(), "RECREATE");
  TParameter<double>("response_rise_10_90_ns", p.responseRiseNs).Write();
  TParameter<double>("response_fwhm_ns", p.responseFwhmNs).Write();
  TParameter<double>("gamma_shape", p.gammaShape).Write();
  TParameter<double>("gamma_tau_ns", p.gammaTauNs).Write();
  TParameter<double>("cfd_fraction", p.fraction).Write();
  TParameter<double>("sample_step_ns", p.sampleStepNs).Write();

  int configuration = 0;
  int eventNumber = 0;
  double t30Trig0 = 0.0, t30Trig1 = 0.0, deltaT = 0.0;
  double photonsTrig0 = 0.0, photonsTrig1 = 0.0;
  double peakTrig0 = 0.0, peakTrig1 = 0.0;
  double riseTrig0 = 0.0, riseTrig1 = 0.0;
  double arrivalMean0 = 0.0, arrivalMean1 = 0.0;
  double arrivalRms0 = 0.0, arrivalRms1 = 0.0;
  double arrivalFirst0 = 0.0, arrivalFirst1 = 0.0;
  double arrivalLast0 = 0.0, arrivalLast1 = 0.0;
  double arrivalQ10_0 = 0.0, arrivalQ10_1 = 0.0;
  double arrivalQ50_0 = 0.0, arrivalQ50_1 = 0.0;
  double arrivalQ90_0 = 0.0, arrivalQ90_1 = 0.0;
  double energyDeposit = 0.0, generatedOpticalPhotons = 0.0;
  double photonT30_0 = 0.0, photonT30_1 = 0.0, photonDeltaT = 0.0;
  TTree metrics("event_metrics", "Reconstructed waveform metrics; configuration 1=LG, 0=noLG");
  metrics.Branch("configuration", &configuration);
  metrics.Branch("event", &eventNumber);
  metrics.Branch("t30_trig0_ns", &t30Trig0);
  metrics.Branch("t30_trig1_ns", &t30Trig1);
  metrics.Branch("delta_t_ns", &deltaT);
  metrics.Branch("photons_trig0", &photonsTrig0);
  metrics.Branch("photons_trig1", &photonsTrig1);
  metrics.Branch("peak_trig0_au", &peakTrig0);
  metrics.Branch("peak_trig1_au", &peakTrig1);
  metrics.Branch("rise_10_90_trig0_ns", &riseTrig0);
  metrics.Branch("rise_10_90_trig1_ns", &riseTrig1);
  metrics.Branch("arrival_mean_trig0_ns", &arrivalMean0);
  metrics.Branch("arrival_mean_trig1_ns", &arrivalMean1);
  metrics.Branch("arrival_rms_trig0_ns", &arrivalRms0);
  metrics.Branch("arrival_rms_trig1_ns", &arrivalRms1);
  metrics.Branch("arrival_first_trig0_ns", &arrivalFirst0);
  metrics.Branch("arrival_first_trig1_ns", &arrivalFirst1);
  metrics.Branch("arrival_last_trig0_ns", &arrivalLast0);
  metrics.Branch("arrival_last_trig1_ns", &arrivalLast1);
  metrics.Branch("arrival_q10_trig0_ns", &arrivalQ10_0);
  metrics.Branch("arrival_q10_trig1_ns", &arrivalQ10_1);
  metrics.Branch("arrival_q50_trig0_ns", &arrivalQ50_0);
  metrics.Branch("arrival_q50_trig1_ns", &arrivalQ50_1);
  metrics.Branch("arrival_q90_trig0_ns", &arrivalQ90_0);
  metrics.Branch("arrival_q90_trig1_ns", &arrivalQ90_1);
  metrics.Branch("energy_deposit", &energyDeposit);
  metrics.Branch("generated_optical_photons", &generatedOpticalPhotons);
  metrics.Branch("photon_cfd30_trig0_ns", &photonT30_0);
  metrics.Branch("photon_cfd30_trig1_ns", &photonT30_1);
  metrics.Branch("photon_delta_t_ns", &photonDeltaT);
  for (const Dataset* dataset : {&lg, &noLg}) {
    configuration = dataset == &lg ? 1 : 0;
    for (const auto& event : dataset->events) {
      eventNumber = event.event;
      t30Trig0 = event.trigger[0].t30;
      t30Trig1 = event.trigger[1].t30;
      deltaT = t30Trig0 - t30Trig1;
      photonsTrig0 = event.trigger[0].photonCount;
      photonsTrig1 = event.trigger[1].photonCount;
      peakTrig0 = event.trigger[0].peak;
      peakTrig1 = event.trigger[1].peak;
      riseTrig0 = event.trigger[0].t90 - event.trigger[0].t10;
      riseTrig1 = event.trigger[1].t90 - event.trigger[1].t10;
      arrivalMean0 = event.trigger[0].arrivalMean;
      arrivalMean1 = event.trigger[1].arrivalMean;
      arrivalRms0 = event.trigger[0].arrivalRms;
      arrivalRms1 = event.trigger[1].arrivalRms;
      arrivalFirst0 = event.trigger[0].arrivalFirst;
      arrivalFirst1 = event.trigger[1].arrivalFirst;
      arrivalLast0 = event.trigger[0].arrivalLast;
      arrivalLast1 = event.trigger[1].arrivalLast;
      arrivalQ10_0 = event.trigger[0].arrivalQ10;
      arrivalQ10_1 = event.trigger[1].arrivalQ10;
      arrivalQ50_0 = event.trigger[0].arrivalQ50;
      arrivalQ50_1 = event.trigger[1].arrivalQ50;
      arrivalQ90_0 = event.trigger[0].arrivalQ90;
      arrivalQ90_1 = event.trigger[1].arrivalQ90;
      energyDeposit = event.energyDeposit;
      generatedOpticalPhotons = event.generatedOpticalPhotons;
      photonT30_0 = event.trigger[0].photonT30;
      photonT30_1 = event.trigger[1].photonT30;
      photonDeltaT = photonT30_0 - photonT30_1;
      metrics.Fill();
    }
  }
  metrics.Write();

  for (const Dataset* dataset : {&lg, &noLg}) {
    for (int trigger = 0; trigger < 2; ++trigger) {
      const Pulse& pulse = dataset->example[trigger];
      TGraph graph(pulse.time.size(), pulse.time.data(), pulse.voltage.data());
      graph.SetName(("waveform_event0_" + dataset->label + "_trig" +
                     std::to_string(trigger)).c_str());
      graph.SetTitle("Event 0 reconstructed waveform;Time [ns];Response [a.u.]");
      graph.Write();
    }
  }
  output.Close();
}

void writeSummary(const fs::path& path,
                  const Dataset& lg,
                  const Dataset& noLg,
                  const Parameters& p) {
  std::ofstream out(path);
  out << std::fixed << std::setprecision(6);
  out << "# Light-guide PMT waveform comparison\n\n"
      << "## Analysis configuration\n\n"
      << "- Single-photon response: causal gamma pulse calibrated to the PDF values\n"
      << "- Response 10-90% rise: " << p.responseRiseNs << " ns\n"
      << "- Response FWHM: " << p.responseFwhmNs << " ns\n"
      << "- Gamma shape: " << p.gammaShape << "\n"
      << "- Gamma time scale: " << p.gammaTauNs << " ns\n"
      << "- CFD fraction: " << 100.0 * p.fraction << "% of each event peak\n"
      << "- Waveform sampling: " << p.sampleStepNs << " ns\n"
      << "- Events: " << lg.events.size() << " LG, " << noLg.events.size() << " noLG\n\n"
      << "> Every stored photon-arrival bin and its count are used in the waveform convolution. The ROOT `waveForm` maps contain zero samples in all events, so an ADC-baseline RMS cannot be extracted. PMT transit-time spread, electronics noise, gain fluctuations, and digitizer jitter are not stored and are not added.\n\n";

  out << "## Timing results\n\n"
      << "| Configuration | Quantity | Mean [ns] | RMS sigma [ps] | central 68% half-width [ps] |\n"
      << "|---|---|---:|---:|---:|\n";
  for (const Dataset* dataset : {&lg, &noLg}) {
    for (int trigger = 0; trigger < 2; ++trigger) {
      Summary s = summarize(select(*dataset, 0, trigger));
      out << '|' << dataset->label << "|t30 trigger " << trigger << '|'
          << s.mean << '|' << 1000.0 * s.sigma << '|'
          << 1000.0 * s.half68 << "|\n";
    }
    Summary delta = summarize(select(*dataset, 1));
    out << '|' << dataset->label << "|t30 trig0 - trig1|"
        << delta.mean << '|' << 1000.0 * delta.sigma << '|'
        << 1000.0 * delta.half68 << "|\n";
    out << '|' << dataset->label << "|single-PMT estimate sigma(delta)/sqrt(2)|- |"
        << 1000.0 * delta.sigma / std::sqrt(2.0) << "|- |\n";
  }

  out << "\n## Photon-histogram CFD validation\n\n"
      << "This applies the same 30% leading-edge interpolation directly to the 200 ps photon histogram before PMT convolution. It is retained as an implementation cross-check against the PDF's approximately 60 ps (LG) and 80 ps (noLG) simulation.\n\n"
      << "| Configuration | sigma photon t30 trigger 0 [ps] | sigma photon t30 trigger 1 [ps] | sigma photon delta t [ps] |\n"
      << "|---|---:|---:|---:|\n";
  for (const Dataset* dataset : {&lg, &noLg}) {
    Summary p0 = summarize(select(*dataset, 11, 0));
    Summary p1 = summarize(select(*dataset, 11, 1));
    Summary pd = summarize(select(*dataset, 12));
    out << '|' << dataset->label << '|' << 1000.0 * p0.sigma << '|'
        << 1000.0 * p1.sigma << '|' << 1000.0 * pd.sigma << "|\n";
  }

  out << "\n## ROOT-extracted optical timing\n\n"
      << "The arrival RMS below is the within-event RMS of detected-photon arrival times; it is already encoded in each reconstructed waveform and is not added again as Gaussian noise.\n\n"
      << "| Configuration | Trigger | Mean arrival [ns] | Mean within-event arrival RMS [ns] | Mean first photon [ns] | Event sigma(first) [ps] | Mean q10-q90 [ns] |\n"
      << "|---|---:|---:|---:|---:|---:|---:|\n";
  for (const Dataset* dataset : {&lg, &noLg}) {
    for (int trigger = 0; trigger < 2; ++trigger) {
      Summary arrivalMean = summarize(select(*dataset, 5, trigger));
      Summary arrivalRms = summarize(select(*dataset, 6, trigger));
      Summary first = summarize(select(*dataset, 7, trigger));
      Summary central = summarize(select(*dataset, 8, trigger));
      out << '|' << dataset->label << '|' << trigger << '|'
          << arrivalMean.mean << '|' << arrivalRms.mean << '|'
          << first.mean << '|' << 1000.0 * first.sigma << '|'
          << central.mean << "|\n";
    }
  }

  out << "\n## ROOT-extracted event production\n\n"
      << "| Configuration | Mean energy deposit | RMS energy deposit | Mean generated optical photons | RMS generated optical photons |\n"
      << "|---|---:|---:|---:|---:|\n";
  for (const Dataset* dataset : {&lg, &noLg}) {
    Summary edep = summarize(select(*dataset, 9));
    Summary generated = summarize(select(*dataset, 10));
    out << '|' << dataset->label << '|' << edep.mean << '|' << edep.sigma
        << '|' << generated.mean << '|' << generated.sigma << "|\n";
  }

  out << "\n## Photon and pulse-shape results\n\n"
      << "| Configuration | Trigger | Mean photons | Mean peak [a.u.] | Mean 10-90% rise [ns] |\n"
      << "|---|---:|---:|---:|---:|\n";
  for (const Dataset* dataset : {&lg, &noLg}) {
    for (int trigger = 0; trigger < 2; ++trigger) {
      Summary photons = summarize(select(*dataset, 2, trigger));
      Summary peak = summarize(select(*dataset, 3, trigger));
      Summary rise = summarize(select(*dataset, 4, trigger));
      out << '|' << dataset->label << '|' << trigger << '|'
          << photons.mean << '|' << peak.mean << '|' << rise.mean << "|\n";
    }
  }

  Summary lgDelta = summarize(select(lg, 1));
  Summary noLgDelta = summarize(select(noLg, 1));
  out << "\n## Direct comparison\n\n"
      << "- LG sigma(delta t): " << 1000.0 * lgDelta.sigma << " ps\n"
      << "- noLG sigma(delta t): " << 1000.0 * noLgDelta.sigma << " ps\n"
      << "- Ratio LG/noLG: " << lgDelta.sigma / noLgDelta.sigma << "\n";
}

int main(int argc, char** argv) {
  if (argc < 4 || argc > 8) {
    std::cerr << "Usage: " << argv[0]
              << " LG.root noLG.root output_dir [responseRise_ns] [responseFwhm_ns]"
                 " [CFD_fraction] [sample_step_ns]\n";
    return 2;
  }
  try {
    Parameters p;
    if (argc > 4) p.responseRiseNs = std::stod(argv[4]);
    if (argc > 5) p.responseFwhmNs = std::stod(argv[5]);
    if (argc > 6) p.fraction = std::stod(argv[6]);
    if (argc > 7) p.sampleStepNs = std::stod(argv[7]);
    if (!(p.responseRiseNs > 0.0 && p.responseFwhmNs > 0.0 &&
          p.fraction > 0.0 && p.fraction < 1.0 && p.sampleStepNs > 0.0)) {
      throw std::runtime_error("Invalid response or CFD parameters");
    }
    calibrateGammaPulse(p);

    const fs::path outputDir = argv[3];
    fs::create_directories(outputDir);
    gStyle->SetOptStat(0);
    Dataset lg = readDataset("LG", argv[1], p);
    Dataset noLg = readDataset("noLG", argv[2], p);
    if (lg.events.size() != 3000 || noLg.events.size() != 3000) {
      throw std::runtime_error("Expected 3000 events in each final CBDsim tree");
    }

    writeCsv(outputDir / "event_metrics.csv", lg, noLg);
    writeRoot(outputDir / "waveform_analysis.root", lg, noLg, p);
    writeSummary(outputDir / "WAVEFORM_COMPARISON.md", lg, noLg, p);
    writePdfStyleSummary(outputDir / "PDF_STYLE_RESULT_SUMMARY.md", lg, noLg);
    overlayPlot((outputDir / "t30_trigger0.png").string(),
                "30% CFD time, trigger 0", "t_{30} [ns]",
                select(lg, 0, 0), select(noLg, 0, 0));
    overlayPlot((outputDir / "t30_trigger1.png").string(),
                "30% CFD time, trigger 1", "t_{30} [ns]",
                select(lg, 0, 1), select(noLg, 0, 1));
    overlayPlot((outputDir / "delta_t30.png").string(),
                "PMT time difference at 30% CFD", "t_{30,0}-t_{30,1} [ns]",
                select(lg, 1), select(noLg, 1));
    overlayPlot((outputDir / "photon_delta_t30.png").string(),
                "Photon-histogram time difference at 30% CFD",
                "t_{30,0}-t_{30,1} [ns]", select(lg, 12), select(noLg, 12));
    overlayPlot((outputDir / "photon_count_trigger0.png").string(),
                "Detected photons, trigger 0", "Detected photons",
                select(lg, 2, 0), select(noLg, 2, 0));
    waveformPlot((outputDir / "waveform_event0_trigger0.png").string(),
                 lg, noLg, 0, p.fraction);
    waveformPlot((outputDir / "waveform_event0_trigger1.png").string(),
                 lg, noLg, 1, p.fraction);
    responseTemplatePlot((outputDir / "r2076_response_template.png").string(), p);
    pdfStyleTimingPlot((outputDir / "PDF_STYLE_TIMING_RESULT.png").string(),
                       lg, noLg);
    std::cout << "Results written to " << outputDir << '\n';
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
