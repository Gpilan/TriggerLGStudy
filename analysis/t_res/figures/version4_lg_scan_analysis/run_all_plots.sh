#!/usr/bin/env bash
set -eo pipefail
cd /u/user/rmsvlf000/Trigger
source envset.sh

FIG=analysis/t_res/figures/version4_lg_scan_analysis
DATA=analysis/t_res/data

echo "=== [1/3] trend ==="
python analysis/plot_lg_scan_trend_v4.py \
  -o "${FIG}/trend/trigger_size_trend_mpv_sigma_v4plusv3_from1_to4_xmax5.png" \
  --csv "${FIG}/trend/trigger_size_trend_mpv_sigma_v4plusv3_from1_to4_xmax5.csv"

tags=(1x1 1p5x1p5 2x2 2p5x2p5 3x3 3p5x3p5)
for tag in "${tags[@]}"; do
  stem="v4_60GeV_e-_LG_${tag}_3000_0"
  root="${DATA}/${stem}.root"
  echo "=== CFD: ${tag} ==="
  python analysis/plot_cfd_per_trigger_locallinear_peak.py "${root}" \
    --peak-frac 0.3 --plot-bin-ps 5 \
    -o "${FIG}/cfd_locallinear_0p30_to5ps/cfd0p30_waveform_recon_per_trigger_${stem}_1x3_to5ps.png"
  echo "=== MPV: ${tag} ==="
  python analysis/plot_event_optical_mpv_matrix_root.py "${root}" \
    --arrival-bin-ps 200 \
    -o "${FIG}/mpv_matrix_200ps/event_optical_mpv_matrix_${stem}_200ps_waveform.png"
done

echo "=== done ==="
