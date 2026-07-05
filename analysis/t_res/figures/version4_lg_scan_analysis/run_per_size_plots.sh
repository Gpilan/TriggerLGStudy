#!/usr/bin/env bash
set -eo pipefail
cd /u/user/rmsvlf000/Trigger
source envset.sh

FIG=analysis/t_res/figures/version4_lg_scan_analysis
DATA=analysis/t_res/data

tags=(1x1 1p5x1p5 2x2 2p5x2p5 3x3 3p5x3p5)
for tag in "${tags[@]}"; do
  stem="v4_60GeV_e-_LG_${tag}_3000_0"
  root="${DATA}/${stem}.root"
  cfd_out="${FIG}/cfd_locallinear_0p30_to5ps/cfd0p30_waveform_recon_per_trigger_${stem}_1x3_to5ps.png"
  mpv_out="${FIG}/mpv_matrix_200ps/event_optical_mpv_matrix_${stem}_200ps_waveform.png"

  if [[ ! -f "${cfd_out}" ]]; then
    echo "=== CFD: ${tag} ==="
    python analysis/plot_cfd_per_trigger_locallinear_peak.py "${root}" \
      --peak-frac 0.3 --plot-bin-ps 5 \
      -o "${cfd_out}"
  else
    echo "=== CFD skip (exists): ${tag} ==="
  fi

  if [[ ! -f "${mpv_out}" ]]; then
    echo "=== MPV: ${tag} ==="
    python analysis/plot_event_optical_mpv_matrix_root.py "${root}" \
      --arrival-bin-ps 200 \
      -o "${mpv_out}"
  else
    echo "=== MPV skip (exists): ${tag} ==="
  fi
done

echo "=== done ==="
