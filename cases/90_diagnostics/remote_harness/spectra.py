"""Limit-cycle test on the per-step cavitated_fraction of a finished run.

The question the dt ladder has to answer is whether the T1 oscillation is a
2*dt numerical artifact (amplitude should collapse as dt falls) or a resolved
physical motion (amplitude should converge to a non-zero value).

Naive argmax on the spectrum is not enough: over a 0.03 s window the first
non-zero DFT bin sits at 33.3 Hz, which coincidentally equals half the 66.7 Hz
shaft rate, so residual settling drift masquerades as half-frequency whirl.
This script therefore
  * linearly detrends the window before transforming, so the settling drift
    does not dominate the spectrum,
  * reports the amplitude exactly at the Nyquist bin, which is the 2*dt
    signature, separately from the rest,
  * reports the strongest peak below Nyquist while skipping the lowest bins,
  * reports a direct time-domain measure of step-to-step alternation.

Usage: python spectra.py <run_dir> [<run_dir> ...]   (each holds results/)
"""
import sys
from pathlib import Path

import numpy as np

T_LO, T_HI = 0.03, 0.06
SKIP_BINS = 3          # drop DC + the lowest bins the window cannot resolve
SHAFT_HZ = 418.879020479 / (2.0 * np.pi)


def load(run_dir):
    path = Path(run_dir) / "results" / "diagnostics.csv"
    if not path.exists():
        path = Path(run_dir) / "results_remote" / "diagnostics.csv"
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    head = lines[0].split(",")
    rows = [dict(zip(head, l.split(","))) for l in lines[1:]]
    t = np.array([float(r["time"]) for r in rows])
    c = np.array([float(r["cavitated_fraction"]) for r in rows])
    return t, c


def analyse(run_dir):
    t, c = load(run_dir)
    dt = np.median(np.diff(t))
    sel = (t >= T_LO) & (t <= T_HI)
    tw, cw = t[sel], c[sel]
    n = len(tw)
    if n < 32:
        return None

    # linear detrend: settling drift is not oscillation
    coef = np.polyfit(tw, cw, 1)
    resid = cw - np.polyval(coef, tw)

    win = np.hanning(n)
    amp = np.abs(np.fft.rfft(resid * win)) * 2.0 / win.sum()
    freq = np.fft.rfftfreq(n, d=dt)
    f_nyq = 1.0 / (2.0 * dt)

    body = amp[SKIP_BINS:-1]
    fbody = freq[SKIP_BINS:-1]
    k = int(np.argmax(body))

    # time-domain 2*dt alternation: second difference isolates step-to-step flip
    alt = float(np.mean(np.abs(np.diff(resid, n=2)))) / 4.0

    return {
        "run": Path(run_dir).name,
        "dt": dt,
        "n_samples": n,
        "f_nyquist_hz": f_nyq,
        "df_hz": float(freq[1]),
        "raw_p2p": float(np.ptp(cw)),
        "drift_over_window": float(np.polyval(coef, tw[-1]) - np.polyval(coef, tw[0])),
        "detrended_std": float(np.std(resid)),
        "detrended_p2p": float(np.ptp(resid)),
        "amp_at_nyquist": float(amp[-1]),
        "alternation_amp": alt,
        "peak_below_nyq_hz": float(fbody[k]),
        "peak_below_nyq_amp": float(body[k]),
        "peak_over_shaft": float(fbody[k]) / SHAFT_HZ,
    }


for d in sys.argv[1:]:
    r = analyse(d)
    if r is None:
        print(f"{d}: too few samples in [{T_LO},{T_HI}]")
        continue
    print(f"\n=== {r['run']}  dt={r['dt']:.2e}  n={r['n_samples']} "
          f"(bin width {r['df_hz']:.1f} Hz, Nyquist {r['f_nyquist_hz']:.0f} Hz) ===")
    print(f"  cav_frac raw p2p        {r['raw_p2p']:.6e}")
    print(f"  linear drift over window {r['drift_over_window']:+.6e}")
    print(f"  detrended std / p2p     {r['detrended_std']:.6e} / {r['detrended_p2p']:.6e}")
    print(f"  amplitude AT Nyquist    {r['amp_at_nyquist']:.6e}   <-- 2*dt signature")
    print(f"  time-domain alternation {r['alternation_amp']:.6e}")
    print(f"  strongest peak < Nyq    {r['peak_below_nyq_hz']:.1f} Hz "
          f"(= {r['peak_over_shaft']:.3f} x shaft {SHAFT_HZ:.1f} Hz), "
          f"amp {r['peak_below_nyq_amp']:.6e}")
