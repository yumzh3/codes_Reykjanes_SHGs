#!/usr/bin/env python3
"""
Written by Mingzhen Yu with assitance of Claude Fable5.
Multitaper spectral analysis of the core 22 MgO time series.

Tests whether the quasi-periodic fluctuations in the detrended MgO series of
core GeoB25422-1 ("core 22") represent a statistically significant periodicity,
in particular near the 41 kyr obliquity band. In one run this script:

  1. Loads the unsmoothed 2-kyr binned MgO series, fills unsampled bins by
     linear interpolation, and removes a linear trend.
  2. Computes the multitaper (MTM) power spectrum (time-bandwidth product
     NW = 3, K = 4 DPSS tapers).
  3. Assesses significance against NSUR AR(1) surrogate series matched to the
     lag-1 autocorrelation and variance of the detrended data, using:
       - frequency-local 95% and 99% confidence levels, and
       - a global test statistic (the maximum median-normalized spectral peak
         across all frequencies), which accounts for the search across
         frequencies.
  4. Repeats the analysis for two windows: the full record (0-131 ka) and the
     0-115 ka interval used in the cross-correlation analysis (Section 4).
  5. Prints a results table, writes the spectra to CSV (for replotting in any
     software), and draws the two-panel Figure S7.

Output: figS7_MgO_spectrum.png/.pdf  (the figure)
        figS7_spectra.csv            (window, frequency, period, observed PSD,
                                      AR(1) median / 95% / 99% levels)
        results table to stdout      (PSD vs local levels at reference periods;
                                      global p-value per window)

Reproducibility: fixed seed; Monte Carlo error on the global p-value at
NSUR = 10,000 is ~+/-0.002; raise NSUR for final manuscript values.

Usage: python3 Mgo_spectral_analysis.py [data.xlsx]   (default 'SHG_data_for_xcorr.xlsx';
       expects a sheet 'core 22' with 'Age' and 'MgO' columns)
Dependencies: numpy, scipy, openpyxl, matplotlib
"""
import sys
import csv
import numpy as np
import openpyxl
from scipy.signal.windows import dpss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------- configuration --------------------------------
DATA = sys.argv[1] if len(sys.argv) > 1 else 'SHG_data_for_xcorr.xlsx'
NSUR = 10000
SEED = 7
DT = 2.0                       # bin width, kyr
NW, K = 3, 4                   # multitaper time-bandwidth product, no. of tapers
WINDOWS = [('full record (0-131 ka)', 132), ('restricted (0-115 ka)', 116)]
REF_PERIODS = [100, 54, 41, 29, 23]   # kyr, reported in the results table
OBLIQUITY = 41.0               # kyr, marked in the figure


def load_mgo(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    rows = list(wb['core 22'].iter_rows(values_only=True))
    hdr = rows[0]
    body = [r for r in rows[1:] if isinstance(r[0], (int, float))]
    arr = np.array([[c if isinstance(c, (int, float)) else np.nan for c in r]
                    for r in body], dtype=float)
    return arr[:, 0], arr[:, hdr.index('MgO')]


def interp_nan(x):
    x = x.copy()
    n = np.isnan(x)
    if n.any():
        x[n] = np.interp(np.flatnonzero(n), np.flatnonzero(~n), x[~n])
    return x


def mtm_psd(x, dt, nw=NW, k=K):
    """Multitaper power spectral density (DC bin dropped)."""
    n = len(x)
    tapers = dpss(n, nw, k)
    S = 0
    for i in range(k):
        S = S + np.abs(np.fft.rfft(x * tapers[i])) ** 2
    f = np.fft.rfftfreq(n, dt)
    return f[1:], (S / k)[1:]


def ar1_surrogates(x, n_sur, rng):
    """AR(1) surrogates matched to lag-1 autocorrelation and variance of x."""
    x0 = x - x.mean()
    r1 = np.clip(np.corrcoef(x0[:-1], x0[1:])[0, 1], -0.99, 0.99)
    n = len(x)
    e = rng.normal(size=(n_sur, n))
    y = np.zeros((n_sur, n))
    y[:, 0] = e[:, 0]
    for i in range(1, n):
        y[:, i] = r1 * y[:, i - 1] + np.sqrt(1 - r1 ** 2) * e[:, i]
    return y * x0.std(), r1


def main():
    age, mgo = load_mgo(DATA)
    csv_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300)
    for ax, (label, amax), plab in zip(axes, WINDOWS, ['a', 'b']):
        sel = age <= amax
        m = interp_nan(mgo[sel])
        a = age[sel]
        x = m - np.polyval(np.polyfit(a, m, 1), a)          # linear detrend
        f, S = mtm_psd(x, DT)
        rng = np.random.default_rng(SEED)
        sur, r1 = ar1_surrogates(x, NSUR, rng)
        sur = sur - sur.mean(1, keepdims=True)
        Ssur = np.array([mtm_psd(s, DT)[1] for s in sur])
        med = np.median(Ssur, axis=0)
        loc95 = np.percentile(Ssur, 95, axis=0)
        loc99 = np.percentile(Ssur, 99, axis=0)
        # global test: max median-normalized peak, observed vs surrogates
        gmax_null = np.max(Ssur / med, axis=1)
        obs_norm = S / med
        ipk = int(np.argmax(obs_norm))
        p_glob = (np.sum(gmax_null >= obs_norm[ipk]) + 1) / (NSUR + 1)
        # ---- report ----
        print(f'== {label}: N = {sel.sum()}, AR(1) r1 = {r1:.2f}')
        print(f'   strongest normalized peak: period = {1 / f[ipk]:.1f} kyr, '
              f'global p = {p_glob:.3f}')
        for per in REF_PERIODS:
            i = int(np.argmin(np.abs(f - 1 / per)))
            tag = ('> 99%' if S[i] > loc99[i] else
                   '> 95%' if S[i] > loc95[i] else 'n.s.')
            print(f'   period {1 / f[i]:>6.1f} kyr: PSD = {S[i]:.4f}  '
                  f'local95 = {loc95[i]:.4f}  local99 = {loc99[i]:.4f}  [{tag}]')
        for i in range(len(f)):
            csv_rows.append(dict(window=label, freq_cyc_per_kyr=round(float(f[i]), 5),
                                 period_kyr=round(float(1 / f[i]), 2),
                                 psd_obs=round(float(S[i]), 5),
                                 ar1_median=round(float(med[i]), 5),
                                 ar1_p95=round(float(loc95[i]), 5),
                                 ar1_p99=round(float(loc99[i]), 5)))
        # ---- plot ----
        ax.plot(f, S, 'k-', lw=1.2, label='MgO (detrended)')
        ax.plot(f, med, '-', color='tab:red', lw=0.9, label='AR(1) median')
        ax.plot(f, loc95, '--', color='tab:red', lw=0.9, label='AR(1) 95%')
        ax.plot(f, loc99, ':', color='tab:red', lw=0.9, label='AR(1) 99%')
        ax.axvline(1 / OBLIQUITY, color='grey', lw=0.7, alpha=0.6)
        ax.text(1 / OBLIQUITY, ax.get_ylim()[1] * 0.02, f' {OBLIQUITY:.0f} kyr',
                fontsize=10, fontfamily='Arial', color='grey', ha='left')
        ax.set_xlabel('Frequency (cycles kyr$^{-1}$)', fontsize=10, fontfamily='Arial')
        ax.set_title(f'{label}', fontsize=10, fontfamily='Arial', loc='left')
        ax.tick_params(axis='both', which='major', labelsize=10, length=3, width=0.5)
        for label in ax.get_xticklabels() + ax.get_yticklabels(): 
            label.set_fontfamily('Arial')
    axes[0].set_ylabel('Power spectral density', fontsize=10, fontfamily='Arial')
    axes[0].legend(fontsize=10, frameon=True, prop={'family': 'Arial'},
                   edgecolor='black', fancybox=False)
    axes[0].text(-0.14, 0.97, 'a', transform=axes[0].transAxes, fontsize=10, fontweight='bold', fontname='Arial')
    axes[1].text(-0.14, 0.97, 'b', transform=axes[1].transAxes, fontsize=10, fontweight='bold', fontname='Arial')
    plt.tight_layout()
    # plt.savefig('figS7_MgO_spectrum.png', dpi=300)
    plt.savefig('figS7_MgO_spectrum.pdf')
    with open('figS7_spectra.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=csv_rows[0].keys())
        w.writeheader()
        w.writerows(csv_rows)
    print('\nSaved: figS7_MgO_spectrum.png/.pdf, figS7_spectra.csv')


if __name__ == '__main__':
    main()

