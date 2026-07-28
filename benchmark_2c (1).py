#!/usr/bin/env python3
"""
Section 2C benchmark: detection-timing comparison on real 2014 West Africa
Ebola national case-reporting (OCHA ROWCA / HDX, data-ebola-public.xlsx).

For each country we compare, for four detectors, the day on which each first
flags the outbreak's growth-rate acceleration, under three observation schedules
built on the SAME underlying (PCHIP-interpolated) epidemic curve:

  real     - the actual historical reporting dates/values (informative sampling)
  regular  - the same curve resampled at evenly spaced dates (n matched)
  shuffled - the real inter-report gaps, order randomized: same overall
             irregularity as real, but no sparse->dense regime-switch structure
             (isolates the regime-switch mechanism from generic noise)

Detectors: WHO 75th/mean+2SD upper threshold, one-sided upper CUSUM,
2-state Gaussian HMM, and our variance-trend (Kendall-tau on rolling variance)
method. Signal = new cases per day between consecutive reports.

Data: OCHA ROWCA compiled sub-national series, Humanitarian Data Exchange.
Download data-ebola-public.xlsx and set ROWCA_XLSX to its path (default: cwd).
method) and regenerates the Liberia and Sierra Leone figures from the same code.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats, interpolate

XLSX = os.environ.get("ROWCA_XLSX", "data-ebola-public.xlsx")
SHEET = "ROWCA Ebola All Sec Review"
N_SHUFFLES = 20
SHUFFLE_SEED = 42


# ---------- detectors (recovered verbatim from the 2026-07-13 session) ----------
def growth_rate_series(t, cumcases):
    dt = np.diff(t)
    dc = np.diff(cumcases)
    rate = dc / np.maximum(dt, 1)
    return t[1:], rate


def who_upper_detect(rate, t, baseline_n=15):
    baseline = rate[:baseline_n]
    threshold = baseline.mean() + 2 * baseline.std()
    post = rate[baseline_n:]
    above = post > threshold
    return t[baseline_n:][np.argmax(above)] if above.any() else None


def cusum_upper_detect(rate, t, baseline_n=15, h=4.0, k=0.5):
    baseline = rate[:baseline_n]
    mu0, sigma0 = baseline.mean(), baseline.std()
    if sigma0 < 1e-9:
        return None
    z = (rate[baseline_n:] - mu0) / sigma0
    s = 0.0
    for i, zi in enumerate(z):
        s = max(0, s + zi - k)
        if s > h:
            return t[baseline_n:][i]
    return None


def hmm_detect_time(rate, t):
    try:
        from hmmlearn import hmm
        model = hmm.GaussianHMM(n_components=2, covariance_type="diag",
                                n_iter=100, random_state=0)
        model.fit(rate.reshape(-1, 1))
        states = model.predict(rate.reshape(-1, 1))
        high = np.argmax(model.means_.flatten())
        sustained = np.where((states == high) & (np.roll(states, -1) == high))[0]
        return t[sustained[0]] if len(sustained) else None
    except Exception:
        return None


def our_method_detect_time(rate, t, win_frac=0.3):
    n = len(rate)
    win = max(5, int(win_frac * n))
    variances = np.array([rate[i - win:i].var() for i in range(win, n)])
    for i in range(10, len(variances)):
        tau, p = stats.kendalltau(np.arange(i), variances[:i])
        if tau > 0 and p < 0.05:
            return t[win + i]
    return None


# ---------- data prep ----------
def country_series(df, country):
    sub = df[(df["Country"] == country) &
             (df["Localite"] == "National") &
             (df["Category"] == "Cases")].sort_values("Date")
    sub = sub.copy()
    sub["Value"] = pd.to_numeric(sub["Value"], errors="coerce")
    sub = sub.dropna(subset=["Value", "Date"]).drop_duplicates(subset="Date")
    sub = sub[sub["Date"] > pd.Timestamp("2013-06-01")]
    t0 = sub["Date"].min()
    days = (sub["Date"] - t0).dt.days.values.astype(int)
    cases = sub["Value"].values.astype(float)
    # enforce monotonic day ordering
    order = np.argsort(days)
    return days[order], cases[order]


def run_country(df, country):
    days, cases = country_series(df, country)
    interp = interpolate.PchipInterpolator(days, cases)

    # real
    real_t, real_cases = days, cases
    # regular counterfactual
    reg_t = np.linspace(days.min(), days.max(), len(days)).astype(int)
    reg_cases = interp(reg_t)

    def detect_all(t, c):
        rt, rate = growth_rate_series(t, c)
        return {
            "WHO": who_upper_detect(rate, rt),
            "CUSUM": cusum_upper_detect(rate, rt),
            "HMM": hmm_detect_time(rate, rt),
            "Ours": our_method_detect_time(rate, rt),
        }

    real = detect_all(real_t, real_cases)
    regular = detect_all(reg_t, reg_cases)

    # shuffled-gap control for our method (mean over N shuffles)
    rng = np.random.default_rng(SHUFFLE_SEED)
    gaps = np.diff(real_t)
    ours_shuffled = []
    for _ in range(N_SHUFFLES):
        sg = rng.permutation(gaps)
        st = np.concatenate([[real_t[0]], real_t[0] + np.cumsum(sg)]).astype(int)
        st = np.clip(st, days.min(), days.max())
        sc = interp(st)
        rt, rate = growth_rate_series(st, sc)
        d = our_method_detect_time(rate, rt)
        if d is not None:
            ours_shuffled.append(d)
    shuf_mean = np.mean(ours_shuffled) if ours_shuffled else float("nan")
    shuf_rng = (min(ours_shuffled), max(ours_shuffled)) if ours_shuffled else (None, None)

    return real, regular, shuf_mean, shuf_rng, len(ours_shuffled)


def main():
    df = pd.read_excel(XLSX, sheet_name=SHEET)
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")

    for country in ["Guinea", "Liberia", "Sierra Leone"]:
        real, regular, shuf_mean, shuf_rng, nsh = run_country(df, country)
        print(f"=== {country} ===")
        print(f"{'method':<8}{'real':>8}{'regular':>10}")
        for m in ["WHO", "CUSUM", "HMM", "Ours"]:
            print(f"{m:<8}{str(real[m]):>8}{str(regular[m]):>10}")
        print(f"Ours shuffled-gap control: mean day {shuf_mean:.0f} "
              f"(range {shuf_rng[0]}-{shuf_rng[1]}, {nsh}/{N_SHUFFLES} detected)")
        if isinstance(real['Ours'], (int, np.integer)) and isinstance(regular['Ours'], (int, np.integer)):
            print(f"  -> Ours delay real vs regular: "
                  f"{real['Ours'] - regular['Ours']:+d} days")
        print()


if __name__ == "__main__":
    main()
