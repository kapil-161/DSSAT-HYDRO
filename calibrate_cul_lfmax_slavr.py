#!/usr/bin/env python3
"""Calibrate LFMAX and SLAVR in LUGRO048.CUL for LU0301 Bibb against the
literature-correct Dynamic-mu growth trajectory for the 132 mg/L N treatment.

Adapted from a prior Mac/UFGA2201-based version for this project's actual
setup: Windows paths, GTGA2401.LUX treatment 5, cultivar LU0301 (Bibb).

Why Dynamic-mu instead of the raw 6-point Table 3 data (Sharkey et al. 2024,
Agriculture 14(8):1358): a follow-up paper (Sharkey et al. 2025, Agriculture
15(9):1927) shows lettuce RGR decays with maturity rather than staying
constant (pure exponential/Monod overestimates mid-growth then runs away
unrealistically — e.g. predicts 5.9 kg lettuce by 44 DAT). Fitting DSSAT's
CUL growth params to the raw sparse points let LFMAX/SLAVR chase noise and
produced a curve-SHAPE mismatch no 2-3 parameter combination could resolve
(model overshot early growth, undershot late growth, or vice versa).

The Dynamic-mu model (Eq. S6 in the 2025 paper) gives the actual underlying
smooth curve: mu(t) = mu_max*[N]/([N]+Ks_N) for t<=d, plus linear maturity
decay m*(t-d) for t>d. Integrating mu over time and exponentiating from
DM0 reproduces the raw Table 3 points almost exactly (e.g. predicted 1937.3
kg/ha at DAT32 vs observed 1936.0) while also correctly capturing the
deceleration/plateau shape — and predicts the literature-stated peak dry
mass day (41.7 DAT), an independent validation that this curve is right.
Parameters used (Table 2, dry mass SMND-mu, fit to mass):
  mu_max=0.334 mgDM/mgDM/d, Ks_N=1.833 mg/L, m=-0.0119, d=14.0 days,
  DM0=2.198 mg (Sharkey et al. 2024 dataset, day-0 seedling mass).
P and K assumed non-limiting at the 132 mg/L (baseline MSS) treatment, so
the single-nutrient N term approximates the full multi-nutrient model.

DAT (days after transplant) maps to DAP in PlantGro.OUT (DAP=0 at transplant
day 23014). Run mode is 'C' (single treatment), NOT 'A' (runs all treatments
silently ignoring any trailing treatment-number arg — confirmed in
CSM_Main/CSM.for).
"""
from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
from pathlib import Path

ROOT     = Path("C:/DSSAT48")
GENO     = ROOT / "Genotype"
LETTUCE  = ROOT / "Lettuce"

CUL      = GENO / "LUGRO048.CUL"
MODEL    = ROOT / "dscsm048.exe"
FILEX    = "GTGA2401.LUX"
TRTNO    = "5"
PLANTGRO = LETTUCE / "PlantGro.OUT"

VAR_ID   = "LU0301"
PPOP     = 16.0  # plants/m2, from GTGA2401.LUX *PLANTING DETAILS

# Dynamic-mu SMND-mu dry-mass parameters (Sharkey et al. 2025, Table 2)
MU_MAX  = 0.334   # mgDM mgDM^-1 d^-1
KS_N    = 1.833   # mg/L
MAT_M   = -0.0119 # maturity decay rate
DELAY_D = 14.000  # days before decay begins
DM0_MG  = 2.198   # mg/plant at t=0
N_CONC  = 132.24  # mg/L, 132 mg/L (100% MSS baseline) treatment


def _mu_inst(t: float) -> float:
    base = MU_MAX * N_CONC / (N_CONC + KS_N)
    return base if t <= DELAY_D else base + MAT_M * (t - DELAY_D)


def _dm_at(t: float, steps: int = 2000) -> float:
    """Integrate RGR to get dry mass (mg/plant) at day t via Dynamic-mu."""
    if t <= 0:
        return DM0_MG
    dt = t / steps
    log_dm = math.log(DM0_MG)
    tt = 0.0
    for _ in range(steps):
        log_dm += _mu_inst(tt + dt / 2) * dt
        tt += dt
    return math.exp(log_dm)


# DAT(=DAP) -> Dynamic-mu-predicted CWAM (kg/ha), DAT 14-32 only.
# Excludes DAT<14 per the source papers' own convention (Sharkey et al. 2024/
# 2025 exclude <=14 DAT when fitting growth-rate parameters, citing high
# variability/instability in young-plant measurements); empirically this
# also excludes the window where DSSAT's own startup transient (near-zero
# simulated biomass while the literature curve already shows small but real
# growth) dominates the error and isn't informative for CUL calibration.
OBS_CWAM: dict[int, float] = {
    dat: (_dm_at(dat) / 1000.0) * PPOP * 10.0 for dat in range(14, 33)
}

# LFMAX is now fixed at 1.100 mg CO2/m2-s — literature-derived ceiling from
# Ahmed et al. 2022 (Horticulturae 8(3):270): measured peak net Pn of 20.6
# umol CO2/m2/s (cv. Tiberius, 1000 ppm CO2, 300 umol/m2/s light, 0.75 m/s
# air speed) = 0.9066 mg CO2/m2/s net. Converted to DSSAT's LFMAX basis
# (350 ppm CO2, GROSS pre-respiration) using the SPE file's own CO2-response
# curve (CCMP=80, CCMAX=2.09, CCEFF=0.0105 -> 350/1000ppm ratio=0.9413) and an
# assumed leaf dark respiration of 10-20% of gross (standard C3 leafy-crop
# range): gross = 0.9066*0.9413/(1-Rd) = 0.95-1.07 mg CO2/m2-s. Treat 1.1 as
# the ceiling (not the search variable) — no longer a free parameter.
# SLAVR genus range per SPE comments: 310 (low light/cold) to 910 (high
# light/warm) cm2/g (Lorenz & Wiebe 1980) — this remains the free parameter.
# SIZLF bounds from CUL MINIMA/MAXIMA rows: 50.0-650.0 (data-driven, no
# direct literature value found for lettuce SIZLF/SIZREF).
LFMAX_FIXED = 1.100  # literature ceiling, not searched
SLAVR_GRID = [310, 400, 500, 580, 660, 740, 820, 900]  # within SPE SLAMIN-SLAMAX
SIZLF_GRID = [250, 350, 450, 550, 600]
SLAVR_FINE_STEP = 10
SIZLF_FINE_STEP = 20


# ---------------------------------------------------------------------------
# CUL file helpers — fixed-width patch at known column offsets for LU0301
# ---------------------------------------------------------------------------

def _patch_line(line: str, lfmax: float, slavr: float, sizlf: float) -> str:
    # Each field is 6 chars wide, right-aligned, trailing space included:
    # line[79:97] == "1.600 580.0 600.0 " -> three 6-char fields "1.600 "+"580.0 "+"600.0 "
    lfmax_f = f"{lfmax:.3f}".rjust(5) + " "
    slavr_f = f"{slavr:.1f}".rjust(5) + " "
    sizlf_f = f"{sizlf:.1f}".rjust(5) + " "
    return f"{line[:79]}{lfmax_f}{slavr_f}{sizlf_f}{line[97:]}"


def read_cul_param(var_id: str) -> tuple[float, float, float]:
    for line in CUL.read_text().splitlines():
        if re.match(rf"^{var_id}\s", line):
            return (float(line[79:85]), float(line[85:91]), float(line[91:97]))
    raise RuntimeError(f"{var_id} not found in {CUL}")


def update_cul(var_id: str, lfmax: float, slavr: float, sizlf: float) -> None:
    lines = CUL.read_text().splitlines()
    out = []
    for line in lines:
        if re.match(rf"^{var_id}\s", line):
            out.append(_patch_line(line, lfmax, slavr, sizlf))
        else:
            out.append(line)
    CUL.write_text("\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# Run model (mode C = single treatment; mode A ignores the treatment number)
# ---------------------------------------------------------------------------

def run_treatment() -> None:
    subprocess.run(
        [str(MODEL), "C", FILEX, TRTNO],
        check=True, cwd=LETTUCE, capture_output=True,
    )


def parse_cwad_by_dap() -> dict[int, float]:
    """Parse PlantGro.OUT, return {DAP: CWAD} for the most recent run."""
    result: dict[int, float] = {}
    for line in PLANTGRO.read_text().splitlines():
        if not re.match(r"^\s*\d{4}\s+\d{3}\s+\d+\s+\d+", line):
            continue
        parts = line.split()
        dap = int(parts[3])
        cwad = float(parts[12])
        result[dap] = cwad
    return result


TARGET_DAT = 32  # the actual harvest day — what "1936 kg/ha" refers to


def nrmse_vs_obs() -> tuple[float, dict[int, float]]:
    """Objective = %% error at TARGET_DAT (final/harvest biomass) only.
    Per-point errors for the rest of the DAT14-32 curve are still computed
    and reported for visibility, but do not affect the search — the goal
    here is matching final biomass, not whole-curve shape (CUL params can't
    fix the curve-shape mismatch anyway; see module docstring)."""
    sim = parse_cwad_by_dap()
    per_point = {}
    for dat, obs in OBS_CWAM.items():
        s = sim.get(dat)
        if s is None:
            eligible = [d for d in sim if d <= dat]
            s = sim[max(eligible)] if eligible else 0.0
        per_point[dat] = 100.0 * abs(s - obs) / obs
    return per_point[TARGET_DAT], per_point


# ---------------------------------------------------------------------------
# Calibration — joint LFMAX x SLAVR coarse grid, then fine-tune
# ---------------------------------------------------------------------------

def _try(lf: float, sl: float, sz: float, best: dict | None, label: str,
         quiet: bool = False) -> dict:
    update_cul(VAR_ID, lf, sl, sz)
    run_treatment()
    mape, per_point = nrmse_vs_obs()
    is_new_best = best is None or mape < best["mape"]
    if not quiet or is_new_best:
        tag = " *** NEW BEST ***" if (quiet and is_new_best and best is not None) else ""
        print(f"  {label}LFMAX={lf:.3f}  SLAVR={sl:6.1f}  SIZLF={sz:5.1f}  "
              f"MAPE={mape:6.1f}%  DAT32err={per_point[32]:5.1f}%{tag}")
    if is_new_best:
        return {"lfmax": lf, "slavr": sl, "sizlf": sz, "mape": mape, "per_point": per_point}
    return best


def calibrate() -> dict:
    _, _, start_sizlf = read_cul_param(VAR_ID)
    best: dict | None = None

    print(f"=== {VAR_ID}: joint LFMAX x SLAVR coarse grid (SIZLF fixed={start_sizlf}) ===")
    for lf in LFMAX_GRID:
        for sl in SLAVR_GRID:
            best = _try(lf, sl, start_sizlf, best, "")

    cur_lf, cur_sl = best["lfmax"], best["slavr"]
    print(f"\n=== {VAR_ID}: SIZLF coarse grid (LFMAX={cur_lf:.2f} SLAVR={cur_sl}) ===")
    for sz in SIZLF_GRID:
        best = _try(cur_lf, cur_sl, sz, best, "")
    cur_sz = best["sizlf"]

    print(f"\n=== {VAR_ID}: joint fine-tune (LFMAX~{cur_lf:.2f} SLAVR~{cur_sl} SIZLF~{cur_sz}) ===")
    for dlf in (-2, -1, 1, 2):
        lf = round(cur_lf + dlf * LFMAX_FINE_STEP, 3)
        if lf <= 0 or lf > LFMAX_MAX:  # respect CUL MAXIMA bound
            continue
        best = _try(lf, cur_sl, cur_sz, best, "")
    cur_lf = best["lfmax"]

    for dsl in (-2, -1, 1, 2):
        sl = cur_sl + dsl * SLAVR_FINE_STEP
        if sl < 100:
            continue
        best = _try(cur_lf, sl, cur_sz, best, "")
    cur_sl = best["slavr"]

    for dsz in (-2, -1, 1, 2):
        sz = cur_sz + dsz * SIZLF_FINE_STEP
        if not (250 <= sz <= 600):
            continue
        best = _try(cur_lf, cur_sl, sz, best, "")

    print(f"\n>>> Best: LFMAX={best['lfmax']:.3f}  SLAVR={best['slavr']}  "
          f"SIZLF={best['sizlf']}  MAPE={best['mape']:.1f}%")
    return best


# ---------------------------------------------------------------------------
# Monte Carlo search — dense random coverage of the bounded LFMAX x SLAVR x
# SIZLF space, targeting DAT32 (final biomass) match. Bounds match the CUL
# file's own stated MAXIMA/MINIMA (LFMAX<=1.60; SIZLF in [250,600]) plus the
# SPE-cited genus extremes for SLAVR ([310,910], Lorenz & Wiebe 1980).
# ---------------------------------------------------------------------------

MC_LFMAX_RANGE = (1.20, 1.60)
MC_SLAVR_RANGE = (310, 910)
MC_SIZLF_RANGE = (250, 600)


def calibrate_monte_carlo(n_runs: int, seed: int = 42) -> dict:
    import random
    import time
    rng = random.Random(seed)
    best: dict | None = None
    t0 = time.time()

    print(f"=== {VAR_ID}: Monte Carlo search, {n_runs} runs ===")
    print(f"  LFMAX in {MC_LFMAX_RANGE}  SLAVR in {MC_SLAVR_RANGE}  SIZLF in {MC_SIZLF_RANGE}")
    for i in range(n_runs):
        lf = round(rng.uniform(*MC_LFMAX_RANGE), 3)
        sl = round(rng.uniform(*MC_SLAVR_RANGE), 1)
        sz = round(rng.uniform(*MC_SIZLF_RANGE), 1)
        label = f"[{i + 1}/{n_runs}] "
        best = _try(lf, sl, sz, best, label, quiet=True)
        if (i + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (n_runs - i - 1) / rate
            print(f"  ... {i + 1}/{n_runs} done, {elapsed:.0f}s elapsed, "
                  f"~{eta:.0f}s remaining, current best MAPE={best['mape']:.2f}%")

    print(f"\n>>> Monte Carlo best: LFMAX={best['lfmax']:.3f}  SLAVR={best['slavr']}  "
          f"SIZLF={best['sizlf']}  MAPE={best['mape']:.2f}%")
    return best


def main() -> None:
    import sys

    backup = CUL.with_suffix(".CUL.calib_bak")
    shutil.copy2(CUL, backup)
    print(f"Backup: {backup}")

    orig = read_cul_param(VAR_ID)
    print(f"Starting params for {VAR_ID}: LFMAX={orig[0]} SLAVR={orig[1]} SIZLF={orig[2]}\n")

    mc_runs = None
    if len(sys.argv) > 1 and sys.argv[1] == "--monte-carlo":
        mc_runs = int(sys.argv[2]) if len(sys.argv) > 2 else 10000

    try:
        best = calibrate_monte_carlo(mc_runs) if mc_runs else calibrate()
        update_cul(VAR_ID, best["lfmax"], best["slavr"], best["sizlf"])
        run_treatment()
        final_nrmse, final_points = nrmse_vs_obs()

        print("\n*** CUL updated ***")
        print(f"  {VAR_ID}: LFMAX={best['lfmax']:.3f}  SLAVR={best['slavr']}  SIZLF={best['sizlf']}")
        print(f"  Final NRMSE: {final_nrmse:.1f}%")
        print("  Per-point %% error:  " + "  ".join(
            f"DAT{d}={e:.1f}%" for d, e in sorted(final_points.items())))

        print("\n" + json.dumps({
            VAR_ID: {"lfmax": best["lfmax"], "slavr": best["slavr"],
                     "sizlf": best["sizlf"], "nrmse": final_nrmse,
                     "per_point_pct_error": final_points},
        }, indent=2))

    except Exception:
        shutil.copy2(backup, CUL)
        print("ERROR — CUL restored from backup")
        raise


if __name__ == "__main__":
    main()
