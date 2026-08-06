#!/usr/bin/env python3
"""Calibrate LFMAX / SLAVR in LUGRO048.CUL against final observed CWAM
(kg/ha) for all 9 optimum-treatment cultivars (DSSBatch_OPT.v48).

Generalized from a Bibb-only, literature-growth-curve-fitting version:
this one calibrates every cultivar in the optimum-treatment table against
its own single final observed CWAM value, and reads its LFMAX/SLAVR search
bounds directly from LUGRO048.CUL's own MINIMA/MAXIMA rows (VAR# 999991 /
999992) instead of hardcoding per-cultivar literature bounds -- edit those
two rows in the CUL file and this script's search space changes with them.

The 9 (FILEX, TRTNO) pairs are read directly from DSSBatch_OPT.v48 itself
(not duplicated here) so they can never drift out of sync with the batch
file used everywhere else in this project. Only VAR_ID and the observed
CWAM (which aren't in the batch file) are kept in TREATMENT_INFO below, in
the same order as the batch file's lines.

Each grid point runs a single treatment via a throwaway 1-line batch file
(mode B), then reads CWAMS (simulated CWAM) from Evaluate.OUT -- the same
field ($28, 1-indexed / [27] 0-indexed) used by hand all session. Mode B
with a real DSSBatch-style line is used instead of mode C + PlantGro.OUT
so a single, already-proven parsing path covers every experiment.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT    = Path("C:/DSSAT48")
GENO    = ROOT / "Genotype"
LETTUCE = ROOT / "Lettuce"

CUL       = GENO / "LUGRO048.CUL"
MODEL     = ROOT / "dscsm048.exe"
OPT_BATCH = LETTUCE / "DSSBatch_OPT.v48"
TMP_BATCH = LETTUCE / "_calib_tmp.v48"
EVALUATE  = LETTUCE / "Evaluate.OUT"

# Fixed-width columns for LFMAX/SLAVR/SIZLF -- confirmed identical across
# every cultivar row (incl. MINIMA/MAXIMA) in LUGRO048.CUL.
LFMAX_COL = slice(79, 85)
SLAVR_COL = slice(85, 91)
SIZLF_COL = slice(91, 97)

# VAR_ID and observed final CWAM (kg/ha) for each of the 9 optimum
# treatments, in the SAME ORDER as the data lines in DSSBatch_OPT.v48.
# obs values and cultivar assignments per project_optimum_treatments.
TREATMENT_INFO = [
    # label,       VAR_ID,   obs_cwam
    ("Bibb",       "LU0301", 1936.0),
    ("Rex",        "LU0001", 1420.0),
    ("Muir",       "LU0002", 1165.0),
    ("Skyphos",    "LU0003", 1677.0),
    ("BG23-1251",  "LU0201", 1796.9),
    ("Waldmanns",  "LU0202", 1826.8),
    ("SITONIA",    "LU0004", 2862.1),
    ("Salvius",    "LU0006", 1350.3),
    ("MC",         "LU0007", 1551.4),
]

# Coarse-grid resolution for the joint LFMAX x SLAVR search, then a fine
# tune pass around the coarse best. SIZLF is left at its current per
# cultivar value -- it shapes leaf-size distribution, not total biomass.
LFMAX_STEPS = 5
SLAVR_STEPS = 6
FINE_LFMAX_STEP = 0.02
FINE_SLAVR_STEP = 10
FINE_PASSES = 3


def read_bounds() -> dict[str, float]:
    """Read MINIMA/MAXIMA rows from the live CUL file -- bounds are
    whatever is currently in the file, not a hardcoded literature number.
    Editing the MINIMA/MAXIMA rows in LUGRO048.CUL changes what this
    script's search space is, next run."""
    bounds: dict[str, float] = {}
    for line in CUL.read_text().splitlines():
        if line.startswith("999991"):  # MINIMA
            bounds["lfmax_min"] = float(line[LFMAX_COL])
            bounds["slavr_min"] = float(line[SLAVR_COL])
        elif line.startswith("999992"):  # MAXIMA
            bounds["lfmax_max"] = float(line[LFMAX_COL])
            bounds["slavr_max"] = float(line[SLAVR_COL])
    missing = {"lfmax_min", "lfmax_max", "slavr_min", "slavr_max"} - bounds.keys()
    if missing:
        raise RuntimeError(f"MINIMA/MAXIMA rows missing fields: {missing}")
    return bounds


def read_opt_batch_template() -> tuple[str, list[str]]:
    """Split DSSBatch_OPT.v48 into (preamble, data_lines). The preamble is
    everything through the @FILEX header line, kept byte-for-byte -- a
    minimal hand-built 3-line batch (just $BATCH + @FILEX + one data line)
    reproducibly fails with "Error in the format of the batch file" from
    dscsm048.exe, so the real file's full comment block is reused rather
    than reconstructed."""
    all_lines = OPT_BATCH.read_text().splitlines()
    preamble, data_lines = [], []
    in_data = False
    for l in all_lines:
        if l.lstrip().startswith("@FILEX"):
            preamble.append(l)
            in_data = True
        elif in_data:
            if l.strip():
                data_lines.append(l)
        else:
            preamble.append(l)
    if len(data_lines) != len(TREATMENT_INFO):
        raise RuntimeError(
            f"DSSBatch_OPT.v48 has {len(data_lines)} treatment lines but "
            f"TREATMENT_INFO has {len(TREATMENT_INFO)} -- keep them in sync")
    return "\n".join(preamble) + "\n", data_lines


def _patch_line(line: str, lfmax: float, slavr: float, sizlf: float) -> str:
    lfmax_f = f"{lfmax:.3f}".rjust(5) + " "
    slavr_f = f"{slavr:.1f}".rjust(5) + " "
    sizlf_f = f"{sizlf:.1f}".rjust(5) + " "
    return f"{line[:79]}{lfmax_f}{slavr_f}{sizlf_f}{line[97:]}"


def read_cul_param(var_id: str) -> tuple[float, float, float]:
    for line in CUL.read_text().splitlines():
        if re.match(rf"^{var_id}\s", line):
            return (float(line[LFMAX_COL]), float(line[SLAVR_COL]), float(line[SIZLF_COL]))
    raise RuntimeError(f"{var_id} not found in {CUL}")


def update_cul(var_id: str, lfmax: float, slavr: float, sizlf: float) -> None:
    lines = CUL.read_text().splitlines()
    out, found = [], False
    for line in lines:
        if re.match(rf"^{var_id}\s", line):
            out.append(_patch_line(line, lfmax, slavr, sizlf))
            found = True
        else:
            out.append(line)
    if not found:
        raise RuntimeError(f"{var_id} not found in {CUL}")
    CUL.write_text("\n".join(out) + "\n")


def run_and_get_cwam(preamble: str, batch_line: str) -> float:
    """Run one treatment (via a throwaway 1-data-line batch, mode B) and
    return its simulated CWAM (kg/ha) from Evaluate.OUT.

    write_bytes (not write_text) is deliberate: DSSBatch_OPT.v48 is LF-only
    and dscsm048.exe's batch parser fails ("Error in the format of the
    batch file") on the CRLF that Path.write_text() introduces via
    Windows' universal-newline translation.

    Passing the batch file as a bare relative filename (with cwd=LETTUCE)
    is deliberate, not cosmetic: dscsm048.exe reproducibly misparses an
    absolute-path filename argument -- its own error message prints the
    filename truncated ("_calib_tmp." with ".v48" cut off), which is why
    it always reported "Error in the format of the batch file" even for
    byte-identical, correctly-formatted content. Every other working
    invocation in this project (DSSBatch_OPT.v48 etc.) already used a bare
    relative filename; this just matches that convention."""
    TMP_BATCH.write_bytes((preamble + batch_line + "\n").encode("utf-8"))
    result = subprocess.run(
        [str(MODEL), "CRGRO048", "B", TMP_BATCH.name],
        cwd=LETTUCE, capture_output=True,
    )
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, result.args, output=result.stdout, stderr=result.stderr)
    data_lines = [l for l in EVALUATE.read_text().splitlines()
                  if re.match(r"^\s*\d+\s+\S+\s+\d+", l)]
    if not data_lines:
        raise RuntimeError("No data row found in Evaluate.OUT after run")
    parts = data_lines[-1].split()
    return float(parts[27])  # CWAMS (simulated), 0-indexed == awk $28


def _try(preamble: str, batch_line: str, obs: float, lf: float, sl: float, sz: float,
         best: dict | None, var_id: str, label: str) -> dict:
    update_cul(var_id, lf, sl, sz)
    sim = run_and_get_cwam(preamble, batch_line)
    err = 100.0 * abs(sim - obs) / obs
    is_new_best = best is None or err < best["err"]
    tag = " *** NEW BEST ***" if is_new_best and best is not None else ""
    print(f"  {label}LFMAX={lf:.3f}  SLAVR={sl:6.1f}  sim={sim:7.1f}  "
          f"obs={obs:7.1f}  err={err:5.1f}%{tag}")
    if is_new_best:
        return {"lfmax": lf, "slavr": sl, "sizlf": sz, "sim": sim, "err": err}
    return best


def calibrate_one(label: str, var_id: str, obs: float, preamble: str, batch_line: str,
                   bounds: dict[str, float]) -> dict:
    _, _, start_sz = read_cul_param(var_id)
    lf_lo, lf_hi = bounds["lfmax_min"], bounds["lfmax_max"]
    sl_lo, sl_hi = bounds["slavr_min"], bounds["slavr_max"]

    print(f"\n=== {label} ({var_id}): obs={obs} kg/ha  "
          f"bounds LFMAX[{lf_lo},{lf_hi}] SLAVR[{sl_lo},{sl_hi}] ===")

    lf_grid = [round(lf_lo + i * (lf_hi - lf_lo) / (LFMAX_STEPS - 1), 3)
               for i in range(LFMAX_STEPS)]
    sl_grid = [round(sl_lo + i * (sl_hi - sl_lo) / (SLAVR_STEPS - 1), 1)
               for i in range(SLAVR_STEPS)]

    best: dict | None = None
    for lf in lf_grid:
        for sl in sl_grid:
            best = _try(preamble, batch_line, obs, lf, sl, start_sz, best, var_id, "coarse ")

    cur_lf, cur_sl = best["lfmax"], best["slavr"]
    for _ in range(FINE_PASSES):
        improved = False
        for dlf in (-1, 1):
            lf = round(cur_lf + dlf * FINE_LFMAX_STEP, 3)
            if not (lf_lo <= lf <= lf_hi):
                continue
            nb = _try(preamble, batch_line, obs, lf, cur_sl, start_sz, best, var_id, "fine   ")
            if nb is not best:
                best, improved = nb, True
        for dsl in (-1, 1):
            sl = round(cur_sl + dsl * FINE_SLAVR_STEP, 1)
            if not (sl_lo <= sl <= sl_hi):
                continue
            nb = _try(preamble, batch_line, obs, cur_lf, sl, start_sz, best, var_id, "fine   ")
            if nb is not best:
                best, improved = nb, True
        cur_lf, cur_sl = best["lfmax"], best["slavr"]
        if not improved:
            break

    best["sizlf"] = start_sz
    # _try() leaves the CUL file at whatever combination it last probed,
    # which is not necessarily the best one found -- re-apply the winner
    # explicitly so the file on disk always matches what's reported.
    update_cul(var_id, best["lfmax"], best["slavr"], best["sizlf"])
    print(f">>> {label} best: LFMAX={best['lfmax']:.3f}  SLAVR={best['slavr']}  "
          f"SIZLF={start_sz}  sim={best['sim']:.1f}  obs={obs}  err={best['err']:.1f}%")
    return best


def main() -> None:
    import sys

    backup = CUL.with_suffix(".CUL.calib_bak")
    shutil.copy2(CUL, backup)
    print(f"Backup: {backup}")

    bounds = read_bounds()
    print("Bounds from CUL MINIMA/MAXIMA rows:", bounds)

    preamble, batch_lines = read_opt_batch_template()
    only = sys.argv[1] if len(sys.argv) > 1 else None

    results: dict[str, dict] = {}
    try:
        for (label, var_id, obs), batch_line in zip(TREATMENT_INFO, batch_lines):
            if only and label.lower() != only.lower():
                continue
            results[label] = calibrate_one(label, var_id, obs, preamble, batch_line, bounds)

        print("\n=== FINAL SUMMARY ===")
        for label, r in results.items():
            print(f"  {label:12s} LFMAX={r['lfmax']:.3f}  SLAVR={r['slavr']:6.1f}  "
                  f"SIZLF={r['sizlf']:6.1f}  sim={r['sim']:7.1f}  err={r['err']:5.1f}%")

        print("\n" + json.dumps(results, indent=2))

    except Exception:
        shutil.copy2(backup, CUL)
        print("ERROR - CUL restored from backup")
        raise
    finally:
        TMP_BATCH.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
