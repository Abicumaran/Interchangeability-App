"""Core statistical engine for the PROXIMA trueness + Bland–Altman Streamlit app.

The implementation is derived directly from the validated notebook
`PROXIMA_Trueness_and_BlandAltman_v2.ipynb`.

The two outlier branches remain deliberately separate:
1. Trueness regression: validated generalized ESD on raw replicate rows.
2. Paired-specimen Bland–Altman: global/analyte flags, Shapiro–Wilk,
   manual Grubbs Gcrit if normal, or MAD modified-Z if non-normal,
   with at most one linked replicate removed per analyte before donor averaging.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import pearsonr, t
from sklearn.linear_model import HuberRegressor
from sklearn.metrics import mean_squared_error, r2_score


ANALYTES: dict[str, dict[str, object]] = {
    "RBC":   {"mhs": "RBC",     "ref": "RBC_ref",   "normal": (3.85, 5.65),   "unit": "10^12/L"},
    "WBC":   {"mhs": "WBC_2",   "ref": "WBC_ref",   "normal": (3.60, 10.50),  "unit": "10^9/L"},
    "NEUT":  {"mhs": "NEUT_2",  "ref": "NEUT_ref",  "normal": (1.50, 7.70),   "unit": "10^9/L"},
    "LYMPH": {"mhs": "LYMPH_2", "ref": "LYMPH_ref", "normal": (1.10, 4.00),   "unit": "10^9/L"},
    "MXD":   {"mhs": "MXD_2",   "ref": "MXD_ref",   "normal": (0.10, 1.60),   "unit": "10^9/L"},
    "PLT":   {"mhs": "PLT",     "ref": "PLT_ref",   "normal": (160.0, 370.0), "unit": "10^9/L"},
    "MCV":   {"mhs": "MCV",     "ref": "MCV_ref",   "normal": (80.0, 101.0),  "unit": "fL"},
    "RDW":   {"mhs": "RDW",     "ref": "RDW_ref",   "normal": (11.0, 16.0),   "unit": "%"},
    "MCH":   {"mhs": "MCH",     "ref": "MCH_ref",   "normal": (27.0, 34.0),   "unit": "pg"},
    "MCHC":  {"mhs": "MCHC",    "ref": "MCHC_ref",  "normal": (320.0, 360.0), "unit": "g/L"},
    "HGB":   {"mhs": "HGB",     "ref": "HGB_ref",   "normal": (118.0, 172.0), "unit": "g/L"},
    "HCT":   {"mhs": "HCT",     "ref": "HCT_ref",   "normal": (35.0, 49.0),   "unit": "percentage points"},
}

DEFAULT_ANALYSIS_GROUPS: dict[str, list[str]] = {
    "Capillary_Mc_Mc2": ["Mc", "Mc2"],
    "Venous_MK2": ["MK2"],
    "Venous_MK3": ["MK3"],
    "Direct_Capillary_Dir": ["Dir"],
    "TASSO_Tas": ["Tas"],
}


def normalize_bool(value) -> bool | None:
    if pd.isna(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n"}:
        return False
    return None


def parse_blood_sample_id(value: object) -> dict[str, object]:
    """Parse IDs such as 210726D05-Mc-2, retaining D-suffixes such as D14m."""
    text = "" if value is None else str(value).strip()
    match = re.search(r"(?P<donor>D[^-]+)-(?P<specimen>[^-]+)-(?P<replicate>\d+)$", text)
    if not match:
        return {"donor": np.nan, "specimen_type": np.nan, "replicate_number": np.nan, "parse_ok": False}
    gd = match.groupdict()
    return {
        "donor": gd["donor"],
        "specimen_type": gd["specimen"],
        "replicate_number": int(gd["replicate"]),
        "parse_ok": True,
    }


def load_and_prepare(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "batch_id" not in df.columns or "bloodSampleId" not in df.columns or "global_flag" not in df.columns:
        raise ValueError("Input must contain batch_id, bloodSampleId and global_flag columns.")

    parsed = pd.DataFrame([parse_blood_sample_id(v) for v in df["bloodSampleId"]])
    for col in parsed.columns:
        df[col] = parsed[col].values
    df["global_flag_bool"] = df["global_flag"].map(normalize_bool)
    df["batch_id"] = df["batch_id"].astype(str)

    required = []
    for cfg in ANALYTES.values():
        required.extend([str(cfg["mhs"]), str(cfg["ref"])])
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"Missing analyte columns: {missing}")

    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# ---------------------------------------------------------------------------
# Validated APT ESD logic, retained as closely as possible from the source
# notebook. Detection is performed on raw replicate rows before donor means.
# ---------------------------------------------------------------------------

def generalized_esd_with_batch(x, batch_ids, alpha: float = 0.05, max_outliers_frac: float = 0.05) -> dict:
    x = np.asarray(x, dtype=float)
    batch_ids = np.asarray(batch_ids)
    mask = ~np.isnan(x)
    x = x[mask]
    batch_ids = batch_ids[mask]
    n_total = len(x)
    if n_total < 20:
        raise ValueError("ESD requires at least 20 samples")
    h = int(np.floor(max_outliers_frac * n_total))
    if h < 1:
        return {"n_outliers": 0, "outlier_batch_ids": [], "test_statistics": [], "critical_values": []}

    remaining = x.copy()
    remaining_ids = batch_ids.copy()
    esd_stats, lambdas, removed_ids = [], [], []
    for i in range(1, h + 1):
        mean = remaining.mean()
        sd = remaining.std(ddof=1)
        if sd == 0:
            break
        deviations = np.abs(remaining - mean) / sd
        max_idx = int(np.argmax(deviations))
        esd_i = float(deviations[max_idx])
        nu = n_total - i - 1
        p = 1 - alpha / (2 * (n_total - i + 1))
        t_crit = t.ppf(p, nu)
        lambda_i = ((n_total - i) * t_crit) / np.sqrt((nu + t_crit**2) * (n_total - i + 1))
        esd_stats.append(esd_i)
        lambdas.append(float(lambda_i))
        removed_ids.append(remaining_ids[max_idx])
        remaining = np.delete(remaining, max_idx)
        remaining_ids = np.delete(remaining_ids, max_idx)

    n_outliers = 0
    for i, (esd_i, lambda_i) in enumerate(zip(esd_stats, lambdas), start=1):
        if esd_i > lambda_i:
            n_outliers = i
    return {
        "n_outliers": n_outliers,
        "outlier_batch_ids": [str(v) for v in removed_ids[:n_outliers]],
        "test_statistics": esd_stats,
        "critical_values": lambdas,
    }


def sd_cv_pattern(y_true, y_pred, n_bins: int = 5, min_per_bin: int = 5, eps: float = 1e-8):
    y_true = np.asarray(y_true, float).reshape(-1)
    y_pred = np.asarray(y_pred, float).reshape(-1)
    if y_true.size != y_pred.size:
        raise ValueError("y_true and y_pred must have the same length")
    if not (np.isfinite(y_true).all() and np.isfinite(y_pred).all()):
        raise ValueError("There is some NaN/inf")
    d = y_true - y_pred
    m = 0.5 * (y_true + y_pred)
    edges = np.unique(np.quantile(m, np.linspace(0, 1, n_bins + 1)))
    if edges.size < 4:
        raise ValueError("Too few unique values to create bins")
    bin_id = np.digitize(m, edges[1:-1], right=True)
    rows = []
    for k in range(edges.size - 1):
        mask = bin_id == k
        if mask.sum() < min_per_bin:
            continue
        mk = m[mask]
        dk = d[mask]
        sd = dk.std(ddof=1)
        mean_m = mk.mean()
        cv = sd / max(abs(mean_m), eps)
        rows.append({"bin": k, "range": (edges[k], edges[k + 1]), "n": int(mask.sum()),
                     "mean_level_m": float(mean_m), "sd_diff": float(sd), "cv_local": float(cv)})
    tbl = pd.DataFrame(rows)
    min_valid_bins = max(2, n_bins // 2)
    if tbl.shape[0] < min_valid_bins:
        raise ValueError(f"Very few valid bins (<{min_valid_bins})")
    sd_ratio = tbl["sd_diff"].max() / max(tbl["sd_diff"].min(), eps)
    cv_ratio = tbl["cv_local"].max() / max(tbl["cv_local"].min(), eps)
    corr_sd = np.corrcoef(tbl["mean_level_m"], tbl["sd_diff"])[0, 1]
    corr_cv = np.corrcoef(tbl["mean_level_m"], tbl["cv_local"])[0, 1]
    if sd_ratio <= 2.0 and abs(corr_sd) < 0.3 and corr_cv < -0.3:
        label, lbl = "SD constant (additive)", "SD constant"
    elif cv_ratio <= 2.0 and abs(corr_cv) < 0.3 and corr_sd > 0.3:
        label, lbl = "CV constant (percent)", "CV constant"
    else:
        label, lbl = "Mixed or uncertain (review by low/high sections)", "Mixed"
    return tbl, {
        "sd_ratio_maxmin": float(sd_ratio), "cv_ratio_maxmin": float(cv_ratio),
        "corr(level, sd)": float(corr_sd), "corr(level, cv)": float(corr_cv),
        "diagnosis": label, "lbl": lbl, "edges_used": edges.tolist(),
    }, lbl


def change_point_sd_by_level(y_true, y_pred, n_bins: int = 5, min_per_bin: int = 5,
                             robust: bool = True, eps: float = 1e-8):
    y_true = np.asarray(y_true, float).reshape(-1)
    y_pred = np.asarray(y_pred, float).reshape(-1)
    d = y_true - y_pred
    m = 0.5 * (y_true + y_pred)
    edges = np.unique(np.quantile(m, np.linspace(0, 1, n_bins + 1)))
    if edges.size < 4:
        raise ValueError("Not enough unique m values to create bins")
    bin_id = np.digitize(m, edges[1:-1], right=True)
    rows = []
    for k in range(edges.size - 1):
        mask = bin_id == k
        nk = int(mask.sum())
        if nk < min_per_bin:
            continue
        dk = d[mask]
        mk = m[mask]
        if robust:
            med = np.median(dk)
            mad = np.median(np.abs(dk - med))
            scale = 1.4826 * mad
        else:
            scale = dk.std(ddof=1)
        rows.append({"bin": k, "left": float(edges[k]), "right": float(edges[k + 1]),
                     "n": nk, "mean_level_m": float(mk.mean()), "scale": float(scale)})
    tbl = pd.DataFrame(rows).reset_index(drop=True)
    min_valid_bins = max(2, n_bins // 2)
    if tbl.shape[0] < min_valid_bins:
        raise ValueError(f"Too few valid bins (<{min_valid_bins})")
    log_scales = np.log(tbl["scale"].values + eps)
    jumps = np.abs(np.diff(log_scales))
    cut_pos = int(np.argmax(jumps))
    cutoff_m = float(tbl.loc[cut_pos, "right"])
    return cutoff_m, tbl, {
        "method": "robust MAD scale" if robust else "SD scale",
        "cut_between_bins": (int(tbl.loc[cut_pos, "bin"]), int(tbl.loc[cut_pos + 1, "bin"])),
        "jump_score_abs_log": float(jumps[cut_pos]), "cutoff_m": cutoff_m,
        "edges_used": edges.tolist(),
    }


def run_esd(diff, batch_ids, alpha: float = 0.05, max_outliers_frac: float = 0.05):
    if len(diff) < 20:
        return [], 0, {"n_outliers": 0, "outlier_batch_ids": [], "test_statistics": [], "critical_values": []}
    result = generalized_esd_with_batch(diff, batch_ids, alpha=alpha, max_outliers_frac=max_outliers_frac)
    return result["outlier_batch_ids"], result["n_outliers"], result


def run_binned_esd(df: pd.DataFrame, true_col: str, pred_col: str, batch_col: str = "batch_id",
                   n_bins: int = 10, max_outliers_frac: float = 0.05):
    work = df.copy()
    work["bin"] = pd.qcut(work[true_col], q=n_bins, duplicates="drop")
    outlier_ids, details = [], []
    for b in work["bin"].dropna().unique():
        sub = work[work["bin"] == b]
        if len(sub) < 20:
            continue
        diff = sub[pred_col] - sub[true_col]
        out_bin, _, result = run_esd(diff, sub[batch_col], max_outliers_frac=max_outliers_frac)
        outlier_ids.extend(out_bin)
        details.append({"bin": str(b), "result": result})
    return outlier_ids, len(outlier_ids), details


def compute_diff(df: pd.DataFrame, true_col: str, pred_col: str, method: str) -> pd.Series:
    if method == "SD":
        return df[pred_col] - df[true_col]
    if method == "CV":
        return (df[pred_col] - df[true_col]) / df[true_col] * 100
    raise ValueError("Unknown diff method")


def _validated_esd_single_analyte(work: pd.DataFrame, true_col: str, pred_col: str,
                                  n_bins: int = 2, max_outliers_bins_frac: float = 0.08) -> tuple[list[str], dict]:
    original_n = len(work)
    if original_n < 20:
        return [], {"error_model": "Not run: <20 replicate rows", "n_replicates": original_n,
                    "n_outliers": 0, "cutoff": np.nan, "esd_detail": "[]"}

    transformed = work[["batch_id", true_col, pred_col]].dropna().copy()
    if len(transformed) < 20:
        return [], {"error_model": "Not run: <20 complete replicate rows", "n_replicates": len(transformed),
                    "n_outliers": 0, "cutoff": np.nan, "esd_detail": "[]"}

    # Preserve the original validated transformation exactly.
    transformed[true_col] = np.log1p(transformed[true_col] + 10000)
    transformed[pred_col] = np.log1p(transformed[pred_col] + 10000)

    n_size = len(transformed)
    if n_size >= 300:
        n_size_bins, min_per_bin = 10, 30
    elif n_size >= 100:
        n_size_bins, min_per_bin = 5, 30
    elif n_size >= 50:
        n_size_bins, min_per_bin = 4, 15
    else:
        n_size_bins, min_per_bin = 4, 5

    try:
        _, summary, error_model = sd_cv_pattern(
            transformed[true_col], transformed[pred_col], n_bins=n_size_bins, min_per_bin=min_per_bin
        )
    except Exception as exc:
        summary, error_model = {"diagnosis": f"Fallback due to: {exc}"}, "Mixed"

    try:
        cutoff, _, cp_info = change_point_sd_by_level(
            transformed[true_col], transformed[pred_col], n_bins=n_size_bins,
            min_per_bin=min_per_bin, robust=True
        )
    except Exception as exc:
        cutoff, cp_info = float(transformed[true_col].median()), {"error": str(exc)}

    esd_detail: object = []
    if error_model == "SD constant":
        diff = compute_diff(transformed, true_col, pred_col, "SD")
        outlier_ids, _, esd_detail = run_esd(diff, transformed["batch_id"])
    elif error_model == "CV constant":
        diff = compute_diff(transformed, true_col, pred_col, "CV")
        outlier_ids, _, esd_detail = run_esd(diff, transformed["batch_id"])
    else:
        # This follows the original branching. With its +10000 log transform,
        # the near-zero branch is normally not entered.
        if (transformed[true_col] < 1).any():
            outlier_ids, _, esd_detail = run_binned_esd(
                transformed, true_col, pred_col, n_bins=n_bins,
                max_outliers_frac=max_outliers_bins_frac
            )
        else:
            low = transformed[transformed[true_col] <= cutoff]
            high = transformed[transformed[true_col] > cutoff]
            out_low, _, detail_low = run_esd(compute_diff(low, true_col, pred_col, "SD"), low["batch_id"])
            out_high, _, detail_high = run_esd(compute_diff(high, true_col, pred_col, "CV"), high["batch_id"])
            outlier_ids = list(out_low) + list(out_high)
            esd_detail = {"low": detail_low, "high": detail_high}

    diagnostics = {
        "error_model": error_model,
        "error_model_diagnosis": summary.get("diagnosis", ""),
        "n_replicates": len(transformed),
        "n_outliers": len(outlier_ids),
        "cutoff": cutoff,
        "change_point_info": json.dumps(cp_info, default=str),
        "esd_detail": json.dumps(esd_detail, default=str),
    }
    return [str(v) for v in outlier_ids], diagnostics


def detect_esd_outliers(group_df: pd.DataFrame, group_name: str,
                        analytes: Mapping[str, Mapping[str, object]] = ANALYTES,
                        n_bins: int = 2, max_outliers_bins_frac: float = 0.08) -> tuple[pd.DataFrame, pd.DataFrame]:
    outlier_rows, diagnostics_rows = [], []
    for analyte, cfg in analytes.items():
        mhs_col, ref_col = str(cfg["mhs"]), str(cfg["ref"])
        complete = group_df.dropna(subset=[mhs_col, ref_col]).copy()
        ids, diag = _validated_esd_single_analyte(
            complete, ref_col, mhs_col, n_bins=n_bins,
            max_outliers_bins_frac=max_outliers_bins_frac
        )
        diagnostics_rows.append({"analysis_group": group_name, "analyte": analyte,
                                 "mhs_column": mhs_col, "reference_column": ref_col, **diag})
        if ids:
            selected = group_df[group_df["batch_id"].astype(str).isin(ids)].copy()
            for _, row in selected.iterrows():
                ref = row.get(ref_col, np.nan)
                mhs = row.get(mhs_col, np.nan)
                residual = mhs - ref if pd.notna(mhs) and pd.notna(ref) else np.nan
                residual_pct = 100 * residual / ref if pd.notna(residual) and pd.notna(ref) and ref != 0 else np.nan
                outlier_rows.append({
                    "analysis_group": group_name,
                    "analyte": analyte,
                    "batch_id": str(row.get("batch_id", "")),
                    "bloodSampleId": row.get("bloodSampleId", ""),
                    "donor": row.get("donor", ""),
                    "specimen_type": row.get("specimen_type", ""),
                    "replicate_number": row.get("replicate_number", np.nan),
                    "mhs_column": mhs_col,
                    "reference_column": ref_col,
                    "mhs_value": mhs,
                    "reference_value": ref,
                    "residual_native": residual,
                    "residual_percent": residual_pct,
                    "exclusion_reason": "Validated generalized ESD; linked MHS/reference replicate removed for this analyte only",
                })
    outlier_columns = [
        "analysis_group", "analyte", "batch_id", "bloodSampleId", "donor",
        "specimen_type", "replicate_number", "mhs_column", "reference_column",
        "mhs_value", "reference_value", "residual_native", "residual_percent",
        "exclusion_reason",
    ]
    return pd.DataFrame(outlier_rows, columns=outlier_columns), pd.DataFrame(diagnostics_rows)


def build_donor_means(group_df: pd.DataFrame, outliers_df: pd.DataFrame,
                      analytes: Mapping[str, Mapping[str, object]] = ANALYTES) -> tuple[pd.DataFrame, pd.DataFrame]:
    donors = sorted(group_df["donor"].dropna().astype(str).unique())
    result = pd.DataFrame({"Donor": donors})
    total_counts = group_df.groupby("donor").size().rename("Replicates").reset_index().rename(columns={"donor": "Donor"})
    result = result.merge(total_counts, on="Donor", how="left")
    counts = pd.DataFrame({"Donor": donors})

    for analyte, cfg in analytes.items():
        mhs_col, ref_col = str(cfg["mhs"]), str(cfg["ref"])
        out_ids: set[str] = set()
        if not outliers_df.empty:
            out_ids = set(outliers_df.loc[outliers_df["analyte"] == analyte, "batch_id"].astype(str))
        sub = group_df[~group_df["batch_id"].astype(str).isin(out_ids)].dropna(subset=[mhs_col, ref_col]).copy()
        grouped = sub.groupby("donor")[[ref_col, mhs_col]].mean().reset_index().rename(columns={"donor": "Donor"})
        n = sub.groupby("donor").size().rename(f"{analyte}_n").reset_index().rename(columns={"donor": "Donor"})
        result = result.merge(grouped, on="Donor", how="left")
        counts = counts.merge(n, on="Donor", how="left")

    result.insert(0, "#", np.arange(1, len(result) + 1))
    ordered = ["#", "Donor"]
    for analyte, cfg in analytes.items():
        ordered.extend([str(cfg["ref"]), str(cfg["mhs"])])
    ordered.append("Replicates")
    result = result[ordered]
    return result, counts


def build_raw_with_donor_means(group_df: pd.DataFrame, donor_means: pd.DataFrame,
                               analytes: Mapping[str, Mapping[str, object]] = ANALYTES) -> pd.DataFrame:
    """Interleave clean raw rows and the final analyte-specific donor mean row."""
    cols = ["batch_id", "bloodSampleId", "specimen_type", "donor"]
    for cfg in analytes.values():
        cols.extend([str(cfg["ref"]), str(cfg["mhs"])])
    cols.extend(["global_flag", "Replicates", "row_type"])
    blocks = []
    for donor in donor_means["Donor"].astype(str):
        raw = group_df[group_df["donor"].astype(str) == donor].copy()
        raw["Replicates"] = np.nan
        raw["row_type"] = "raw replicate"
        raw = raw.reindex(columns=cols)
        mean_row = donor_means[donor_means["Donor"].astype(str) == donor].iloc[0]
        row = {"batch_id": "", "bloodSampleId": "", "specimen_type": "",
               "donor": f"{donor}_aver", "global_flag": "",
               "Replicates": mean_row["Replicates"], "row_type": "donor mean after analyte-specific ESD"}
        for cfg in analytes.values():
            row[str(cfg["ref"])] = mean_row[str(cfg["ref"])]
            row[str(cfg["mhs"])] = mean_row[str(cfg["mhs"])]
        blocks.append(raw)
        blocks.append(pd.DataFrame([row], columns=cols))
    return pd.concat(blocks, ignore_index=True) if blocks else pd.DataFrame(columns=cols)


def huber_regression(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    model = HuberRegressor()
    model.fit(np.asarray(x).reshape(-1, 1), np.asarray(y))
    return float(model.coef_[0]), float(model.intercept_)


def pearson_ci(r: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(r):
        return np.nan, np.nan
    r = float(np.clip(r, -0.999999999, 0.999999999))
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - z_crit * se)), float(np.tanh(z + z_crit * se))


def bootstrap_huber_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 5000,
                       alpha: float = 0.05, random_state: int = 20260722) -> tuple[float, float, float, float, int]:
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 4:
        return np.nan, np.nan, np.nan, np.nan, 0
    rng = np.random.default_rng(random_state)
    slopes, intercepts = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        xb, yb = x[idx], y[idx]
        if len(np.unique(xb)) < 2:
            continue
        try:
            slope, intercept = huber_regression(xb, yb)
            slopes.append(slope)
            intercepts.append(intercept)
        except Exception:
            continue
    if not slopes:
        return np.nan, np.nan, np.nan, np.nan, 0
    q = [100 * alpha / 2, 100 * (1 - alpha / 2)]
    slo, shi = np.percentile(slopes, q)
    ilo, ihi = np.percentile(intercepts, q)
    return float(slo), float(shi), float(ilo), float(ihi), len(slopes)


def compute_regression_metrics(donor_means: pd.DataFrame, group_name: str,
                               ac_dict: Mapping[str, Mapping[str, float]] | None = None,
                               analytes: Mapping[str, Mapping[str, object]] = ANALYTES,
                               run_huber_bootstrap: bool = True,
                               n_boot: int = 5000) -> pd.DataFrame:
    rows = []
    ac_dict = ac_dict or {}
    for analyte, cfg in analytes.items():
        mhs_col, ref_col = str(cfg["mhs"]), str(cfg["ref"])
        sub = donor_means[["Donor", ref_col, mhs_col]].dropna()
        n = len(sub)
        if n >= 2 and sub[ref_col].nunique() >= 2 and sub[mhs_col].nunique() >= 2:
            x, y = sub[ref_col].to_numpy(float), sub[mhs_col].to_numpy(float)
            r, p = pearsonr(x, y)
            ci_lo, ci_hi = pearson_ci(r, n)
            slope, intercept = huber_regression(x, y)
            r2 = r2_score(x, y)
            rmse = math.sqrt(mean_squared_error(x, y))
            bias = float(np.mean(y - x))
            bias_pct = float(100 * bias / np.mean(x)) if np.mean(x) != 0 else np.nan
            if run_huber_bootstrap:
                bslo, bshi, bilo, bihi, n_success = bootstrap_huber_ci(x, y, n_boot=n_boot)
            else:
                bslo = bshi = bilo = bihi = np.nan
                n_success = 0
        else:
            r = p = ci_lo = ci_hi = slope = intercept = r2 = rmse = bias = bias_pct = np.nan
            bslo = bshi = bilo = bihi = np.nan
            n_success = 0

        normal_min, normal_max = cfg["normal"]
        ac = ac_dict.get(str(cfg["mhs"]), ac_dict.get(analyte, {}))
        ac_corr = ac.get("corr", np.nan)
        ac_slope_low = ac.get("slope_low", np.nan)
        ac_slope_high = ac.get("slope_high", np.nan)
        ac_intercept = ac.get("intercept_limit", np.nan)
        bias_limit = ac.get("bias_limit", np.nan)
        rows.append({
            "analysis_group": group_name,
            "analyte": analyte,
            "unit": cfg["unit"],
            "n_donors": n,
            "pearson_r": r,
            "pearson_p": p,
            "pearson_ci_95_low_fisher_z": ci_lo,
            "pearson_ci_95_high_fisher_z": ci_hi,
            "huber_slope": slope,
            "huber_intercept": intercept,
            "huber_equation": f"y = {slope:.4g}x {intercept:+.4g}" if np.isfinite(slope) else "",
            "huber_slope_boot_ci_95_low": bslo,
            "huber_slope_boot_ci_95_high": bshi,
            "huber_intercept_boot_ci_95_low": bilo,
            "huber_intercept_boot_ci_95_high": bihi,
            "huber_bootstrap_successful_resamples": n_success,
            "r_squared_vs_identity": r2,
            "rmse_native": rmse,
            "mean_bias_native": bias,
            "mean_bias_percent": bias_pct,
            "reference_min": sub[ref_col].min() if n else np.nan,
            "reference_max": sub[ref_col].max() if n else np.nan,
            "mhs_min": sub[mhs_col].min() if n else np.nan,
            "mhs_max": sub[mhs_col].max() if n else np.nan,
            "normal_min": normal_min,
            "normal_max": normal_max,
            "reference_below_normal_n": int((sub[ref_col] < normal_min).sum()) if n else 0,
            "reference_within_normal_n": int(((sub[ref_col] >= normal_min) & (sub[ref_col] <= normal_max)).sum()) if n else 0,
            "reference_above_normal_n": int((sub[ref_col] > normal_max).sum()) if n else 0,
            "AC_corr": ac_corr,
            "AC_slope_low": ac_slope_low,
            "AC_slope_high": ac_slope_high,
            "AC_intercept_limit": ac_intercept,
            "AC_bias_limit_percent": bias_limit,
            "pass_pearson_r": bool(r >= ac_corr) if np.isfinite(r) and np.isfinite(ac_corr) else np.nan,
            "pass_fisher_ci_low": bool(ci_lo >= ac_corr) if np.isfinite(ci_lo) and np.isfinite(ac_corr) else np.nan,
            "pass_huber_slope": bool(ac_slope_low <= slope <= ac_slope_high) if np.isfinite(slope) and np.isfinite(ac_slope_low) and np.isfinite(ac_slope_high) else np.nan,
            "pass_huber_intercept": bool(abs(intercept) <= ac_intercept) if np.isfinite(intercept) and np.isfinite(ac_intercept) else np.nan,
            "pass_bias": bool(abs(bias_pct) <= bias_limit) if np.isfinite(bias_pct) and np.isfinite(bias_limit) else np.nan,
        })
    return pd.DataFrame(rows)


REGRESSION_ANNOTATION_POSITIONS: dict[str, tuple[float, float, str, str]] = {
    "top_left": (0.025, 0.975, "left", "top"),
    "top_center": (0.500, 0.975, "center", "top"),
    "top_right": (0.975, 0.975, "right", "top"),
    "middle_left": (0.025, 0.500, "left", "center"),
    "middle_right": (0.975, 0.500, "right", "center"),
    "bottom_left": (0.025, 0.025, "left", "bottom"),
    "bottom_center": (0.500, 0.025, "center", "bottom"),
    "bottom_right": (0.975, 0.025, "right", "bottom"),
}


def _annotation_rect(position: str, width: float, height: float) -> tuple[float, float, float, float]:
    """Approximate an annotation box in axes-fraction coordinates."""
    x, y, ha, va = REGRESSION_ANNOTATION_POSITIONS[position]
    if ha == "left":
        x0, x1 = x, x + width
    elif ha == "right":
        x0, x1 = x - width, x
    else:
        x0, x1 = x - width / 2.0, x + width / 2.0
    if va == "top":
        y0, y1 = y - height, y
    elif va == "bottom":
        y0, y1 = y, y + height
    else:
        y0, y1 = y - height / 2.0, y + height / 2.0
    return x0, x1, y0, y1


def _rect_overlap(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> float:
    x_overlap = max(0.0, min(a[1], b[1]) - max(a[0], b[0]))
    y_overlap = max(0.0, min(a[3], b[3]) - max(a[2], b[2]))
    return x_overlap * y_overlap


def _score_annotation_rect(
    rect: tuple[float, float, float, float],
    point_axes: np.ndarray,
    red_line_axes_x: Sequence[float],
    trend_axes: np.ndarray | None = None,
) -> float:
    """Penalize boxes crossing observations, normal-range lines or the fitted line."""
    x0, x1, y0, y1 = rect
    score = 0.0
    if point_axes.size:
        inside = (
            (point_axes[:, 0] >= x0) & (point_axes[:, 0] <= x1) &
            (point_axes[:, 1] >= y0) & (point_axes[:, 1] <= y1)
        )
        score += float(inside.sum()) * 125.0
        # Extra clearance around observations so labels do not visually touch markers.
        padded = (
            (point_axes[:, 0] >= x0 - 0.025) & (point_axes[:, 0] <= x1 + 0.025) &
            (point_axes[:, 1] >= y0 - 0.025) & (point_axes[:, 1] <= y1 + 0.025)
        )
        score += float((padded & ~inside).sum()) * 18.0
    for red_x in red_line_axes_x:
        if np.isfinite(red_x) and x0 - 0.012 <= red_x <= x1 + 0.012:
            score += 85.0
    if trend_axes is not None and trend_axes.size:
        inside_trend = (
            (trend_axes[:, 0] >= x0) & (trend_axes[:, 0] <= x1) &
            (trend_axes[:, 1] >= y0) & (trend_axes[:, 1] <= y1)
        )
        score += float(inside_trend.sum()) * 1.2
    # Discourage clipping at the axes edge.
    if x0 < 0 or x1 > 1 or y0 < 0 or y1 > 1:
        score += 10_000.0
    return score


def _choose_regression_annotation_positions(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    normal_min: float,
    normal_max: float,
    x_line: np.ndarray | None,
    y_line: np.ndarray | None,
    equation_position: str = "auto",
    metrics_position: str = "auto",
) -> tuple[str, str]:
    """Choose reproducible empty-space positions for equation and metrics labels."""
    valid = set(REGRESSION_ANNOTATION_POSITIONS)
    if equation_position != "auto" and equation_position not in valid:
        raise ValueError(f"Unknown equation annotation position: {equation_position}")
    if metrics_position != "auto" and metrics_position not in valid:
        raise ValueError(f"Unknown metrics annotation position: {metrics_position}")

    # Ensure axes limits and transforms are final before converting to axes fractions.
    ax.relim()
    ax.autoscale_view()
    if len(x):
        point_display = ax.transData.transform(np.column_stack([x, y]))
        point_axes = ax.transAxes.inverted().transform(point_display)
    else:
        point_axes = np.empty((0, 2), dtype=float)

    y_mid = float(np.nanmean(y)) if len(y) and np.isfinite(np.nanmean(y)) else 0.0
    red_display = ax.transData.transform(np.array([[normal_min, y_mid], [normal_max, y_mid]], dtype=float))
    red_axes_x = ax.transAxes.inverted().transform(red_display)[:, 0]

    trend_axes = None
    if x_line is not None and y_line is not None and len(x_line):
        trend_display = ax.transData.transform(np.column_stack([x_line, y_line]))
        trend_axes = ax.transAxes.inverted().transform(trend_display)

    eq_candidates = ["top_left", "bottom_left", "top_center", "bottom_center", "middle_left", "top_right", "bottom_right", "middle_right"]
    metric_candidates = ["top_right", "bottom_right", "top_left", "bottom_left", "middle_right", "top_center", "bottom_center", "middle_left"]
    eq_rects = {p: _annotation_rect(p, width=0.34, height=0.085) for p in eq_candidates}
    metric_rects = {p: _annotation_rect(p, width=0.38, height=0.245) for p in metric_candidates}

    eq_choices = [equation_position] if equation_position != "auto" else eq_candidates
    metric_choices = [metrics_position] if metrics_position != "auto" else metric_candidates
    best: tuple[float, str, str] | None = None
    for eq_rank, eq_pos in enumerate(eq_choices):
        for metric_rank, metric_pos in enumerate(metric_choices):
            if eq_pos == metric_pos:
                continue
            eq_rect, metric_rect = eq_rects[eq_pos], metric_rects[metric_pos]
            score = _score_annotation_rect(eq_rect, point_axes, red_axes_x, trend_axes)
            score += _score_annotation_rect(metric_rect, point_axes, red_axes_x, trend_axes)
            score += _rect_overlap(eq_rect, metric_rect) * 20_000.0
            # Stable tie-breaker retains familiar top-left equation/top-right metrics when clear.
            score += eq_rank * 0.05 + metric_rank * 0.05
            candidate = (score, eq_pos, metric_pos)
            if best is None or candidate < best:
                best = candidate
    if best is None:
        return "top_left", "top_right"
    return best[1], best[2]


def plot_regression_analyte(donor_means: pd.DataFrame, metrics_row: pd.Series,
                            analytes: Mapping[str, Mapping[str, object]] = ANALYTES,
                            output_path: str | Path | None = None,
                            ax=None,
                            equation_position: str = "auto",
                            metrics_position: str = "auto"):
    analyte = str(metrics_row["analyte"])
    cfg = analytes[analyte]
    ref_col, mhs_col = str(cfg["ref"]), str(cfg["mhs"])
    sub = donor_means[["Donor", ref_col, mhs_col]].dropna()
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
    else:
        fig = ax.figure
    x, y = sub[ref_col].to_numpy(float), sub[mhs_col].to_numpy(float)
    ax.scatter(x, y, s=34, alpha=0.9, edgecolors="black", linewidths=0.35, color="#1F77B4", zorder=3)
    x_line = None
    y_line = None
    if len(x) >= 2 and np.isfinite(metrics_row["huber_slope"]):
        x_line = np.linspace(min(x.min(), metrics_row["normal_min"]), max(x.max(), metrics_row["normal_max"]), 200)
        y_line = metrics_row["huber_slope"] * x_line + metrics_row["huber_intercept"]
        ax.plot(x_line, y_line, linestyle=":", linewidth=1.6, color="#1F77B4", zorder=2)
    ax.axvline(metrics_row["normal_min"], linewidth=1.6, color="red", zorder=1)
    ax.axvline(metrics_row["normal_max"], linewidth=1.6, color="red", zorder=1)
    ax.set_title(analyte, fontweight="bold")
    ax.set_xlabel("Sysmex reference")
    ax.set_ylabel("MHS")
    ax.grid(True, alpha=0.28, zorder=0)

    eq_pos, metric_pos = _choose_regression_annotation_positions(
        ax, x, y,
        float(metrics_row["normal_min"]), float(metrics_row["normal_max"]),
        x_line, y_line,
        equation_position=equation_position,
        metrics_position=metrics_position,
    )
    annotation_box = {
        "boxstyle": "round,pad=0.24",
        "facecolor": "white",
        "edgecolor": "#777777",
        "linewidth": 0.45,
        "alpha": 0.88,
    }
    eq_x, eq_y, eq_ha, eq_va = REGRESSION_ANNOTATION_POSITIONS[eq_pos]
    metric_x, metric_y, metric_ha, metric_va = REGRESSION_ANNOTATION_POSITIONS[metric_pos]
    eq = metrics_row["huber_equation"]
    ax.text(eq_x, eq_y, eq, transform=ax.transAxes, ha=eq_ha, va=eq_va,
            fontsize=8, bbox=annotation_box, zorder=6)
    ax.text(metric_x, metric_y,
            f"R = {metrics_row['pearson_r']:.3f}\n95% CI = {metrics_row['pearson_ci_95_low_fisher_z']:.3f} to {metrics_row['pearson_ci_95_high_fisher_z']:.3f}\n"
            f"n = {int(metrics_row['n_donors'])}\nRef min/max = {metrics_row['reference_min']:.3g}/{metrics_row['reference_max']:.3g}",
            transform=ax.transAxes, ha=metric_ha, va=metric_va, fontsize=8,
            bbox=annotation_box, zorder=6)
    if own_fig:
        fig.tight_layout()
        if output_path:
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=220, bbox_inches="tight")
        plt.close(fig)
    return ax


def save_regression_plots(donor_means: pd.DataFrame, metrics: pd.DataFrame,
                          group_name: str, output_dir: str | Path,
                          equation_position: str = "auto",
                          metrics_position: str = "auto") -> dict[str, Path]:
    output_dir = Path(output_dir)
    individual_dir = output_dir / "plots_individual"
    individual_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for _, row in metrics.iterrows():
        path = individual_dir / f"{group_name}_{row['analyte']}_regression.png"
        plot_regression_analyte(
            donor_means, row, output_path=path,
            equation_position=equation_position,
            metrics_position=metrics_position,
        )
        paths[str(row["analyte"])] = path

    fig, axes = plt.subplots(4, 3, figsize=(16, 18))
    for ax, (_, row) in zip(axes.ravel(), metrics.iterrows()):
        plot_regression_analyte(
            donor_means, row, ax=ax,
            equation_position=equation_position,
            metrics_position=metrics_position,
        )
    # Hide unused axes when a subset of analytes is selected.
    for ax in axes.ravel()[len(metrics):]:
        ax.set_visible(False)
    fig.suptitle(f"{group_name}: MHS versus matched Sysmex reference after replicate-level filtering", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    panel_path = output_dir / f"{group_name}_regression_panel.png"
    fig.savefig(panel_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths["panel"] = panel_path
    return paths


def build_cross_specimen_table(group_a: Mapping[str, object], group_b: Mapping[str, object],
                               name_a: str = "Capillary", name_b: str = "MK2",
                               analytes: Mapping[str, Mapping[str, object]] = ANALYTES) -> pd.DataFrame:
    a = group_a["donor_means"].copy()
    b = group_b["donor_means"].copy()
    a_cols = {"Donor": "Donor"}
    b_cols = {"Donor": "Donor"}
    for analyte, cfg in analytes.items():
        a_cols[str(cfg["ref"])] = f"{analyte}_ref_{name_a}"
        a_cols[str(cfg["mhs"])] = f"{analyte}_MHS_{name_a}"
        b_cols[str(cfg["ref"])] = f"{analyte}_ref_{name_b}"
        b_cols[str(cfg["mhs"])] = f"{analyte}_MHS_{name_b}"
    return a[list(a_cols)].rename(columns=a_cols).merge(
        b[list(b_cols)].rename(columns=b_cols), on="Donor", how="inner"
    )


def compute_cross_specimen_metrics(matched: pd.DataFrame, name_a: str = "Capillary", name_b: str = "MK2",
                                  analytes: Mapping[str, Mapping[str, object]] = ANALYTES) -> pd.DataFrame:
    rows = []
    for analyte in analytes:
        for source in ["MHS", "ref"]:
            xcol = f"{analyte}_{source}_{name_b}"
            ycol = f"{analyte}_{source}_{name_a}"
            sub = matched[[xcol, ycol]].dropna()
            n = len(sub)
            if n >= 2 and sub[xcol].nunique() >= 2 and sub[ycol].nunique() >= 2:
                x, y = sub[xcol].to_numpy(float), sub[ycol].to_numpy(float)
                r, p = pearsonr(x, y)
                lo, hi = pearson_ci(r, n)
                slope, intercept = huber_regression(x, y)
            else:
                r = p = lo = hi = slope = intercept = np.nan
            rows.append({"analyte": analyte, "source": "MHS" if source == "MHS" else "Sysmex reference",
                         "x_specimen": name_b, "y_specimen": name_a, "n_matched_donors": n,
                         "pearson_r": r, "pearson_p": p, "pearson_ci_95_low_fisher_z": lo,
                         "pearson_ci_95_high_fisher_z": hi, "huber_slope": slope,
                         "huber_intercept": intercept})
    return pd.DataFrame(rows)


def export_excel_xlsxwriter(output_path: str | Path, group_outputs: Mapping[str, Mapping[str, object]],
                            global_exclusions: pd.DataFrame, config_df: pd.DataFrame) -> None:
    """Create a portable Excel workbook in the same broad layout as 'for Abi.xlsx'."""
    output_path = Path(output_path)
    with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
        workbook = writer.book
        header_fmt = workbook.add_format({"bold": True, "font_color": "white", "bg_color": "#1F4E78", "border": 1, "align": "center"})
        subheader_fmt = workbook.add_format({"bold": True, "bg_color": "#D9EAD3", "border": 1, "align": "center"})
        cell_fmt = workbook.add_format({"border": 1, "num_format": "0.000"})
        pass_fmt = workbook.add_format({"bg_color": "#D9EAD3", "border": 1})
        fail_fmt = workbook.add_format({"bg_color": "#F4CCCC", "border": 1})

        all_metrics = pd.concat([v["metrics"] for v in group_outputs.values()], ignore_index=True) if group_outputs else pd.DataFrame()
        all_metrics.to_excel(writer, sheet_name="Summary_All", index=False)
        ws = writer.sheets["Summary_All"]
        ws.freeze_panes(1, 0)
        ws.set_row(0, 30, header_fmt)
        ws.set_column(0, 0, 23)
        ws.set_column(1, 1, 12)
        ws.set_column(2, 2, 18)
        ws.set_column(3, 3, 11)
        ws.set_column(4, 10, 22)
        ws.set_column(11, len(all_metrics.columns) - 1, 24)
        for col_name in ["pass_pearson_r", "pass_fisher_ci_low", "pass_huber_slope", "pass_huber_intercept", "pass_bias"]:
            if col_name in all_metrics.columns:
                col = all_metrics.columns.get_loc(col_name)
                ws.conditional_format(1, col, max(1, len(all_metrics)), col, {"type": "cell", "criteria": "==", "value": True, "format": pass_fmt})
                ws.conditional_format(1, col, max(1, len(all_metrics)), col, {"type": "cell", "criteria": "==", "value": False, "format": fail_fmt})

        for group_name, bundle in group_outputs.items():
            short = re.sub(r"[^A-Za-z0-9]", "", group_name)[:18]
            donor = bundle["donor_means"]
            metrics = bundle["metrics"]
            raw = bundle["raw_clean"]
            raw_with_means = bundle["raw_with_means"]
            outliers = bundle["outliers"]
            diagnostics = bundle["diagnostics"]
            counts = bundle["replicate_counts"]
            plots = bundle["plots"]

            donor_sheet = f"Reg_{short}"[:31]
            donor.to_excel(writer, sheet_name=donor_sheet, index=False, startrow=0)
            ws_d = writer.sheets[donor_sheet]
            ws_d.freeze_panes(1, 2)
            ws_d.set_row(0, 24, header_fmt)
            ws_d.set_column(0, 0, 5)
            ws_d.set_column(1, 1, 12)
            ws_d.set_column(2, len(donor.columns) - 2, 13, cell_fmt)
            ws_d.set_column(len(donor.columns) - 1, len(donor.columns) - 1, 11)
            start = len(donor) + 3
            metric_display = metrics[["analyte", "n_donors", "pearson_r", "pearson_ci_95_low_fisher_z", "pearson_ci_95_high_fisher_z",
                                      "huber_slope", "huber_intercept", "huber_equation", "reference_min", "reference_max",
                                      "mhs_min", "mhs_max", "normal_min", "normal_max"]]
            metric_display.to_excel(writer, sheet_name=donor_sheet, index=False, startrow=start)
            ws_d.set_row(start, 24, subheader_fmt)
            ws_d.set_column(0, 0, 12)
            ws_d.set_column(1, 14, 19)
            plot_start_row = start + len(metric_display) + 3
            for i, analyte in enumerate(ANALYTES):
                img = plots.get(analyte)
                if img and Path(img).exists():
                    row = plot_start_row + (i // 3) * 25
                    col = (i % 3) * 9
                    ws_d.insert_image(row, col, str(img), {"x_scale": 0.72, "y_scale": 0.72})

            metrics.to_excel(writer, sheet_name=f"Metrics_{short}"[:31], index=False)
            raw.to_excel(writer, sheet_name=f"Raw_{short}"[:31], index=False)
            raw_with_means.to_excel(writer, sheet_name=f"RawMeans_{short}"[:31], index=False)
            outliers.to_excel(writer, sheet_name=f"Outliers_{short}"[:31], index=False)
            diagnostics.to_excel(writer, sheet_name=f"ESDdiag_{short}"[:31], index=False)
            counts.to_excel(writer, sheet_name=f"RepCounts_{short}"[:31], index=False)
            for sname in [f"Metrics_{short}"[:31], f"Raw_{short}"[:31], f"RawMeans_{short}"[:31], f"Outliers_{short}"[:31], f"ESDdiag_{short}"[:31], f"RepCounts_{short}"[:31]]:
                sheet = writer.sheets[sname]
                sheet.freeze_panes(1, 0)
                sheet.set_row(0, 24, header_fmt)
                sheet.set_column(0, 80, 18)

        global_exclusions.to_excel(writer, sheet_name="Global_Exclusions", index=False)
        config_df.to_excel(writer, sheet_name="Configuration", index=False)
        for sname in ["Global_Exclusions", "Configuration"]:
            ws2 = writer.sheets[sname]
            ws2.freeze_panes(1, 0)
            ws2.set_row(0, 24, header_fmt)
            ws2.set_column(0, 80, 18)


def run_pipeline(input_csv: str | Path, output_root: str | Path,
                 analysis_groups: Mapping[str, Sequence[str]] = DEFAULT_ANALYSIS_GROUPS,
                 groups_to_run: Sequence[str] = ("Capillary_Mc_Mc2", "Venous_MK2"),
                 ac_json: str | Path | None = None,
                 run_huber_bootstrap: bool = True,
                 n_boot: int = 5000,
                 regression_equation_position: str = "auto",
                 regression_metrics_position: str = "auto") -> dict[str, object]:
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    df = load_and_prepare(input_csv)

    invalid_parse = df[~df["parse_ok"]].copy()
    global_exclusions = df[df["global_flag_bool"] != False].copy()
    global_exclusions["exclusion_reason"] = np.where(
        global_exclusions["global_flag_bool"] == True,
        "global_flag = TRUE",
        "global_flag missing or unparseable"
    )
    global_exclusions.to_csv(output_root / "01_global_flag_exclusions.csv", index=False)
    invalid_parse.to_csv(output_root / "00_unparsed_sample_ids.csv", index=False)

    clean = df[df["global_flag_bool"] == False].copy()
    clean.to_csv(output_root / "02_all_global_flag_false_rows.csv", index=False)

    ac_dict = {}
    if ac_json:
        with open(ac_json, "r", encoding="utf-8") as f:
            ac_dict = json.load(f)

    group_outputs: dict[str, dict[str, object]] = {}
    all_outliers = []
    all_metrics = []
    for group_name in groups_to_run:
        specimen_types = list(analysis_groups[group_name])
        group_df = clean[clean["specimen_type"].isin(specimen_types)].copy()
        if group_df.empty:
            print(f"Skipping {group_name}: no rows for specimen types {specimen_types}")
            continue
        group_dir = output_root / group_name
        group_dir.mkdir(parents=True, exist_ok=True)
        group_df.to_csv(group_dir / "raw_clean_replicates.csv", index=False)
        outliers, diagnostics = detect_esd_outliers(group_df, group_name)
        outliers.to_csv(group_dir / "esd_outliers_removed.csv", index=False)
        diagnostics.to_csv(group_dir / "esd_diagnostics.csv", index=False)
        donor_means, replicate_counts = build_donor_means(group_df, outliers)
        donor_means.to_csv(group_dir / "donor_means_after_outlier_removal.csv", index=False)
        raw_with_means = build_raw_with_donor_means(group_df, donor_means)
        raw_with_means.to_csv(group_dir / "raw_data_with_donor_means.csv", index=False)
        replicate_counts.to_csv(group_dir / "donor_analyte_replicate_counts.csv", index=False)
        metrics = compute_regression_metrics(
            donor_means, group_name, ac_dict=ac_dict,
            run_huber_bootstrap=run_huber_bootstrap, n_boot=n_boot
        )
        metrics.to_csv(group_dir / "regression_metrics_huber_fisherz.csv", index=False)
        plots = save_regression_plots(
            donor_means, metrics, group_name, group_dir,
            equation_position=regression_equation_position,
            metrics_position=regression_metrics_position,
        )
        group_outputs[group_name] = {
            "raw_clean": group_df,
            "raw_with_means": raw_with_means,
            "outliers": outliers,
            "diagnostics": diagnostics,
            "donor_means": donor_means,
            "replicate_counts": replicate_counts,
            "metrics": metrics,
            "plots": plots,
            "specimen_types": specimen_types,
        }
        all_outliers.append(outliers)
        all_metrics.append(metrics)

    if all_outliers:
        pd.concat(all_outliers, ignore_index=True).to_csv(output_root / "03_all_esd_outliers_removed.csv", index=False)
    if all_metrics:
        pd.concat(all_metrics, ignore_index=True).to_csv(output_root / "04_all_regression_metrics.csv", index=False)

    if "Capillary_Mc_Mc2" in group_outputs and "Venous_MK2" in group_outputs:
        cross_dir = output_root / "CrossSpecimen_Capillary_vs_MK2"
        cross_dir.mkdir(parents=True, exist_ok=True)
        matched = build_cross_specimen_table(
            group_outputs["Capillary_Mc_Mc2"], group_outputs["Venous_MK2"],
            name_a="Capillary", name_b="MK2"
        )
        matched.to_csv(cross_dir / "matched_donor_means.csv", index=False)
        cross_metrics = compute_cross_specimen_metrics(matched, name_a="Capillary", name_b="MK2")
        cross_metrics.to_csv(cross_dir / "cross_specimen_regression_metrics.csv", index=False)
    else:
        matched = pd.DataFrame()
        cross_metrics = pd.DataFrame()

    config_rows = []
    for analyte, cfg in ANALYTES.items():
        config_rows.append({"section": "analyte", "name": analyte, "mhs_column": cfg["mhs"], "reference_column": cfg["ref"],
                            "normal_min": cfg["normal"][0], "normal_max": cfg["normal"][1], "unit": cfg["unit"]})
    for name, types in analysis_groups.items():
        config_rows.append({"section": "analysis_group", "name": name, "mhs_column": ", ".join(types),
                            "reference_column": "", "normal_min": np.nan, "normal_max": np.nan, "unit": ""})
    config_df = pd.DataFrame(config_rows)
    config_df.to_csv(output_root / "05_configuration.csv", index=False)

    excel_path = output_root / "PROXIMA_Trueness_Regression_Output.xlsx"
    export_excel_xlsxwriter(excel_path, group_outputs, global_exclusions, config_df)
    return {"input": df, "clean": clean, "global_exclusions": global_exclusions,
            "groups": group_outputs, "cross_specimen_matched": matched,
            "cross_specimen_metrics": cross_metrics, "excel_path": excel_path}



# =============================================================================
# PROXIMA paired-specimen Bland–Altman extension
# =============================================================================
from collections import defaultdict
from dataclasses import dataclass
from statistics import mean
from scipy.stats import shapiro

# Contextual acceptance criteria reproduced from the supplied Algocyte report.
# IMPORTANT: these are trueness benchmarks, not universally established
# capillary-versus-venous acceptance criteria. Review/edit before final use.
AC_LIMIT_PCT = {
    "RBC": 7.0, "WBC": 12.5, "PLT": 17.5, "HCT": 8.6, "HGB": 6.0,
    "MCV": 6.6, "RDW": 13.5, "MCH": 9.5, "MCHC": 6.8,
    "NEUT": 11.6, "LYMPH": 14.0, "MXD": 25.0,
}
CLIA_LIMIT_PCT = {"RBC": 4.0, "WBC": 10.0, "PLT": 25.0, "HCT": 4.0, "HGB": 4.0}

SHAPIRO_ALPHA = 0.05
OUTLIER_ALPHA = 0.05
OUTLIER_MAD_Z = 3.5
BA_BOOTSTRAP_ITERATIONS = 10_000
DONOR_PROFILE_BOOTSTRAPS = 5_000
BA_RANDOM_SEED = 20260723

METHOD_DESCRIPTIONS = {
    "M02": "Direct signed capillary-versus-venous MHS comparison",
    "M11": "Signed Sysmex reference-adjusted capillary-versus-venous comparison (preferred reference-adjusted option)",
    "M05": "Absolute-gap magnitude sensitivity analysis; not signed and not a substitute for M02/M11",
}


def natural_donor_key(value):
    text = str(value)
    match = re.fullmatch(r"([A-Za-z]*)(\d+)([A-Za-z]*)", text)
    if not match:
        return (text.lower(), 0, "")
    prefix, number, suffix = match.groups()
    return (prefix.lower(), int(number), suffix.lower())


def grubbs_critical(n: int, alpha: float = OUTLIER_ALPHA) -> float:
    """Two-sided Grubbs G critical value computed manually for the observed n."""
    if n < 3:
        return np.nan
    tcrit = float(stats.t.ppf(1.0 - alpha / (2.0 * n), n - 2))
    return float(((n - 1) / np.sqrt(n)) * np.sqrt(tcrit**2 / (n - 2 + tcrit**2)))


def detect_one_grubbs_or_mad_outlier(residuals: Sequence[float]) -> dict[str, object]:
    """At most one outlier: Shapiro -> Grubbs if normal, MAD modified-Z if non-normal."""
    values = np.asarray(residuals, dtype=float)
    finite_idx = np.flatnonzero(np.isfinite(values))
    v = values[finite_idx]
    n = len(v)
    result: dict[str, object] = {
        "n_residuals": n,
        "shapiro_W": np.nan,
        "shapiro_p": np.nan,
        "normality_branch": "insufficient",
        "method": "none",
        "candidate_index": None,
        "statistic": np.nan,
        "critical_threshold": np.nan,
        "removed": False,
        "residual_median": float(np.median(v)) if n else np.nan,
        "residual_MAD": np.nan,
        "residual_mean": float(np.mean(v)) if n else np.nan,
        "residual_SD": float(np.std(v, ddof=1)) if n >= 2 else np.nan,
    }
    if n < 3 or not np.isfinite(v).all() or np.std(v, ddof=1) <= 0:
        return result
    sw = shapiro(v)
    result["shapiro_W"] = float(sw.statistic)
    result["shapiro_p"] = float(sw.pvalue)
    if sw.pvalue >= SHAPIRO_ALPHA:
        mu = float(np.mean(v)); sd = float(np.std(v, ddof=1))
        scores = np.abs(v - mu) / sd
        j = int(np.argmax(scores))
        g = float(scores[j]); gcrit = grubbs_critical(n, OUTLIER_ALPHA)
        result.update({
            "normality_branch": "normal", "method": "Grubbs",
            "candidate_index": int(finite_idx[j]), "statistic": g,
            "critical_threshold": gcrit, "removed": bool(g >= gcrit),
        })
    else:
        med = float(np.median(v)); mad = float(np.median(np.abs(v - med)))
        result.update({"normality_branch": "non-normal", "method": "MAD modified Z",
                       "residual_MAD": mad, "critical_threshold": OUTLIER_MAD_Z})
        if mad > 0:
            modified_z = 0.6745 * (v - med) / mad
            scores = np.abs(modified_z)
            j = int(np.argmax(scores)); z = float(scores[j])
            result.update({"candidate_index": int(finite_idx[j]), "statistic": z,
                           "removed": bool(z >= OUTLIER_MAD_Z)})
    return result


def _flag_is_false(series: pd.Series) -> pd.Series:
    return series.map(normalize_bool).eq(False)


def build_ba_filtered_replicates(
    df: pd.DataFrame,
    specimen_a_types: Sequence[str],
    specimen_b_types: Sequence[str],
    specimen_a_label: str,
    specimen_b_label: str,
    analytes: Mapping[str, Mapping[str, object]] = ANALYTES,
) -> dict[str, object]:
    """Global/analyte flag filtering and one analyte-level raw-replicate outlier."""
    work = df.copy().reset_index(drop=True)
    work["_row_id"] = np.arange(len(work), dtype=int)
    global_keep = work["global_flag_bool"].eq(False) & work["parse_ok"]
    global_exclusions = work.loc[~global_keep].copy()
    global_exclusions["exclusion_reason"] = np.where(
        ~work.loc[~global_keep, "parse_ok"], "Unparsed bloodSampleId", "global_flag was not FALSE"
    )
    base = work.loc[global_keep].copy()
    allowed_types = set(specimen_a_types) | set(specimen_b_types)
    base = base[base["specimen_type"].isin(allowed_types)].copy()
    base["BA_specimen_side"] = np.where(base["specimen_type"].isin(specimen_a_types), specimen_a_label, specimen_b_label)

    filtered_by_analyte: dict[str, pd.DataFrame] = {}
    audit_rows: list[dict[str, object]] = []
    removed_rows: list[dict[str, object]] = []
    analyte_flag_exclusions: list[dict[str, object]] = []

    for analyte, cfg in analytes.items():
        mhs_col, ref_col = str(cfg["mhs"]), str(cfg["ref"])
        flag_col = f"{analyte}_flag" if f"{analyte}_flag" in base.columns else None
        candidate = base.copy()
        if flag_col:
            flag_keep = _flag_is_false(candidate[flag_col])
            excluded = candidate.loc[~flag_keep].copy()
            for _, row in excluded.iterrows():
                analyte_flag_exclusions.append({
                    "analyte": analyte, "flag_column": flag_col,
                    "batch_id": row.get("batch_id"), "bloodSampleId": row.get("bloodSampleId"),
                    "donor": row.get("donor"), "specimen_type": row.get("specimen_type"),
                    "replicate_number": row.get("replicate_number"), "flag_value": row.get(flag_col),
                    "exclusion_reason": "analyte-specific flag was not FALSE",
                })
            candidate = candidate.loc[flag_keep].copy()
        candidate = candidate.dropna(subset=[mhs_col, ref_col]).copy()
        candidate = candidate[candidate[ref_col] != 0].copy()
        candidate["residual_pct"] = 100.0 * (candidate[mhs_col] - candidate[ref_col]) / candidate[ref_col]
        decision = detect_one_grubbs_or_mad_outlier(candidate["residual_pct"].to_numpy(float))
        selected = None
        if decision["candidate_index"] is not None:
            selected = candidate.iloc[int(decision["candidate_index"])]
        removed_row_id = int(selected["_row_id"]) if selected is not None and decision["removed"] else None
        filtered = candidate[candidate["_row_id"] != removed_row_id].copy()
        filtered_by_analyte[analyte] = filtered

        audit = {
            "analyte": analyte, "mhs_column": mhs_col, "reference_column": ref_col,
            "n_valid_replicates_before": len(candidate), "n_valid_replicates_after": len(filtered),
            **decision,
            "candidate_row_id": int(selected["_row_id"]) if selected is not None else None,
            "candidate_batch_id": selected.get("batch_id") if selected is not None else None,
            "candidate_bloodSampleId": selected.get("bloodSampleId") if selected is not None else None,
            "candidate_donor": selected.get("donor") if selected is not None else None,
            "candidate_specimen_type": selected.get("specimen_type") if selected is not None else None,
            "candidate_specimen_side": selected.get("BA_specimen_side") if selected is not None else None,
            "candidate_replicate_number": selected.get("replicate_number") if selected is not None else None,
            "candidate_MHS_value": float(selected[mhs_col]) if selected is not None else np.nan,
            "candidate_Sysmex_value": float(selected[ref_col]) if selected is not None else np.nan,
            "candidate_residual_pct": float(selected["residual_pct"]) if selected is not None else np.nan,
            "removal_scope": "one linked MHS/Sysmex raw replicate for this analyte only" if decision["removed"] else "none",
        }
        audit_rows.append(audit)
        if decision["removed"]:
            removed_rows.append(audit)

    return {
        "base_rows": base,
        "global_exclusions": global_exclusions,
        "analyte_flag_exclusions": pd.DataFrame(analyte_flag_exclusions),
        "filtered_by_analyte": filtered_by_analyte,
        "outlier_audit": pd.DataFrame(audit_rows),
        "removed_outliers": pd.DataFrame(removed_rows),
    }


def build_ba_donor_rows(
    filtered_by_analyte: Mapping[str, pd.DataFrame],
    specimen_a_label: str,
    specimen_b_label: str,
    analytes: Mapping[str, Mapping[str, object]] = ANALYTES,
) -> tuple[dict[str, list[dict[str, object]]], pd.DataFrame]:
    donor_rows: dict[str, list[dict[str, object]]] = {}
    long_rows: list[dict[str, object]] = []
    for analyte, cfg in analytes.items():
        mhs_col, ref_col = str(cfg["mhs"]), str(cfg["ref"])
        sub = filtered_by_analyte[analyte].copy()
        rows: list[dict[str, object]] = []
        for donor in sorted(sub["donor"].dropna().astype(str).unique(), key=natural_donor_key):
            a = sub[(sub["donor"].astype(str) == donor) & (sub["BA_specimen_side"] == specimen_a_label)]
            b = sub[(sub["donor"].astype(str) == donor) & (sub["BA_specimen_side"] == specimen_b_label)]
            if a.empty or b.empty:
                continue
            rec = {
                "analyte": analyte, "unit": cfg["unit"], "donor": donor,
                "a_label": specimen_a_label, "b_label": specimen_b_label,
                "a_rows": a[[mhs_col, ref_col]].to_numpy(float),
                "b_rows": b[[mhs_col, ref_col]].to_numpy(float),
                "mhs_a": float(a[mhs_col].mean()), "sys_a": float(a[ref_col].mean()),
                "mhs_b": float(b[mhs_col].mean()), "sys_b": float(b[ref_col].mean()),
                "a_replicate_n": int(len(a)), "b_replicate_n": int(len(b)),
            }
            rows.append(rec)
            long_rows.append({k: v for k, v in rec.items() if k not in {"a_rows", "b_rows"}})
        donor_rows[analyte] = rows
    return donor_rows, pd.DataFrame(long_rows)


def ba_method_values(record: Mapping[str, object]) -> dict[str, float]:
    a, b, ra, rb = float(record["mhs_a"]), float(record["mhs_b"]), float(record["sys_a"]), float(record["sys_b"])
    return {
        "M02_pct": 100.0 * (a - b) / ((a + b) / 2.0),
        "M11_pct": 100.0 * (a - ra) / ra - 100.0 * (b - rb) / rb,
        "M05_pct": 100.0 * (abs(rb - ra) - abs(b - a)) / ((rb + ra) / 2.0),
        "S02_pct": 100.0 * (ra - rb) / ((ra + rb) / 2.0),
        "M02_native": a - b,
        "M11_native": (a - ra) - (b - rb),
        "M05_native": abs(rb - ra) - abs(b - a),
        "S02_native": ra - rb,
        "M02_x": (a + b) / 2.0,
        "M11_x": (ra + rb) / 2.0,
        "M05_x": (ra + rb) / 2.0,
    }


def ba_summary(values: Sequence[float], rng: np.random.Generator, bootstrap_iterations: int = BA_BOOTSTRAP_ITERATIONS) -> dict[str, float]:
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    n = len(v)
    if n < 3:
        return {k: np.nan for k in [
            "mean_bias", "sd", "mean_ci_low", "mean_ci_high", "loa_low", "loa_high",
            "loa_low_ci_low", "loa_low_ci_high", "loa_high_ci_low", "loa_high_ci_high",
            "bootstrap_mean_ci_low", "bootstrap_mean_ci_high", "bootstrap_mag", "bootstrap_mag_ci_low",
            "bootstrap_mag_ci_high", "bootstrap_loa_low_ci_low", "bootstrap_loa_low_ci_high",
            "bootstrap_loa_high_ci_low", "bootstrap_loa_high_ci_high", "empirical_loa_low",
            "empirical_loa_high", "shapiro_W", "shapiro_p"
        ]} | {"n": n}
    mean_bias = float(np.mean(v)); sd = float(np.std(v, ddof=1))
    tcrit = float(stats.t.ppf(0.975, n - 1)); mean_hw = tcrit * sd / np.sqrt(n)
    loa_low = mean_bias - 1.96 * sd; loa_high = mean_bias + 1.96 * sd
    loa_se = sd * np.sqrt(1.0/n + 1.96**2/(2.0*(n-1))); loa_hw = tcrit * loa_se
    sw_w = sw_p = np.nan
    if 3 <= n <= 5000 and sd > 0:
        sw = shapiro(v); sw_w = float(sw.statistic); sw_p = float(sw.pvalue)
    idx = rng.integers(0, n, size=(bootstrap_iterations, n))
    boot = v[idx]
    boot_means = np.mean(boot, axis=1)
    boot_sds = np.std(boot, axis=1, ddof=1)
    boot_loa_low = boot_means - 1.96 * boot_sds
    boot_loa_high = boot_means + 1.96 * boot_sds
    empirical_low, empirical_high = np.percentile(v, [2.5, 97.5])
    return {
        "n": n, "mean_bias": mean_bias, "sd": sd,
        "mean_ci_low": float(mean_bias - mean_hw), "mean_ci_high": float(mean_bias + mean_hw),
        "loa_low": float(loa_low), "loa_high": float(loa_high),
        "loa_low_ci_low": float(loa_low - loa_hw), "loa_low_ci_high": float(loa_low + loa_hw),
        "loa_high_ci_low": float(loa_high - loa_hw), "loa_high_ci_high": float(loa_high + loa_hw),
        "bootstrap_mean_ci_low": float(np.percentile(boot_means, 2.5)),
        "bootstrap_mean_ci_high": float(np.percentile(boot_means, 97.5)),
        "bootstrap_mag": float(abs(mean_bias)),
        "bootstrap_mag_ci_low": float(np.percentile(np.abs(boot_means), 2.5)),
        "bootstrap_mag_ci_high": float(np.percentile(np.abs(boot_means), 97.5)),
        "bootstrap_loa_low_ci_low": float(np.percentile(boot_loa_low, 2.5)),
        "bootstrap_loa_low_ci_high": float(np.percentile(boot_loa_low, 97.5)),
        "bootstrap_loa_high_ci_low": float(np.percentile(boot_loa_high, 2.5)),
        "bootstrap_loa_high_ci_high": float(np.percentile(boot_loa_high, 97.5)),
        "empirical_loa_low": float(empirical_low), "empirical_loa_high": float(empirical_high),
        "shapiro_W": sw_w, "shapiro_p": sw_p,
    }


def status_within(summary: Mapping[str, float], limit: float | None) -> dict[str, object]:
    if limit is None or int(summary.get("n", 0)) < 3:
        return {"mean": None, "point": None, "full": None}
    return {
        "mean": bool(summary["mean_ci_low"] >= -limit and summary["mean_ci_high"] <= limit),
        "point": bool(summary["loa_low"] >= -limit and summary["loa_high"] <= limit),
        "full": bool(summary["loa_low_ci_low"] >= -limit and summary["loa_high_ci_high"] <= limit),
    }


def _pf(value: object) -> str:
    if value is None or pd.isna(value): return "N/A"
    return "Pass" if bool(value) else "Fail"


def compute_ba_results(
    donor_rows: Mapping[str, Sequence[Mapping[str, object]]],
    outlier_audit: pd.DataFrame,
    criteria_confirmed: bool,
    analytes: Mapping[str, Mapping[str, object]] = ANALYTES,
) -> dict[str, pd.DataFrame]:
    rng = np.random.default_rng(BA_RANDOM_SEED)
    percent_rows: list[dict[str, object]] = []
    native_rows: list[dict[str, object]] = []
    combined_rows: list[dict[str, object]] = []
    donor_values: list[dict[str, object]] = []
    audit_lookup = outlier_audit.set_index("analyte").to_dict("index") if not outlier_audit.empty else {}
    for analyte, cfg in analytes.items():
        records = list(donor_rows.get(analyte, []))
        calculations = []
        for record in records:
            vals = ba_method_values(record)
            calculations.append(vals)
            donor_values.append({
                "analyte": analyte, "unit": cfg["unit"], "donor": record["donor"],
                "a_replicate_n": record["a_replicate_n"], "b_replicate_n": record["b_replicate_n"],
                "mhs_a": record["mhs_a"], "sys_a": record["sys_a"], "mhs_b": record["mhs_b"], "sys_b": record["sys_b"],
                **vals,
            })
        normal_low, normal_high = map(float, cfg["normal"])
        midpoint = (normal_low + normal_high) / 2.0
        audit = audit_lookup.get(analyte, {})
        for method in ["M02", "M11", "M05"]:
            pct = ba_summary([c[f"{method}_pct"] for c in calculations], rng)
            native = ba_summary([c[f"{method}_native"] for c in calculations], rng)
            ac_limit = AC_LIMIT_PCT.get(analyte); clia_limit = CLIA_LIMIT_PCT.get(analyte)
            ac = status_within(pct, ac_limit); clia = status_within(pct, clia_limit)
            pct_row = {
                "analyte": analyte, "unit": cfg["unit"], "method": method,
                "method_description": METHOD_DESCRIPTIONS[method], **pct,
                "AC_limit_pct": ac_limit, "AC_mean_CI": _pf(ac["mean"]), "AC_point_LoA": _pf(ac["point"]), "AC_full_LoA": _pf(ac["full"]),
                "CLIA_limit_pct": clia_limit, "CLIA_mean_CI": _pf(clia["mean"]), "CLIA_point_LoA": _pf(clia["point"]), "CLIA_full_LoA": _pf(clia["full"]),
                "criteria_status": "USER VERIFIED" if criteria_confirmed else "UNVERIFIED / PROVISIONAL CONTEXT ONLY",
                "outlier_removed": bool(audit.get("removed", False)), "outlier_method": audit.get("method"),
                "outlier_donor": audit.get("candidate_donor") if audit.get("removed") else None,
                "outlier_specimen": audit.get("candidate_specimen_type") if audit.get("removed") else None,
            }
            native_row = {"analyte": analyte, "unit": cfg["unit"], "method": method,
                          "method_description": METHOD_DESCRIPTIONS[method], **native}
            percent_rows.append(pct_row); native_rows.append(native_row)
            combined_rows.append({
                "analyte": analyte, "unit": cfg["unit"], "method": method,
                "method_description": METHOD_DESCRIPTIONS[method], "n": pct["n"],
                "normal_low": normal_low, "normal_high": normal_high, "normal_midpoint": midpoint,
                "50pct_of_normal_midpoint": 0.5 * midpoint,
                "mean_bias_pct": pct["mean_bias"], "mean_bias_pct_CI_low": pct["mean_ci_low"], "mean_bias_pct_CI_high": pct["mean_ci_high"],
                "LoA_pct_low": pct["loa_low"], "LoA_pct_high": pct["loa_high"],
                "LoA_pct_low_endpoint_CI_low": pct["loa_low_ci_low"], "LoA_pct_low_endpoint_CI_high": pct["loa_low_ci_high"],
                "LoA_pct_high_endpoint_CI_low": pct["loa_high_ci_low"], "LoA_pct_high_endpoint_CI_high": pct["loa_high_ci_high"],
                "exact_native_mean": native["mean_bias"], "exact_native_mean_CI_low": native["mean_ci_low"], "exact_native_mean_CI_high": native["mean_ci_high"],
                "exact_native_LoA_low": native["loa_low"], "exact_native_LoA_high": native["loa_high"],
                "approx_midpoint_scaled_mean_units": pct["mean_bias"] * midpoint / 100.0 if np.isfinite(pct["mean_bias"]) else np.nan,
                "approx_midpoint_scaled_LoA_low_units": pct["loa_low"] * midpoint / 100.0 if np.isfinite(pct["loa_low"]) else np.nan,
                "approx_midpoint_scaled_LoA_high_units": pct["loa_high"] * midpoint / 100.0 if np.isfinite(pct["loa_high"]) else np.nan,
                "AC_limit_pct": ac_limit, "AC_approx_units_at_midpoint": ac_limit * midpoint / 100.0 if ac_limit is not None else np.nan,
                "AC_mean_CI": _pf(ac["mean"]), "AC_point_LoA": _pf(ac["point"]), "AC_full_LoA": _pf(ac["full"]),
                "CLIA_limit_pct": clia_limit, "CLIA_approx_units_at_midpoint": clia_limit * midpoint / 100.0 if clia_limit is not None else np.nan,
                "CLIA_mean_CI": _pf(clia["mean"]), "CLIA_point_LoA": _pf(clia["point"]), "CLIA_full_LoA": _pf(clia["full"]),
                "criteria_status": pct_row["criteria_status"],
            })
    return {
        "percent_results": pd.DataFrame(percent_rows),
        "native_results": pd.DataFrame(native_rows),
        "combined_context": pd.DataFrame(combined_rows),
        "donor_values": pd.DataFrame(donor_values),
    }


def bootstrap_within_donor_profile(record: Mapping[str, object], rng: np.random.Generator,
                                   iterations: int = DONOR_PROFILE_BOOTSTRAPS) -> dict[str, float]:
    a_rows = np.asarray(record["a_rows"], dtype=float); b_rows = np.asarray(record["b_rows"], dtype=float)
    point = ba_method_values(record)
    ia = rng.integers(0, len(a_rows), size=(iterations, len(a_rows)))
    ib = rng.integers(0, len(b_rows), size=(iterations, len(b_rows)))
    a = a_rows[ia].mean(axis=1); b = b_rows[ib].mean(axis=1)
    mhs_a, sys_a = a[:, 0], a[:, 1]; mhs_b, sys_b = b[:, 0], b[:, 1]
    boot = {
        "M11_pct": 100*(mhs_a-sys_a)/sys_a - 100*(mhs_b-sys_b)/sys_b,
        "S02_pct": 100*(sys_a-sys_b)/((sys_a+sys_b)/2),
        "M11_native": (mhs_a-sys_a)-(mhs_b-sys_b),
        "S02_native": sys_a-sys_b,
    }
    out = {}
    for key, values in boot.items():
        values = values[np.isfinite(values)]
        lo, hi = np.percentile(values, [2.5, 97.5])
        out[key] = point[key]; out[f"{key}_ci_low"] = float(lo); out[f"{key}_ci_high"] = float(hi)
    return out


def _asym_yerr(point, low, high):
    p = np.asarray(point, float); lo = np.asarray(low, float); hi = np.asarray(high, float)
    return np.vstack([np.maximum(p-lo, 0), np.maximum(hi-p, 0)])


def make_donor_profile_records(donor_rows: Mapping[str, Sequence[Mapping[str, object]]]) -> pd.DataFrame:
    rng = np.random.default_rng(BA_RANDOM_SEED + 111)
    records = []
    for analyte, rows in donor_rows.items():
        if not rows: continue
        ref_scale = float(np.mean([(r["sys_a"] + r["sys_b"])/2 for r in rows]))
        native_ac = AC_LIMIT_PCT[analyte] / 100.0 * ref_scale
        for record in rows:
            vals = bootstrap_within_donor_profile(record, rng)
            records.append({"analyte": analyte, "donor": record["donor"], **vals,
                            "AC_pct": AC_LIMIT_PCT[analyte], "AC_native_approx": native_ac})
    return pd.DataFrame(records)


def save_donor_profile_plots(profile_df: pd.DataFrame, output_dir: Path,
                             analytes: Mapping[str, Mapping[str, object]] = ANALYTES,
                             plot_mode: str = "full") -> dict[str, Path]:
    output_dir = Path(output_dir); pct_dir = output_dir / "donor_profiles_percent"; nat_dir = output_dir / "donor_profiles_native"
    pct_dir.mkdir(parents=True, exist_ok=True); nat_dir.mkdir(parents=True, exist_ok=True)
    mhs_color="#1f77b4"; mhs_ci="#add8e6"; sys_color="#2ca02c"; sys_ci="#90ee90"; ac_color="#d62728"
    paths = {}
    def draw(ax, analyte, native=False):
        rows = profile_df[profile_df["analyte"] == analyte].sort_values("donor", key=lambda s: s.map(natural_donor_key))
        donors = rows["donor"].tolist(); x=np.arange(len(donors),dtype=float); off=.14
        if not native:
            m=rows["M11_pct"].to_numpy(float); ml=rows["M11_pct_ci_low"].to_numpy(float); mh=rows["M11_pct_ci_high"].to_numpy(float)
            s=rows["S02_pct"].to_numpy(float); sl=rows["S02_pct_ci_low"].to_numpy(float); sh=rows["S02_pct_ci_high"].to_numpy(float)
            lim=AC_LIMIT_PCT[analyte]; ylabel="Mean Bland-Altman bias (%)"
            labels=("Reference-adjusted MHS capillary vs venous", "Sysmex capillary vs venous", "Contextual AC trueness limits")
        else:
            m=rows["M11_native"].to_numpy(float); ml=rows["M11_native_ci_low"].to_numpy(float); mh=rows["M11_native_ci_high"].to_numpy(float)
            s=rows["S02_native"].to_numpy(float); sl=rows["S02_native_ci_low"].to_numpy(float); sh=rows["S02_native_ci_high"].to_numpy(float)
            lim=float(rows["AC_native_approx"].iloc[0]); ylabel=f"Native mean difference ({analytes[analyte]['unit']})"
            labels=("MHS matched raw-error difference", "Sysmex capillary - venous", "Approximate native AC limits")
        ax.errorbar(x-off,m,yerr=_asym_yerr(m,ml,mh),fmt="o",markersize=4,capsize=2,color=mhs_color,ecolor=mhs_ci,elinewidth=1,label=labels[0])
        ax.errorbar(x+off,s,yerr=_asym_yerr(s,sl,sh),fmt="o",markersize=4,capsize=2,color=sys_color,ecolor=sys_ci,elinewidth=1,label=labels[1])
        ax.axhline(lim,color=ac_color,linestyle=":",linewidth=1.2,label=labels[2]); ax.axhline(-lim,color=ac_color,linestyle=":",linewidth=1.2)
        ax.axhline(0,color="black",linewidth=.7,alpha=.55); ax.set_title(analyte); ax.set_ylabel(ylabel)
        ax.set_xticks(x); ax.set_xticklabels(donors,rotation=45,ha="right"); ax.set_xlabel("Donor"); ax.grid(axis="y",alpha=.18); ax.legend(loc="best",fontsize=8)
    if plot_mode == "full":
      for analyte in analytes:
        if profile_df[profile_df["analyte"] == analyte].empty: continue
        for native, folder, suffix in [(False,pct_dir,"percent"),(True,nat_dir,"native")]:
            fig,ax=plt.subplots(figsize=(10,5)); draw(ax,analyte,native); fig.tight_layout()
            path=folder/f"{analyte}_M11_vs_Sysmex_{suffix}_bias_FILTERED.png"; fig.savefig(path,dpi=160); plt.close(fig)
            paths[f"{analyte}_{suffix}"]=path
    for native, name in [(False,"Appendix_B1_Percent_Donor_Profiles_FILTERED.png"),(True,"Appendix_B2_Native_Donor_Profiles_FILTERED.png")]:
        fig,axes=plt.subplots(4,3,figsize=(22,14.5)); letters="ABCDEFGHIJKL"
        for i,(ax,analyte) in enumerate(zip(axes.ravel(), analytes)):
            if profile_df[profile_df["analyte"] == analyte].empty:
                ax.axis("off"); continue
            draw(ax,analyte,native)
            ax.text(-0.04,1.04,letters[i],transform=ax.transAxes,fontsize=18,fontweight="bold",va="top",ha="left",bbox=dict(boxstyle="round,pad=0.08",facecolor="white",edgecolor="black",linewidth=1))
        fig.tight_layout(pad=1.1); path=output_dir/name; fig.savefig(path,dpi=150); plt.close(fig); paths[name]=path
    return paths


def save_standard_bland_altman_plots(donor_values: pd.DataFrame, combined: pd.DataFrame, output_dir: Path,
                                      plot_mode: str = "full") -> dict[str, Path]:
    output_dir=Path(output_dir); output_dir.mkdir(parents=True,exist_ok=True); paths={}
    for unit_kind in ["pct","native"]:
        for method in ["M02","M11","M05"]:
            subdir=output_dir/f"{method}_{unit_kind}"; subdir.mkdir(parents=True,exist_ok=True)
            fig_panel, axes = plt.subplots(4,3,figsize=(18,18))
            for ax, analyte in zip(axes.ravel(), ANALYTES):
                vals=donor_values[donor_values["analyte"]==analyte].copy()
                summ=combined[(combined["analyte"]==analyte)&(combined["method"]==method)]
                if vals.empty or summ.empty:
                    ax.axis("off"); continue
                row=summ.iloc[0]
                x=vals[f"{method}_x"].to_numpy(float)
                y=vals[f"{method}_{unit_kind}"].to_numpy(float)
                mean_bias=float(row["mean_bias_pct"] if unit_kind=="pct" else row["exact_native_mean"])
                loa_low=float(row["LoA_pct_low"] if unit_kind=="pct" else row["exact_native_LoA_low"])
                loa_high=float(row["LoA_pct_high"] if unit_kind=="pct" else row["exact_native_LoA_high"])
                ax.scatter(x,y,s=34,color="#1f77b4",alpha=.9)
                ax.axhline(mean_bias,color="black",linewidth=1.2,label="Mean bias")
                ax.axhline(loa_low,color="#d62728",linestyle="--",linewidth=1.2,label="95% LoA")
                ax.axhline(loa_high,color="#d62728",linestyle="--",linewidth=1.2)
                if unit_kind=="pct":
                    lim=AC_LIMIT_PCT.get(analyte)
                    if lim is not None:
                        ax.axhline(lim,color="#d62728",linestyle=":",linewidth=1.0,label="Contextual AC")
                        ax.axhline(-lim,color="#d62728",linestyle=":",linewidth=1.0)
                ax.set_title(f"{analyte} — {method}",fontweight="bold")
                ax.set_xlabel("Mean MHS level" if method=="M02" else "Mean Sysmex reference level")
                ax.set_ylabel("Difference (%)" if unit_kind=="pct" else f"Difference ({ANALYTES[analyte]['unit']})")
                ax.grid(True,alpha=.22)
                ax.text(.02,.98,f"n={len(y)}\nMean={mean_bias:.3g}\nLoA={loa_low:.3g} to {loa_high:.3g}",transform=ax.transAxes,va="top",fontsize=8)
                ax.legend(loc="best",fontsize=7)
                if plot_mode == "full":
                    fig, ax2 = plt.subplots(figsize=(7,5))
                    ax2.scatter(x,y,s=40,color="#1f77b4",alpha=.9)
                    ax2.axhline(mean_bias,color="black",linewidth=1.3,label="Mean bias")
                    ax2.axhline(loa_low,color="#d62728",linestyle="--",linewidth=1.3,label="95% LoA")
                    ax2.axhline(loa_high,color="#d62728",linestyle="--",linewidth=1.3)
                    if unit_kind=="pct" and AC_LIMIT_PCT.get(analyte) is not None:
                        lim=AC_LIMIT_PCT[analyte]; ax2.axhline(lim,color="#d62728",linestyle=":",linewidth=1,label="Contextual AC"); ax2.axhline(-lim,color="#d62728",linestyle=":",linewidth=1)
                    ax2.set_title(f"{analyte} — {METHOD_DESCRIPTIONS[method]}",fontweight="bold")
                    ax2.set_xlabel("Mean MHS level" if method=="M02" else "Mean Sysmex reference level")
                    ax2.set_ylabel("Difference (%)" if unit_kind=="pct" else f"Difference ({ANALYTES[analyte]['unit']})")
                    ax2.grid(True,alpha=.22); ax2.legend(loc="best",fontsize=8)
                    ax2.text(.02,.98,f"n={len(y)}\nMean={mean_bias:.3g}\n95% LoA={loa_low:.3g} to {loa_high:.3g}",transform=ax2.transAxes,va="top",fontsize=9)
                    fig.tight_layout(); path=subdir/f"{analyte}_{method}_Bland_Altman_{unit_kind}.png"; fig.savefig(path,dpi=140); plt.close(fig); paths[f"{analyte}_{method}_{unit_kind}"]=path
            fig_panel.suptitle(f"{method} Bland–Altman plots — {'percentage' if unit_kind=='pct' else 'native units'}",fontsize=15,fontweight="bold")
            fig_panel.tight_layout(rect=[0,0,1,.98]); panel_path=output_dir/f"{method}_Bland_Altman_panel_{unit_kind}.png"; fig_panel.savefig(panel_path,dpi=140); plt.close(fig_panel); paths[f"{method}_panel_{unit_kind}"]=panel_path
    return paths


def export_combined_workbook(path: Path, regression_result: Mapping[str, object], ba_bundle: Mapping[str, object], criteria_df: pd.DataFrame) -> None:
    path=Path(path)
    with pd.ExcelWriter(path,engine="xlsxwriter") as writer:
        wb=writer.book
        header=wb.add_format({"bold":True,"font_color":"white","bg_color":"#1F4E78","border":1,"align":"center","valign":"vcenter","text_wrap":True})
        pass_fmt=wb.add_format({"bg_color":"#D9EAD3","border":1}); fail_fmt=wb.add_format({"bg_color":"#F4CCCC","border":1}); warn_fmt=wb.add_format({"bg_color":"#FFF2CC","border":1})
        reg_metrics=pd.concat([b["metrics"] for b in regression_result["groups"].values()],ignore_index=True) if regression_result.get("groups") else pd.DataFrame()
        sheets={
            "Regression_Summary":reg_metrics,
            "BA_Combined_Context":ba_bundle["combined_context"],
            "BA_Percent":ba_bundle["percent_results"],
            "BA_Native":ba_bundle["native_results"],
            "BA_Donor_Values":ba_bundle["donor_values"],
            "BA_Donor_Means":ba_bundle["donor_means"],
            "BA_Outlier_Audit":ba_bundle["outlier_audit"],
            "BA_Removed_Outliers":ba_bundle["removed_outliers"],
            "BA_AnalyteFlag_Excluded":ba_bundle["analyte_flag_exclusions"],
            "Global_Exclusions":ba_bundle["global_exclusions"],
            "Acceptance_Criteria":criteria_df,
        }
        for name,df in sheets.items():
            df.to_excel(writer,sheet_name=name,index=False)
            ws=writer.sheets[name]; ws.freeze_panes(1,0); ws.set_row(0,32,header); ws.autofilter(0,0,max(1,len(df)),max(0,len(df.columns)-1)); ws.set_column(0,max(0,len(df.columns)-1),18)
            if name in {"BA_Combined_Context","BA_Percent"}:
                for cname in [c for c in df.columns if c.endswith("CI") or c.endswith("LoA")]:
                    col=df.columns.get_loc(cname)
                    ws.conditional_format(1,col,max(1,len(df)),col,{"type":"text","criteria":"containing","value":"Pass","format":pass_fmt})
                    ws.conditional_format(1,col,max(1,len(df)),col,{"type":"text","criteria":"containing","value":"Fail","format":fail_fmt})
                if "criteria_status" in df.columns:
                    col=df.columns.get_loc("criteria_status"); ws.conditional_format(1,col,max(1,len(df)),col,{"type":"text","criteria":"containing","value":"UNVERIFIED","format":warn_fmt})
        for group_name,bundle in regression_result.get("groups",{}).items():
            short=re.sub(r"[^A-Za-z0-9]","",group_name)[:20]
            for prefix,key in [("RegDonor_","donor_means"),("RegOutlier_","outliers"),("RegESD_","diagnostics"),("RegRepN_","replicate_counts")]:
                df=bundle[key]; sname=(prefix+short)[:31]; df.to_excel(writer,sheet_name=sname,index=False); ws=writer.sheets[sname]; ws.freeze_panes(1,0); ws.set_row(0,28,header); ws.set_column(0,max(0,len(df.columns)-1),17)


def run_ba_pipeline(
    input_csv: str | Path,
    output_root: str | Path,
    specimen_a_types: Sequence[str] = ("Mc", "Mc2"),
    specimen_b_types: Sequence[str] = ("MK2",),
    specimen_a_label: str = "Capillary",
    specimen_b_label: str = "MK2",
    criteria_confirmed: bool = False,
    plot_mode: str = "full",
) -> dict[str, object]:
    output_root=Path(output_root); output_root.mkdir(parents=True,exist_ok=True)
    df=load_and_prepare(input_csv)
    filt=build_ba_filtered_replicates(df,specimen_a_types,specimen_b_types,specimen_a_label,specimen_b_label)
    donor_rows, donor_means=build_ba_donor_rows(filt["filtered_by_analyte"],specimen_a_label,specimen_b_label)
    results=compute_ba_results(donor_rows,filt["outlier_audit"],criteria_confirmed)
    profile_df=make_donor_profile_records(donor_rows)
    profile_paths=save_donor_profile_plots(profile_df,output_root/"plots_donor_profiles",plot_mode=plot_mode)
    ba_paths=save_standard_bland_altman_plots(results["donor_values"],results["combined_context"],output_root/"plots_bland_altman",plot_mode=plot_mode)
    criteria_df=pd.DataFrame([{
        "analyte":a,"unit":ANALYTES[a]["unit"],"normal_min":ANALYTES[a]["normal"][0],"normal_max":ANALYTES[a]["normal"][1],
        "AC_limit_pct":AC_LIMIT_PCT.get(a),"CLIA_limit_pct":CLIA_LIMIT_PCT.get(a),
        "criteria_confirmed_by_user":criteria_confirmed,
        "note":"Review/edit before final use; contextual trueness benchmark, not established capillary-versus-venous acceptance criterion."
    } for a in ANALYTES])
    bundle={**filt,**results,"donor_rows":donor_rows,"donor_means":donor_means,"profile_records":profile_df,"profile_paths":profile_paths,"ba_plot_paths":ba_paths,"criteria":criteria_df}
    # CSV outputs
    frames={
        "BA_PERCENT_RESULTS.csv":bundle["percent_results"],"BA_NATIVE_RESULTS.csv":bundle["native_results"],
        "BA_PERCENT_NATIVE_CONTEXT.csv":bundle["combined_context"],"BA_DONOR_VALUES.csv":bundle["donor_values"],
        "BA_DONOR_MEANS.csv":bundle["donor_means"],"BA_OUTLIER_AUDIT.csv":bundle["outlier_audit"],
        "BA_REMOVED_OUTLIERS.csv":bundle["removed_outliers"],"BA_ANALYTE_FLAG_EXCLUSIONS.csv":bundle["analyte_flag_exclusions"],
        "GLOBAL_FLAG_EXCLUSIONS.csv":bundle["global_exclusions"],"BA_DONOR_PROFILE_DATA.csv":bundle["profile_records"],
        "AC_CLIA_ACCEPTANCE_CRITERIA_REVIEW.csv":criteria_df,
    }
    for filename,frame in frames.items(): frame.to_csv(output_root/filename,index=False)
    with (output_root/"BA_METHOD_NOTES.txt").open("w",encoding="utf-8") as h:
        h.write("M02 is direct signed MHS capillary-versus-venous.\nM11 is signed Sysmex reference-adjusted and is the preferred reference-adjusted option.\nM05 is absolute-gap magnitude-only sensitivity analysis and is not signed.\nOutlier screening is on raw linked replicate residuals before donor averaging: Shapiro-Wilk; Grubbs with manually computed Gcrit if normal; MAD modified Z if non-normal; maximum one replicate per analyte.\n")
    return bundle


# =============================================================================
# Streamlit-app integration helpers
# =============================================================================
import contextlib
import copy
import tempfile
import threading
import zipfile

DEFAULT_ANALYTE_CONFIG = copy.deepcopy(ANALYTES)
DEFAULT_AC_LIMIT_PCT = copy.deepcopy(AC_LIMIT_PCT)
DEFAULT_CLIA_LIMIT_PCT = copy.deepcopy(CLIA_LIMIT_PCT)
_RUNTIME_LOCK = threading.RLock()


def default_analyte_table() -> pd.DataFrame:
    """Return the notebook defaults as an editable app configuration table."""
    rows = []
    for analyte, cfg in DEFAULT_ANALYTE_CONFIG.items():
        rows.append({
            "Analyte": analyte,
            "MHS column": str(cfg["mhs"]),
            "Reference column": str(cfg["ref"]),
            "Flag column": f"{analyte}_flag",
            "Normal low": float(cfg["normal"][0]),
            "Normal high": float(cfg["normal"][1]),
            "Unit": str(cfg["unit"]),
            "AC limit %": float(DEFAULT_AC_LIMIT_PCT.get(analyte, np.nan)),
            "CLIA limit %": float(DEFAULT_CLIA_LIMIT_PCT.get(analyte, np.nan)),
        })
    return pd.DataFrame(rows)


def _finite_or_none(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def build_runtime_analytes(config_rows: Sequence[Mapping[str, object]]) -> tuple[dict[str, dict[str, object]], dict[str, float], dict[str, float]]:
    """Convert the app's edited analyte table into notebook runtime dictionaries."""
    runtime: dict[str, dict[str, object]] = {}
    ac: dict[str, float] = {}
    clia: dict[str, float] = {}
    for row in config_rows:
        analyte = str(row["Analyte"]).strip().upper()
        if analyte not in DEFAULT_ANALYTE_CONFIG:
            raise ValueError(f"Unsupported analyte: {analyte}")
        low = float(row["Normal low"])
        high = float(row["Normal high"])
        if not np.isfinite(low) or not np.isfinite(high) or high <= low:
            raise ValueError(f"{analyte}: normal high must be greater than normal low.")
        runtime[analyte] = {
            # Internal canonical output names remain identical to the notebook.
            "mhs": str(DEFAULT_ANALYTE_CONFIG[analyte]["mhs"]),
            "ref": str(DEFAULT_ANALYTE_CONFIG[analyte]["ref"]),
            "normal": (low, high),
            "unit": str(row["Unit"]).strip() or str(DEFAULT_ANALYTE_CONFIG[analyte]["unit"]),
            "source_mhs": str(row["MHS column"]),
            "source_ref": str(row["Reference column"]),
            "source_flag": None if str(row.get("Flag column", "None")) in {"", "None", "<none>"} else str(row.get("Flag column")),
        }
        ac_value = _finite_or_none(row.get("AC limit %"))
        clia_value = _finite_or_none(row.get("CLIA limit %"))
        if ac_value is None or ac_value <= 0:
            raise ValueError(f"{analyte}: AC limit must be a positive percentage.")
        ac[analyte] = ac_value
        if clia_value is not None and clia_value > 0:
            clia[analyte] = clia_value
    if not runtime:
        raise ValueError("Select at least one analyte.")
    return runtime, ac, clia


def canonicalize_input_dataframe(
    data: pd.DataFrame,
    runtime_analytes: Mapping[str, Mapping[str, object]],
    *,
    batch_id_col: str | None,
    global_flag_col: str | None,
    sample_id_col: str | None = None,
    donor_col: str | None = None,
    specimen_col: str | None = None,
    replicate_col: str | None = None,
    treat_missing_global_flag_as_false: bool = False,
) -> pd.DataFrame:
    """Create the canonical notebook input while preserving the uploaded columns."""
    if data.empty:
        raise ValueError("The uploaded dataset is empty.")
    out = data.copy().reset_index(drop=True)

    if batch_id_col is None:
        out["batch_id"] = [f"row_{i+1:06d}" for i in range(len(out))]
    else:
        if batch_id_col not in out.columns:
            raise ValueError(f"Batch ID column not found: {batch_id_col}")
        out["batch_id"] = out[batch_id_col].astype(str)

    if global_flag_col is None:
        if not treat_missing_global_flag_as_false:
            raise ValueError("Select a global flag column or explicitly allow all rows to be treated as global_flag = FALSE.")
        out["global_flag"] = False
    else:
        if global_flag_col not in out.columns:
            raise ValueError(f"Global flag column not found: {global_flag_col}")
        out["global_flag"] = out[global_flag_col]

    if sample_id_col is not None:
        if sample_id_col not in out.columns:
            raise ValueError(f"Sample ID column not found: {sample_id_col}")
        out["bloodSampleId"] = out[sample_id_col].astype(str)
    else:
        required = {"donor": donor_col, "specimen": specimen_col, "replicate": replicate_col}
        missing = [name for name, col in required.items() if not col or col not in out.columns]
        if missing:
            raise ValueError(
                "To synthesize bloodSampleId, select donor, specimen-type and replicate columns. "
                f"Missing: {', '.join(missing)}"
            )
        donor_text = out[str(donor_col)].astype(str).str.strip()
        specimen_text = out[str(specimen_col)].astype(str).str.strip()
        replicate_numeric = pd.to_numeric(out[str(replicate_col)], errors="coerce")
        fallback = out.groupby([donor_text, specimen_text], dropna=False).cumcount() + 1
        replicate_text = replicate_numeric.where(replicate_numeric.notna(), fallback).astype(int).astype(str)
        out["bloodSampleId"] = donor_text + "-" + specimen_text + "-" + replicate_text

    # Copy mapped analyte columns to the canonical names used by the notebook.
    for analyte, cfg in runtime_analytes.items():
        source_mhs = str(cfg["source_mhs"])
        source_ref = str(cfg["source_ref"])
        if source_mhs not in out.columns:
            raise ValueError(f"{analyte}: MHS column not found: {source_mhs}")
        if source_ref not in out.columns:
            raise ValueError(f"{analyte}: reference column not found: {source_ref}")
        out[str(cfg["mhs"])] = out[source_mhs]
        out[str(cfg["ref"])] = out[source_ref]
        canonical_flag = f"{analyte}_flag"
        source_flag = cfg.get("source_flag")
        if source_flag is None:
            if canonical_flag in out.columns:
                out = out.drop(columns=[canonical_flag])
        else:
            if str(source_flag) not in out.columns:
                raise ValueError(f"{analyte}: flag column not found: {source_flag}")
            out[canonical_flag] = out[str(source_flag)]

    return out


@contextlib.contextmanager
def notebook_runtime_configuration(
    runtime_analytes: Mapping[str, Mapping[str, object]],
    ac_limits: Mapping[str, float],
    clia_limits: Mapping[str, float],
    *,
    shapiro_alpha: float = 0.05,
    outlier_alpha: float = 0.05,
    mad_z_threshold: float = 3.5,
    ba_bootstrap_iterations: int = 10_000,
    donor_profile_bootstraps: int = 5_000,
    random_seed: int = 20260723,
):
    """Safely apply per-session configuration to notebook globals.

    The original notebook functions use module-level dictionaries. Streamlit can
    serve multiple sessions in one Python process, so analyses are serialized and
    all globals are restored after each run.
    """
    global SHAPIRO_ALPHA, OUTLIER_ALPHA, OUTLIER_MAD_Z
    global BA_BOOTSTRAP_ITERATIONS, DONOR_PROFILE_BOOTSTRAPS, BA_RANDOM_SEED

    with _RUNTIME_LOCK:
        saved_analytes = copy.deepcopy(ANALYTES)
        saved_ac = copy.deepcopy(AC_LIMIT_PCT)
        saved_clia = copy.deepcopy(CLIA_LIMIT_PCT)
        saved_scalars = (
            SHAPIRO_ALPHA, OUTLIER_ALPHA, OUTLIER_MAD_Z,
            BA_BOOTSTRAP_ITERATIONS, DONOR_PROFILE_BOOTSTRAPS, BA_RANDOM_SEED,
        )
        try:
            ANALYTES.clear()
            ANALYTES.update({
                name: {
                    "mhs": cfg["mhs"], "ref": cfg["ref"],
                    "normal": tuple(cfg["normal"]), "unit": cfg["unit"],
                }
                for name, cfg in runtime_analytes.items()
            })
            AC_LIMIT_PCT.clear()
            AC_LIMIT_PCT.update({str(k): float(v) for k, v in ac_limits.items()})
            CLIA_LIMIT_PCT.clear()
            CLIA_LIMIT_PCT.update({str(k): float(v) for k, v in clia_limits.items()})
            SHAPIRO_ALPHA = float(shapiro_alpha)
            OUTLIER_ALPHA = float(outlier_alpha)
            OUTLIER_MAD_Z = float(mad_z_threshold)
            BA_BOOTSTRAP_ITERATIONS = int(ba_bootstrap_iterations)
            DONOR_PROFILE_BOOTSTRAPS = int(donor_profile_bootstraps)
            BA_RANDOM_SEED = int(random_seed)
            yield
        finally:
            ANALYTES.clear()
            ANALYTES.update(saved_analytes)
            AC_LIMIT_PCT.clear()
            AC_LIMIT_PCT.update(saved_ac)
            CLIA_LIMIT_PCT.clear()
            CLIA_LIMIT_PCT.update(saved_clia)
            (
                SHAPIRO_ALPHA, OUTLIER_ALPHA, OUTLIER_MAD_Z,
                BA_BOOTSTRAP_ITERATIONS, DONOR_PROFILE_BOOTSTRAPS, BA_RANDOM_SEED,
            ) = saved_scalars


def _safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return cleaned or "analysis"


def run_complete_notebook_pipeline(
    canonical_data: pd.DataFrame,
    run_root: str | Path,
    *,
    runtime_analytes: Mapping[str, Mapping[str, object]],
    ac_limits: Mapping[str, float],
    clia_limits: Mapping[str, float],
    regression_groups: Mapping[str, Sequence[str]],
    groups_to_run: Sequence[str],
    ba_specimen_a_types: Sequence[str],
    ba_specimen_b_types: Sequence[str],
    ba_specimen_a_label: str,
    ba_specimen_b_label: str,
    criteria_confirmed: bool,
    run_huber_bootstrap: bool = False,
    n_huber_bootstraps: int = 1000,
    plot_mode: str = "full",
    shapiro_alpha: float = 0.05,
    outlier_alpha: float = 0.05,
    mad_z_threshold: float = 3.5,
    ba_bootstrap_iterations: int = 10_000,
    donor_profile_bootstraps: int = 5_000,
    random_seed: int = 20260723,
    cross_specimen_groups: tuple[str, str] | None = None,
    regression_equation_position: str = "auto",
    regression_metrics_position: str = "auto",
    **compatibility_kwargs: object,
) -> dict[str, object]:
    """Run the notebook from beginning to end and save the same output families.

    ``compatibility_kwargs`` protects deployed Streamlit apps from mixed-file
    upgrades. Known historical aliases are accepted; unknown arguments still
    raise a clear error rather than being silently ignored.
    """
    alias_map = {
        "equation_position": "regression_equation_position",
        "metrics_position": "regression_metrics_position",
    }
    for old_name, new_name in alias_map.items():
        if old_name in compatibility_kwargs:
            value = str(compatibility_kwargs.pop(old_name))
            if new_name == "regression_equation_position":
                regression_equation_position = value
            else:
                regression_metrics_position = value
    if compatibility_kwargs:
        unexpected = ", ".join(sorted(compatibility_kwargs))
        raise TypeError(
            "run_complete_notebook_pipeline() received unsupported keyword "
            f"argument(s): {unexpected}"
        )
    run_root = Path(run_root)
    run_root.mkdir(parents=True, exist_ok=True)
    output_dir = run_root / "PROXIMA_COMBINED_RESULTS"
    trueness_dir = output_dir / "01_trueness_regression"
    ba_dir = output_dir / "02_bland_altman"
    output_dir.mkdir(parents=True, exist_ok=True)

    input_csv = run_root / "canonical_uploaded_input.csv"
    canonical_data.to_csv(input_csv, index=False)

    with notebook_runtime_configuration(
        runtime_analytes, ac_limits, clia_limits,
        shapiro_alpha=shapiro_alpha,
        outlier_alpha=outlier_alpha,
        mad_z_threshold=mad_z_threshold,
        ba_bootstrap_iterations=ba_bootstrap_iterations,
        donor_profile_bootstraps=donor_profile_bootstraps,
        random_seed=random_seed,
    ):
        regression_result = run_pipeline(
            input_csv=input_csv,
            output_root=trueness_dir,
            analysis_groups=regression_groups,
            groups_to_run=groups_to_run,
            ac_json=None,
            run_huber_bootstrap=run_huber_bootstrap,
            n_boot=n_huber_bootstraps,
            regression_equation_position=regression_equation_position,
            regression_metrics_position=regression_metrics_position,
        )

        if cross_specimen_groups is not None:
            group_a, group_b = cross_specimen_groups
            if group_a in regression_result["groups"] and group_b in regression_result["groups"]:
                cross_dir = trueness_dir / f"CrossSpecimen_{_safe_name(group_a)}_vs_{_safe_name(group_b)}"
                cross_dir.mkdir(parents=True, exist_ok=True)
                matched = build_cross_specimen_table(
                    regression_result["groups"][group_a],
                    regression_result["groups"][group_b],
                    name_a=_safe_name(group_a),
                    name_b=_safe_name(group_b),
                )
                cross_metrics = compute_cross_specimen_metrics(
                    matched,
                    name_a=_safe_name(group_a),
                    name_b=_safe_name(group_b),
                )
                matched.to_csv(cross_dir / "matched_donor_means.csv", index=False)
                cross_metrics.to_csv(cross_dir / "cross_specimen_regression_metrics.csv", index=False)
                regression_result["cross_specimen_matched"] = matched
                regression_result["cross_specimen_metrics"] = cross_metrics

        ba_result = run_ba_pipeline(
            input_csv=input_csv,
            output_root=ba_dir,
            specimen_a_types=list(ba_specimen_a_types),
            specimen_b_types=list(ba_specimen_b_types),
            specimen_a_label=ba_specimen_a_label,
            specimen_b_label=ba_specimen_b_label,
            criteria_confirmed=criteria_confirmed,
            plot_mode=plot_mode,
        )

        combined_xlsx = output_dir / "PROXIMA_Trueness_and_BlandAltman_Output.xlsx"
        export_combined_workbook(combined_xlsx, regression_result, ba_result, ba_result["criteria"])

        main_cols = [
            "analyte", "unit", "method", "method_description", "n",
            "normal_low", "normal_high", "50pct_of_normal_midpoint",
            "mean_bias_pct", "mean_bias_pct_CI_low", "mean_bias_pct_CI_high",
            "exact_native_mean", "exact_native_mean_CI_low", "exact_native_mean_CI_high",
            "LoA_pct_low", "LoA_pct_high", "exact_native_LoA_low", "exact_native_LoA_high",
            "AC_limit_pct", "AC_mean_CI", "AC_point_LoA", "AC_full_LoA",
            "CLIA_limit_pct", "CLIA_mean_CI", "CLIA_point_LoA", "CLIA_full_LoA",
            "criteria_status",
        ]
        method_summary = ba_result["combined_context"][main_cols].copy()
        summary_csv = output_dir / "PROXIMA_BA_MAIN_PERCENT_NATIVE_SUMMARY.csv"
        method_summary.to_csv(summary_csv, index=False)

        app_config = {
            "selected_analytes": list(runtime_analytes),
            "regression_groups": {k: list(v) for k, v in regression_groups.items()},
            "groups_to_run": list(groups_to_run),
            "BA_specimen_A_types": list(ba_specimen_a_types),
            "BA_specimen_B_types": list(ba_specimen_b_types),
            "BA_specimen_A_label": ba_specimen_a_label,
            "BA_specimen_B_label": ba_specimen_b_label,
            "acceptance_criteria_confirmed": bool(criteria_confirmed),
            "trueness_outlier_method": "Validated generalized ESD on raw replicate rows before donor averaging",
            "BA_outlier_method": "Shapiro-Wilk -> manual Grubbs Gcrit if normal; MAD modified-Z if non-normal; maximum one replicate per analyte before donor averaging",
            "run_huber_bootstrap": bool(run_huber_bootstrap),
            "n_huber_bootstraps": int(n_huber_bootstraps),
            "plot_mode": plot_mode,
            "shapiro_alpha": float(shapiro_alpha),
            "outlier_alpha": float(outlier_alpha),
            "mad_z_threshold": float(mad_z_threshold),
            "BA_bootstrap_iterations": int(ba_bootstrap_iterations),
            "donor_profile_bootstraps": int(donor_profile_bootstraps),
            "random_seed": int(random_seed),
            "regression_equation_position": regression_equation_position,
            "regression_metrics_position": regression_metrics_position,
        }
        with (output_dir / "APP_RUN_CONFIGURATION.json").open("w", encoding="utf-8") as handle:
            json.dump(app_config, handle, indent=2)

        zip_path = run_root / "PROXIMA_COMBINED_RESULTS.zip"
        if zip_path.exists():
            zip_path.unlink()
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in output_dir.rglob("*"):
                if path.is_file():
                    archive.write(path, path.relative_to(run_root))

    return {
        "canonical_input_path": input_csv,
        "output_dir": output_dir,
        "zip_path": zip_path,
        "combined_xlsx": combined_xlsx,
        "summary_csv": summary_csv,
        "regression": regression_result,
        "ba": ba_result,
    }


def read_binary(path: str | Path) -> bytes:
    return Path(path).read_bytes()


def collect_output_files(output_dir: str | Path, suffixes: Sequence[str] | None = None) -> dict[str, Path]:
    output_dir = Path(output_dir)
    suffix_set = {s.lower() for s in suffixes} if suffixes else None
    files = {}
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        if suffix_set is not None and path.suffix.lower() not in suffix_set:
            continue
        files[str(path.relative_to(output_dir))] = path
    return files
