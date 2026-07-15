
from __future__ import annotations

import io
import math
import zipfile
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


ALPHA_DEFAULT = 0.05

NORMAL_RANGES = {
    "RBC": (3.85, 5.65),
    "WBC": (3.60, 10.50),
    "PLT": (160.0, 370.0),
    "HCT": (35.0, 49.0),
    "HGB": (118.0, 172.0),
    "MCV": (80.0, 101.0),
    "RDW": (11.0, 16.0),
    "MCH": (27.0, 34.0),
    "MCHC": (320.0, 360.0),
    "NEUT": (1.50, 7.70),
    "LYMPH": (1.10, 4.00),
    "MXD": (0.10, 1.60),
}

MHS_SCREENSHOT_BIAS = {
    "RBC": 4.0,
    "WBC": 11.9,
    "PLT": 15.7,
    "HCT": 4.4,
    "HGB": 4.6,
    "MCV": 1.2,
    "RDW": 3.7,
    "MCH": 0.3,
    "MCHC": 0.4,
    "NEUT": 7.3,
    "LYMPH": 11.3,
    "MXD": 11.4,
}

INITIAL_EP09_BIAS = {
    "RBC": 7.0,
    "WBC": 12.5,
    "PLT": 17.5,
    "HCT": 8.6,
    "HGB": 6.0,
    "MCV": 6.6,
    "RDW": 13.5,
    "MCH": 9.5,
    "MCHC": 6.8,
    "NEUT": 11.6,
    "LYMPH": 14.0,
    "MXD": 25.0,
}

DEFAULT_COLUMN_MAP = {
    "RBC": ("RBC", "RBC_ref"),
    "WBC": ("WBC_2", "WBC_ref"),
    "PLT": ("PLT", "PLT_ref"),
    "HCT": ("HCT", "HCT_ref"),
    "HGB": ("HGB", "HGB_ref"),
    "MCV": ("MCV", "MCV_ref"),
    "RDW": ("RDW", "RDW_ref"),
    "MCH": ("MCH", "MCH_ref"),
    "MCHC": ("MCHC", "MCHC_ref"),
    "NEUT": ("NEUT_2", "NEUT_ref"),
    "LYMPH": ("LYMPH_2", "LYMPH_ref"),
    "MXD": ("MXD_2", "MXD_ref"),
}


def as_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def flag_is_true(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "1.0", "yes", "y"})
    )


def shapiro_p(values: Iterable[float]) -> float:
    v = pd.Series(values).replace([np.inf, -np.inf], np.nan).dropna().astype(float)
    if len(v) < 3 or len(v) > 5000 or v.std(ddof=1) <= 0:
        return np.nan
    return float(stats.shapiro(v).pvalue)


def mean_ci(values: Iterable[float], alpha: float = ALPHA_DEFAULT) -> dict[str, float]:
    v = pd.Series(values).replace([np.inf, -np.inf], np.nan).dropna().astype(float)
    n = int(len(v))
    result = {
        "n": n,
        "mean": np.nan,
        "sd": np.nan,
        "se": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "loa_low": np.nan,
        "loa_high": np.nan,
        "median": np.nan,
        "iqr": np.nan,
        "shapiro_p": np.nan,
    }
    if n == 0:
        return result

    result["mean"] = float(v.mean())
    result["median"] = float(v.median())
    result["iqr"] = float(v.quantile(0.75) - v.quantile(0.25))
    result["shapiro_p"] = shapiro_p(v)

    if n >= 2:
        sd = float(v.std(ddof=1))
        se = sd / math.sqrt(n)
        tcrit = float(stats.t.ppf(1 - alpha / 2, n - 1))
        result.update(
            {
                "sd": sd,
                "se": se,
                "ci_low": result["mean"] - tcrit * se,
                "ci_high": result["mean"] + tcrit * se,
                "loa_low": result["mean"] - 1.96 * sd,
                "loa_high": result["mean"] + 1.96 * sd,
            }
        )
    return result


def log_ratio_summary(
    method_a: Iterable[float],
    method_b: Iterable[float],
    alpha: float = ALPHA_DEFAULT,
) -> tuple[pd.Series, dict[str, float]]:
    a = pd.Series(method_a, dtype=float)
    b = pd.Series(method_b, dtype=float)
    valid = a.notna() & b.notna() & (a > 0) & (b > 0)
    logs = np.log(a.loc[valid] / b.loc[valid])
    raw = mean_ci(logs, alpha=alpha)

    transformed = dict(raw)
    for key in ["mean", "ci_low", "ci_high", "loa_low", "loa_high", "median"]:
        value = raw.get(key, np.nan)
        transformed[key] = (
            100.0 * (math.exp(value) - 1.0)
            if pd.notna(value) and math.isfinite(value)
            else np.nan
        )
    biases = pd.Series(np.nan, index=a.index, dtype=float)
    biases.loc[valid] = 100.0 * (np.exp(logs) - 1.0)
    return biases, transformed


def holm_adjust(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(
        [np.nan if value is None else float(value) for value in pvalues],
        dtype=float,
    )
    adjusted = np.full_like(p, np.nan, dtype=float)
    valid = np.where(np.isfinite(p))[0]
    if len(valid) == 0:
        return adjusted

    ordered = valid[np.argsort(p[valid])]
    m = len(ordered)
    running = 0.0
    for rank, index in enumerate(ordered):
        candidate = min(1.0, (m - rank) * p[index])
        running = max(running, candidate)
        adjusted[index] = running
    return adjusted


def grubbs_critical(n: int, alpha: float = ALPHA_DEFAULT) -> float:
    if n < 3:
        return np.nan
    tcrit = stats.t.ppf(1 - alpha / (2 * n), n - 2)
    return float(
        ((n - 1) / math.sqrt(n))
        * math.sqrt(tcrit**2 / (n - 2 + tcrit**2))
    )


def detect_grubbs(
    values: Iterable[float],
    alpha: float = ALPHA_DEFAULT,
) -> list[int]:
    v = pd.Series(values, dtype=float)
    valid = v.dropna()
    n = len(valid)
    if n < 7 or valid.std(ddof=1) <= 0:
        return []
    standardized = (valid - valid.mean()).abs() / valid.std(ddof=1)
    index = standardized.idxmax()
    if float(standardized.loc[index]) > grubbs_critical(n, alpha):
        return [int(index)]
    return []


def detect_mad(
    values: Iterable[float],
    threshold: float = 3.5,
    max_outliers: int = 1,
) -> list[int]:
    v = pd.Series(values, dtype=float)
    valid = v.dropna()
    if len(valid) < 4:
        return []
    median = float(valid.median())
    mad = float(stats.median_abs_deviation(valid, scale="normal", nan_policy="omit"))
    if not np.isfinite(mad) or mad <= 0:
        return []
    scores = ((valid - median).abs() / mad).sort_values(ascending=False)
    return [int(i) for i in scores[scores > threshold].index[:max_outliers]]


def generalized_esd(
    values: Iterable[float],
    max_outliers: int = 1,
    alpha: float = ALPHA_DEFAULT,
) -> list[int]:
    """Rosner generalized ESD. Returns original integer indices."""
    series = pd.Series(values, dtype=float)
    work = series.dropna().copy()
    n = len(work)
    if n < 8:
        return []

    max_outliers = max(1, min(int(max_outliers), n - 3))
    removed: list[int] = []
    statistics: list[float] = []
    critical_values: list[float] = []

    for i in range(1, max_outliers + 1):
        if len(work) < 3 or work.std(ddof=1) <= 0:
            break

        deviations = (work - work.mean()).abs()
        index = int(deviations.idxmax())
        statistic = float(deviations.loc[index] / work.std(ddof=1))

        current_n = n - i + 1
        p = 1 - alpha / (2 * current_n)
        tcrit = float(stats.t.ppf(p, current_n - 2))
        critical = float(
            ((current_n - 1) * tcrit)
            / math.sqrt((current_n - 2 + tcrit**2) * current_n)
        )

        removed.append(index)
        statistics.append(statistic)
        critical_values.append(critical)
        work = work.drop(index=index)

    count = 0
    for i, (statistic, critical) in enumerate(
        zip(statistics, critical_values), start=1
    ):
        if statistic > critical:
            count = i

    return removed[:count]


def recommend_outlier_method(values: Iterable[float]) -> str:
    v = pd.Series(values).replace([np.inf, -np.inf], np.nan).dropna()
    n = len(v)
    p = shapiro_p(v)
    if n < 7:
        return "MAD" if n >= 4 else "None"
    if np.isfinite(p) and p >= 0.05:
        return "Grubbs"
    if n >= 15:
        return "MAD (or generalized ESD if multiple outliers are expected)"
    return "MAD"



def detect_outliers_detailed(
    values: Iterable[float],
    method: str,
    alpha: float = ALPHA_DEFAULT,
    mad_threshold: float = 3.5,
    max_outliers: int = 1,
) -> list[dict[str, Any]]:
    series = pd.Series(values, dtype=float)
    method_clean = method.lower()
    if method_clean.startswith("auto"):
        method_clean = recommend_outlier_method(series).lower()

    if "grubbs" in method_clean:
        valid = series.dropna()
        n = len(valid)
        if n < 7 or valid.std(ddof=1) <= 0:
            return []
        scores = (valid - valid.mean()).abs() / valid.std(ddof=1)
        index = int(scores.idxmax())
        statistic = float(scores.loc[index])
        critical = grubbs_critical(n, alpha)
        if statistic > critical:
            return [
                {
                    "index": index,
                    "method_used": "Grubbs",
                    "statistic": statistic,
                    "critical_value": critical,
                }
            ]
        return []

    if "esd" in method_clean:
        work = series.dropna().copy()
        n = len(work)
        if n < 8:
            return []

        max_outliers = max(1, min(int(max_outliers), n - 3))
        candidates: list[dict[str, Any]] = []

        for i in range(1, max_outliers + 1):
            if len(work) < 3 or work.std(ddof=1) <= 0:
                break

            deviations = (work - work.mean()).abs()
            index = int(deviations.idxmax())
            statistic = float(deviations.loc[index] / work.std(ddof=1))

            current_n = n - i + 1
            p = 1 - alpha / (2 * current_n)
            tcrit = float(stats.t.ppf(p, current_n - 2))
            critical = float(
                ((current_n - 1) * tcrit)
                / math.sqrt((current_n - 2 + tcrit**2) * current_n)
            )
            candidates.append(
                {
                    "index": index,
                    "method_used": "Generalized ESD",
                    "statistic": statistic,
                    "critical_value": critical,
                }
            )
            work = work.drop(index=index)

        accepted_count = 0
        for i, candidate in enumerate(candidates, start=1):
            if candidate["statistic"] > candidate["critical_value"]:
                accepted_count = i
        return candidates[:accepted_count]

    if "mad" in method_clean:
        valid = series.dropna()
        if len(valid) < 4:
            return []
        median = float(valid.median())
        mad = float(
            stats.median_abs_deviation(
                valid,
                scale="normal",
                nan_policy="omit",
            )
        )
        if not np.isfinite(mad) or mad <= 0:
            return []
        scores = ((valid - median).abs() / mad).sort_values(ascending=False)
        output = []
        for index, statistic in scores.items():
            if float(statistic) <= mad_threshold:
                break
            output.append(
                {
                    "index": int(index),
                    "method_used": "MAD",
                    "statistic": float(statistic),
                    "critical_value": float(mad_threshold),
                }
            )
            if len(output) >= max_outliers:
                break
        return output

    return []


def detect_outliers(
    values: Iterable[float],
    method: str,
    alpha: float = ALPHA_DEFAULT,
    mad_threshold: float = 3.5,
    max_outliers: int = 1,
) -> list[int]:
    return [
        int(record["index"])
        for record in detect_outliers_detailed(
            values,
            method=method,
            alpha=alpha,
            mad_threshold=mad_threshold,
            max_outliers=max_outliers,
        )
    ]

def paired_tests(method_a: pd.Series, method_b: pd.Series) -> dict[str, Any]:
    frame = pd.DataFrame({"a": method_a, "b": method_b}).dropna()
    n = len(frame)
    differences = frame["a"] - frame["b"]
    normality = shapiro_p(differences)

    t_p = np.nan
    wilcoxon_p = np.nan
    if n >= 2:
        try:
            t_p = float(stats.ttest_rel(frame["a"], frame["b"]).pvalue)
        except Exception:
            pass
        try:
            wilcoxon_p = float(
                stats.wilcoxon(
                    frame["a"],
                    frame["b"],
                    zero_method="wilcox",
                    alternative="two-sided",
                ).pvalue
            )
        except Exception:
            pass

    recommended = (
        "Paired t-test"
        if np.isfinite(normality) and normality >= 0.05
        else "Wilcoxon signed-rank"
    )
    recommended_p = t_p if recommended == "Paired t-test" else wilcoxon_p
    return {
        "paired_test_n": n,
        "difference_shapiro_p": normality,
        "paired_t_p": t_p,
        "wilcoxon_p": wilcoxon_p,
        "recommended_test": recommended,
        "recommended_test_p": recommended_p,
    }


def p_f_from_ci(ci_low: float, ci_high: float, allowable_bias: float) -> str:
    if not all(np.isfinite(v) for v in [ci_low, ci_high, allowable_bias]):
        return "NOT ESTIMABLE"
    return (
        "PASS"
        if ci_low >= -allowable_bias and ci_high <= allowable_bias
        else "FAIL"
    )


def mean_bias_equivalence_n(
    mean_bias: float,
    sd_bias: float,
    allowable_bias: float,
    power: float = 0.80,
    alpha: float = ALPHA_DEFAULT,
) -> float:
    if not all(
        np.isfinite(v)
        for v in [mean_bias, sd_bias, allowable_bias, power, alpha]
    ):
        return np.nan
    margin = allowable_bias - abs(mean_bias)
    if margin <= 0 or sd_bias <= 0:
        return np.inf
    z_alpha = stats.norm.ppf(1 - alpha)
    z_power = stats.norm.ppf(power)
    return float(math.ceil(((z_alpha + z_power) * sd_bias / margin) ** 2))


def mean_bias_ci_precision_n(
    sd_bias: float,
    target_halfwidth: float,
    alpha: float = ALPHA_DEFAULT,
) -> float:
    if not np.isfinite(sd_bias) or not np.isfinite(target_halfwidth):
        return np.nan
    if sd_bias <= 0 or target_halfwidth <= 0:
        return np.inf
    z = stats.norm.ppf(1 - alpha / 2)
    return float(math.ceil((z * sd_bias / target_halfwidth) ** 2))


def loa_ci_halfwidth(
    sd_bias: float,
    n: int,
    alpha: float = ALPHA_DEFAULT,
) -> float:
    if not np.isfinite(sd_bias) or sd_bias <= 0 or n < 3:
        return np.nan
    tcrit = stats.t.ppf(1 - alpha / 2, n - 1)
    return float(
        tcrit
        * sd_bias
        * math.sqrt(1 / n + (1.96**2) / (2 * (n - 1)))
    )


def loa_precision_n(
    sd_bias: float,
    target_halfwidth: float,
    alpha: float = ALPHA_DEFAULT,
    max_n: int = 20000,
) -> float:
    if not np.isfinite(sd_bias) or not np.isfinite(target_halfwidth):
        return np.nan
    if sd_bias <= 0 or target_halfwidth <= 0:
        return np.inf
    for n in range(3, max_n + 1):
        if loa_ci_halfwidth(sd_bias, n, alpha) <= target_halfwidth:
            return float(n)
    return np.inf


def median_log_recalibrated_bias(
    method_a: pd.Series,
    method_b: pd.Series,
    leave_one_out: bool = True,
) -> tuple[pd.Series, float]:
    a = pd.Series(method_a, dtype=float)
    b = pd.Series(method_b, dtype=float)
    valid = a.notna() & b.notna() & (a > 0) & (b > 0)
    ratios = (a.loc[valid] / b.loc[valid]).astype(float)
    output = pd.Series(np.nan, index=a.index, dtype=float)

    if len(ratios) == 0:
        return output, np.nan

    overall_factor = float(np.exp(np.median(np.log(ratios))))

    for index in ratios.index:
        if leave_one_out and len(ratios) >= 4:
            training = ratios.drop(index=index)
            factor = float(np.exp(np.median(np.log(training))))
        else:
            factor = overall_factor
        output.loc[index] = 100.0 * ((ratios.loc[index] / factor) - 1.0)

    return output, overall_factor


def _extract_bias_method(
    pairs: pd.DataFrame,
    method: str,
    leave_one_out_median: bool = True,
) -> tuple[pd.Series, dict[str, Any]]:
    a = pairs["method_a"]
    b = pairs["method_b"]
    ref_a = pairs.get("ref_a", pd.Series(np.nan, index=pairs.index))
    ref_b = pairs.get("ref_b", pd.Series(np.nan, index=pairs.index))

    metadata: dict[str, Any] = {}

    if method == "Raw percent BA":
        return 100.0 * (a - b) / b.replace(0, np.nan), metadata

    if method == "Direct log ratio":
        valid = (a > 0) & (b > 0)
        output = pd.Series(np.nan, index=pairs.index, dtype=float)
        output.loc[valid] = 100.0 * (
            np.exp(np.log(a.loc[valid] / b.loc[valid])) - 1.0
        )
        return output, metadata

    if method == "Median log recalibration":
        output, factor = median_log_recalibrated_bias(
            a,
            b,
            leave_one_out=leave_one_out_median,
        )
        metadata["median_ratio_factor"] = factor
        return output, metadata

    if method == "Sysmex absolute correction":
        output = 100.0 * (((a - b) - (ref_a - ref_b)) / b.replace(0, np.nan))
        return output, metadata

    if method == "Sysmex log-ratio correction":
        valid = (a > 0) & (b > 0) & (ref_a > 0) & (ref_b > 0)
        output = pd.Series(np.nan, index=pairs.index, dtype=float)
        output.loc[valid] = 100.0 * (
            np.exp(
                np.log(a.loc[valid] / b.loc[valid])
                - np.log(ref_a.loc[valid] / ref_b.loc[valid])
            )
            - 1.0
        )
        return output, metadata

    raise ValueError(f"Unknown method: {method}")


def _preaverage_outlier_filter(
    data: pd.DataFrame,
    donor_col: str,
    type_col: str,
    value_col: str,
    arm_labels: tuple[str, str],
    method: str,
    alpha: float,
    mad_threshold: float,
    max_outliers: int,
    sample_id_col: str | None = None,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    filtered = data.copy()
    filtered["_keep_for_analysis"] = True
    log: list[dict[str, Any]] = []

    for arm in arm_labels:
        arm_mask = filtered[type_col].astype(str).eq(str(arm))
        arm_data = filtered.loc[arm_mask].copy()
        if arm_data.empty:
            continue

        donor_median = arm_data.groupby(donor_col)[value_col].transform("median")
        denominator = donor_median.replace(0, np.nan)
        residual = 100.0 * (arm_data[value_col] - donor_median) / denominator

        details = detect_outliers_detailed(
            residual.reset_index(drop=True),
            method=method,
            alpha=alpha,
            mad_threshold=mad_threshold,
            max_outliers=max_outliers,
        )
        for detail in details:
            position = int(detail["index"])
            real_index = arm_data.index[position]
            filtered.loc[real_index, "_keep_for_analysis"] = False
            log.append(
                {
                    "stage": "Pre-averaging replicate residual",
                    "arm": arm,
                    "donor": filtered.loc[real_index, donor_col],
                    "sample_id": (
                        filtered.loc[real_index, sample_id_col]
                        if sample_id_col and sample_id_col in filtered.columns
                        else real_index
                    ),
                    "row_index": real_index,
                    "value": filtered.loc[real_index, value_col],
                    "outlier_metric": residual.iloc[position],
                    "method": detail["method_used"],
                    "outlier_statistic": detail["statistic"],
                    "critical_value_or_threshold": detail["critical_value"],
                }
            )

    return filtered.loc[filtered["_keep_for_analysis"]].copy(), log


def prepare_long_pairs(
    data: pd.DataFrame,
    *,
    donor_col: str,
    type_col: str,
    method_a_label: str,
    method_b_label: str,
    value_col: str,
    reference_col: str | None,
    normal_low: float | None,
    normal_high: float | None,
    filter_to_normal: bool,
    require_reference_normal: bool,
    global_flag_col: str | None,
    sample_id_col: str | None,
    outlier_enabled: bool,
    outlier_stage: str,
    outlier_method: str,
    outlier_alpha: float,
    mad_threshold: float,
    max_outliers: int,
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, Any]]:
    columns = [donor_col, type_col, value_col]
    if reference_col and reference_col in data.columns:
        columns.append(reference_col)
    if global_flag_col and global_flag_col in data.columns:
        columns.append(global_flag_col)
    if sample_id_col and sample_id_col in data.columns:
        columns.append(sample_id_col)

    d = data[columns].copy()
    d[value_col] = as_numeric(d[value_col])
    if reference_col and reference_col in d.columns:
        d[reference_col] = as_numeric(d[reference_col])

    d = d[
        d[type_col].astype(str).isin(
            [str(method_a_label), str(method_b_label)]
        )
    ].copy()
    if global_flag_col and global_flag_col in d.columns:
        d = d.loc[~flag_is_true(d[global_flag_col])].copy()

    rows_initial = len(d)
    d = d.dropna(subset=[donor_col, type_col, value_col]).copy()

    if filter_to_normal and normal_low is not None and normal_high is not None:
        d = d[d[value_col].between(normal_low, normal_high, inclusive="both")].copy()
        if (
            require_reference_normal
            and reference_col
            and reference_col in d.columns
        ):
            d = d[
                d[reference_col].between(
                    normal_low,
                    normal_high,
                    inclusive="both",
                )
            ].copy()

    outlier_log: list[dict[str, Any]] = []
    if outlier_enabled and outlier_stage == "Pre-averaging replicate residual":
        d, outlier_log = _preaverage_outlier_filter(
            d,
            donor_col=donor_col,
            type_col=type_col,
            value_col=value_col,
            arm_labels=(method_a_label, method_b_label),
            method=outlier_method,
            alpha=outlier_alpha,
            mad_threshold=mad_threshold,
            max_outliers=max_outliers,
            sample_id_col=sample_id_col,
        )

    aggregations: dict[str, tuple[str, str]] = {
        "value": (value_col, "mean"),
        "n_replicates": (value_col, "count"),
    }
    if reference_col and reference_col in d.columns:
        aggregations["reference"] = (reference_col, "mean")

    grouped = (
        d.groupby([donor_col, type_col], as_index=False)
        .agg(**aggregations)
    )

    value_wide = grouped.pivot(
        index=donor_col,
        columns=type_col,
        values="value",
    )
    count_wide = grouped.pivot(
        index=donor_col,
        columns=type_col,
        values="n_replicates",
    )

    pairs = pd.DataFrame(index=value_wide.index)
    pairs["donor"] = value_wide.index
    pairs["method_a"] = value_wide.get(method_a_label)
    pairs["method_b"] = value_wide.get(method_b_label)
    pairs["n_a_replicates"] = count_wide.get(method_a_label)
    pairs["n_b_replicates"] = count_wide.get(method_b_label)

    if "reference" in grouped.columns:
        ref_wide = grouped.pivot(
            index=donor_col,
            columns=type_col,
            values="reference",
        )
        pairs["ref_a"] = ref_wide.get(method_a_label)
        pairs["ref_b"] = ref_wide.get(method_b_label)
    else:
        pairs["ref_a"] = np.nan
        pairs["ref_b"] = np.nan

    pairs = pairs.reset_index(drop=True)
    pairs = pairs.dropna(subset=["method_a", "method_b"]).copy()

    audit = {
        "rows_initial": rows_initial,
        "rows_after_filters": len(d),
        "paired_donors_before_pair_outlier": len(pairs),
        "normal_filter_applied": filter_to_normal,
        "reference_normal_required": require_reference_normal,
    }
    return pairs, outlier_log, audit


def prepare_wide_pairs(
    data: pd.DataFrame,
    *,
    donor_col: str | None,
    method_a_col: str,
    method_b_col: str,
    ref_a_col: str | None,
    ref_b_col: str | None,
    normal_low: float | None,
    normal_high: float | None,
    filter_to_normal: bool,
    require_reference_normal: bool,
    global_flag_col: str | None,
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, Any]]:
    d = data.copy()

    if global_flag_col and global_flag_col in d.columns:
        d = d.loc[~flag_is_true(d[global_flag_col])].copy()

    pairs = pd.DataFrame(index=d.index)
    pairs["donor"] = (
        d[donor_col].astype(str)
        if donor_col and donor_col in d.columns
        else d.index.astype(str)
    )
    pairs["method_a"] = as_numeric(d[method_a_col])
    pairs["method_b"] = as_numeric(d[method_b_col])
    pairs["ref_a"] = (
        as_numeric(d[ref_a_col])
        if ref_a_col and ref_a_col in d.columns
        else np.nan
    )
    pairs["ref_b"] = (
        as_numeric(d[ref_b_col])
        if ref_b_col and ref_b_col in d.columns
        else np.nan
    )
    pairs["n_a_replicates"] = 1
    pairs["n_b_replicates"] = 1
    pairs = pairs.dropna(subset=["method_a", "method_b"]).copy()

    if filter_to_normal and normal_low is not None and normal_high is not None:
        mask = (
            pairs["method_a"].between(normal_low, normal_high, inclusive="both")
            & pairs["method_b"].between(normal_low, normal_high, inclusive="both")
        )
        if require_reference_normal and pairs[["ref_a", "ref_b"]].notna().all(axis=1).any():
            mask &= (
                pairs["ref_a"].between(normal_low, normal_high, inclusive="both")
                & pairs["ref_b"].between(normal_low, normal_high, inclusive="both")
            )
        pairs = pairs.loc[mask].copy()

    audit = {
        "rows_initial": len(data),
        "rows_after_filters": len(pairs),
        "paired_donors_before_pair_outlier": len(pairs),
        "normal_filter_applied": filter_to_normal,
        "reference_normal_required": require_reference_normal,
    }
    return pairs.reset_index(drop=True), [], audit



def analyze_pairs(
    pairs: pd.DataFrame,
    *,
    analyte: str,
    allowable_bias: float,
    normal_low: float | None,
    normal_high: float | None,
    outlier_enabled: bool,
    outlier_stage: str,
    outlier_method: str,
    outlier_alpha: float,
    mad_threshold: float,
    max_outliers: int,
    include_direct_log: bool,
    include_median_log: bool,
    median_log_loo: bool,
    include_sysmex_absolute: bool,
    include_sysmex_log: bool,
    power_enabled: bool,
    power: float,
    sd_distance_enabled: bool,
) -> dict[str, Any]:
    original = pairs.copy().reset_index(drop=True)
    original["raw_percent_bias"] = (
        100.0
        * (original["method_a"] - original["method_b"])
        / original["method_b"].replace(0, np.nan)
    )
    original["analysis_included"] = True
    original["outlier_flag"] = False

    outlier_log: list[dict[str, Any]] = []
    outlier_indices: list[int] = []
    if outlier_enabled and outlier_stage == "Donor-pair BA bias":
        details = detect_outliers_detailed(
            original["raw_percent_bias"],
            method=outlier_method,
            alpha=outlier_alpha,
            mad_threshold=mad_threshold,
            max_outliers=max_outliers,
        )
        outlier_indices = [int(detail["index"]) for detail in details]
        for detail in details:
            index = int(detail["index"])
            original.loc[index, "analysis_included"] = False
            original.loc[index, "outlier_flag"] = True
            outlier_log.append(
                {
                    "stage": "Donor-pair BA bias",
                    "arm": "paired",
                    "donor": original.loc[index, "donor"],
                    "sample_id": original.loc[index, "donor"],
                    "row_index": index,
                    "value": original.loc[index, "raw_percent_bias"],
                    "outlier_metric": original.loc[index, "raw_percent_bias"],
                    "method": detail["method_used"],
                    "outlier_statistic": detail["statistic"],
                    "critical_value_or_threshold": detail["critical_value"],
                }
            )

    work = original.loc[original["analysis_included"]].copy().reset_index(drop=True)

    raw_diff = work["method_a"] - work["method_b"]
    raw_pct = (
        100.0
        * raw_diff
        / work["method_b"].replace(0, np.nan)
    )

    shapiro_diff = shapiro_p(raw_diff)
    shapiro_pct = shapiro_p(raw_pct)
    try:
        levene_p = float(
            stats.levene(
                work["method_a"].dropna(),
                work["method_b"].dropna(),
                center="median",
            ).pvalue
        )
    except Exception:
        levene_p = np.nan

    tests = paired_tests(work["method_a"], work["method_b"])

    recommendation_parts = []
    if np.isfinite(shapiro_diff) and shapiro_diff < 0.05:
        recommendation_parts.append(
            "Paired differences are non-normal; report Wilcoxon and inspect a log-ratio or median-based sensitivity analysis."
        )
    else:
        recommendation_parts.append(
            "Paired differences are compatible with normality; the t-based mean-bias CI is reasonable."
        )

    if (work["method_a"] <= 0).any() or (work["method_b"] <= 0).any():
        recommendation_parts.append(
            "Log-ratio methods are unavailable for zero or negative measurements."
        )
    else:
        log_difference = np.log(work["method_a"] / work["method_b"])
        log_p = shapiro_p(log_difference)
        if np.isfinite(log_p) and (
            not np.isfinite(shapiro_pct) or log_p > shapiro_pct
        ):
            recommendation_parts.append(
                "The log-ratio scale improves normality relative to raw percent bias."
            )

    recommendation_parts.append(
        f"Suggested outlier method: {recommend_outlier_method(raw_pct)}."
    )

    methods = ["Raw percent BA"]
    if include_direct_log:
        methods.append("Direct log ratio")
    if include_median_log:
        methods.append("Median log recalibration")

    refs_available = work[["ref_a", "ref_b"]].notna().all(axis=1).sum() >= 2
    if include_sysmex_absolute and refs_available:
        methods.append("Sysmex absolute correction")
    if include_sysmex_log and refs_available:
        methods.append("Sysmex log-ratio correction")

    summary_rows: list[dict[str, Any]] = []
    method_biases_analysis: dict[str, pd.Series] = {}
    method_biases_plot: dict[str, pd.Series] = {}

    for method in methods:
        analysis_biases, metadata = _extract_bias_method(
            work,
            method,
            leave_one_out_median=median_log_loo,
        )
        method_biases_analysis[method] = analysis_biases

        if method == "Direct log ratio":
            valid = (
                work["method_a"].notna()
                & work["method_b"].notna()
                & (work["method_a"] > 0)
                & (work["method_b"] > 0)
            )
            log_values = np.log(
                work.loc[valid, "method_a"]
                / work.loc[valid, "method_b"]
            )
            log_stats = mean_ci(log_values)
            percent_stats = mean_ci(analysis_biases)
            stats_summary = dict(percent_stats)
            for key in ["mean", "ci_low", "ci_high", "loa_low", "loa_high", "median"]:
                value = log_stats.get(key, np.nan)
                stats_summary[key] = (
                    100.0 * (math.exp(value) - 1.0)
                    if pd.notna(value) and math.isfinite(value)
                    else np.nan
                )
            stats_summary["n"] = log_stats["n"]
            stats_summary["shapiro_p"] = log_stats["shapiro_p"]

        elif method == "Sysmex log-ratio correction":
            valid = (
                work["method_a"].notna()
                & work["method_b"].notna()
                & work["ref_a"].notna()
                & work["ref_b"].notna()
                & (work["method_a"] > 0)
                & (work["method_b"] > 0)
                & (work["ref_a"] > 0)
                & (work["ref_b"] > 0)
            )
            log_values = (
                np.log(
                    work.loc[valid, "method_a"]
                    / work.loc[valid, "method_b"]
                )
                - np.log(
                    work.loc[valid, "ref_a"]
                    / work.loc[valid, "ref_b"]
                )
            )
            log_stats = mean_ci(log_values)
            percent_stats = mean_ci(analysis_biases)
            stats_summary = dict(percent_stats)
            for key in ["mean", "ci_low", "ci_high", "loa_low", "loa_high", "median"]:
                value = log_stats.get(key, np.nan)
                stats_summary[key] = (
                    100.0 * (math.exp(value) - 1.0)
                    if pd.notna(value) and math.isfinite(value)
                    else np.nan
                )
            stats_summary["n"] = log_stats["n"]
            stats_summary["shapiro_p"] = log_stats["shapiro_p"]

        else:
            stats_summary = mean_ci(analysis_biases)

        # For visualization, retain excluded donor-pair outliers. This is descriptive:
        # the summary statistics above remain based only on analysis-included pairs.
        plot_biases, _ = _extract_bias_method(
            original,
            method,
            leave_one_out_median=median_log_loo,
        )
        method_biases_plot[method] = plot_biases

        target_halfwidth = allowable_bias
        equivalence_n = (
            mean_bias_equivalence_n(
                stats_summary["mean"],
                stats_summary["sd"],
                allowable_bias,
                power=power,
            )
            if power_enabled
            else np.nan
        )
        ci_precision_n = (
            mean_bias_ci_precision_n(
                stats_summary["sd"],
                target_halfwidth,
            )
            if power_enabled
            else np.nan
        )
        loa_n = (
            loa_precision_n(
                stats_summary["sd"],
                target_halfwidth,
            )
            if power_enabled
            else np.nan
        )

        summary_rows.append(
            {
                "Analyte": analyte,
                "Method": method,
                "n paired": stats_summary["n"],
                "Mean BA bias, %": stats_summary["mean"],
                "Median BA bias, %": stats_summary["median"],
                "BA bias SD, %": stats_summary["sd"],
                "BA 95% CI low, %": stats_summary["ci_low"],
                "BA 95% CI high, %": stats_summary["ci_high"],
                "BA LoA low, %": stats_summary["loa_low"],
                "BA LoA high, %": stats_summary["loa_high"],
                "Bias Shapiro p": stats_summary["shapiro_p"],
                "Allowable trueness bias, %": allowable_bias,
                "P/F": p_f_from_ci(
                    stats_summary["ci_low"],
                    stats_summary["ci_high"],
                    allowable_bias,
                ),
                "Mean-bias equivalence n": equivalence_n,
                "Mean-bias CI precision n": ci_precision_n,
                "LoA precision n": loa_n,
                "Median ratio factor":
                    metadata.get("median_ratio_factor", np.nan),
            }
        )

    summary = pd.DataFrame(summary_rows)

    pair_output = original.copy()
    for method, biases in method_biases_plot.items():
        pair_output[f"{method} bias, %"] = biases.values

    pair_output["pair_average"] = (
        pair_output["method_a"] + pair_output["method_b"]
    ) / 2.0
    pair_output["raw_difference"] = (
        pair_output["method_a"] - pair_output["method_b"]
    )

    sd_table = pd.DataFrame()
    if (
        sd_distance_enabled
        and normal_low is not None
        and normal_high is not None
        and np.isfinite(normal_low)
        and np.isfinite(normal_high)
        and normal_high > normal_low
    ):
        included_pairs = pair_output.loc[pair_output["analysis_included"]].copy()
        source = (
            (included_pairs["ref_a"] + included_pairs["ref_b"]) / 2.0
            if included_pairs[["ref_a", "ref_b"]].notna().all(axis=1).sum() >= 2
            else included_pairs["pair_average"]
        )
        midpoint = (normal_low + normal_high) / 2.0
        sd_proxy = (normal_high - normal_low) / 4.0
        signed = (source - midpoint) / sd_proxy
        absolute = signed.abs()

        included_pairs["SD-distance source"] = source
        included_pairs["Signed SD distance"] = signed
        included_pairs["Absolute SD distance"] = absolute

        labels = ["<0.5 SD", "0.5–1 SD", "1–2 SD", ">2 SD"]
        bins = [-np.inf, 0.5, 1.0, 2.0, np.inf]
        included_pairs["SD distance bin"] = pd.cut(
            absolute,
            bins=bins,
            labels=labels,
            right=False,
        )

        raw_bias_col = "Raw percent BA bias, %"

        def safe_mean_abs(series: pd.Series) -> float:
            values = pd.to_numeric(series, errors="coerce").dropna()
            return float(np.abs(values).mean()) if len(values) else np.nan

        sd_table = (
            included_pairs.groupby("SD distance bin", observed=False)
            .agg(
                n=("donor", "count"),
                mean_absolute_bias_pct=(raw_bias_col, safe_mean_abs),
                mean_signed_bias_pct=(raw_bias_col, "mean"),
                min_signed_sd=("Signed SD distance", "min"),
                max_signed_sd=("Signed SD distance", "max"),
            )
            .reset_index()
        )
        sd_table.insert(0, "Analyte", analyte)
        sd_table["Estimated SD"] = sd_proxy

        # Add SD-distance columns back to the retained rows in the full pair output.
        for column in [
            "SD-distance source",
            "Signed SD distance",
            "Absolute SD distance",
            "SD distance bin",
        ]:
            pair_output.loc[
                pair_output["analysis_included"], column
            ] = included_pairs[column].to_numpy()

    diagnostics = {
        "Analyte": analyte,
        "n before donor-pair outlier removal": len(pairs),
        "n after outlier removal": len(work),
        "Raw difference Shapiro p": shapiro_diff,
        "Raw percent-bias Shapiro p": shapiro_pct,
        "Levene p": levene_p,
        "Recommended outlier method": recommend_outlier_method(raw_pct),
        "Recommendation": " ".join(recommendation_parts),
        **tests,
    }

    return {
        "summary": summary,
        "pairs": pair_output,
        "diagnostics": pd.DataFrame([diagnostics]),
        "outliers": pd.DataFrame(outlier_log),
        "sd_distance": sd_table,
        "method_biases": method_biases_analysis,
    }

def plot_ba_percent(
    pairs: pd.DataFrame,
    *,
    analyte: str,
    method: str,
    allowable_bias: float,
    outlier_donors: set[str] | None = None,
) -> plt.Figure:
    bias_col = f"{method} bias, %"
    if bias_col not in pairs.columns:
        raise KeyError(bias_col)

    columns = ["donor", "pair_average", bias_col]
    if "analysis_included" in pairs.columns:
        columns.append("analysis_included")
    frame = pairs[columns].dropna(subset=["donor", "pair_average", bias_col]).copy()
    included = (
        frame["analysis_included"].astype(bool)
        if "analysis_included" in frame.columns
        else pd.Series(True, index=frame.index)
    )
    stats_summary = mean_ci(frame.loc[included, bias_col])
    outlier_donors = outlier_donors or set()

    fig, ax = plt.subplots(figsize=(9.2, 5.8), dpi=150)
    regular = ~frame["donor"].astype(str).isin(outlier_donors)
    ax.scatter(
        frame.loc[regular, "pair_average"],
        frame.loc[regular, bias_col],
        s=42,
        alpha=0.85,
        label="Paired donor",
    )
    if (~regular).any():
        ax.scatter(
            frame.loc[~regular, "pair_average"],
            frame.loc[~regular, bias_col],
            s=70,
            marker="X",
            label="Detected outlier",
        )

    ax.axhline(
        stats_summary["mean"],
        linewidth=1.5,
        label=f"Mean bias {stats_summary['mean']:.2f}%",
    )
    ax.axhline(
        stats_summary["ci_low"],
        linestyle="--",
        linewidth=1.3,
        label=f"95% CI low {stats_summary['ci_low']:.2f}%",
    )
    ax.axhline(
        stats_summary["ci_high"],
        linestyle="--",
        linewidth=1.3,
        label=f"95% CI high {stats_summary['ci_high']:.2f}%",
    )
    if np.isfinite(allowable_bias):
        ax.axhline(
            allowable_bias,
            linestyle=":",
            linewidth=1.4,
            label=f"+allowable bias {allowable_bias:g}%",
        )
        ax.axhline(-allowable_bias, linestyle=":", linewidth=1.4)

    for _, row in frame.iterrows():
        if str(row["donor"]) in outlier_donors:
            ax.annotate(
                str(row["donor"]),
                (row["pair_average"], row[bias_col]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
            )

    ax.axhline(0, linewidth=0.9)
    ax.set_title(f"{analyte} Bland–Altman percent bias — {method}")
    ax.set_xlabel("Average of paired methods")
    ax.set_ylabel("Percent bias: 100 × (Method A − Method B) / Method B")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return fig


def plot_forest(summary: pd.DataFrame, method: str) -> plt.Figure:
    subset = summary[summary["Method"] == method].copy()
    subset = subset.sort_values("Analyte").reset_index(drop=True)

    fig, ax = plt.subplots(
        figsize=(10.5, max(4.5, 0.55 * len(subset) + 2.0)),
        dpi=150,
    )
    y = np.arange(len(subset))

    mean = subset["Mean BA bias, %"].to_numpy(float)
    low = subset["BA 95% CI low, %"].to_numpy(float)
    high = subset["BA 95% CI high, %"].to_numpy(float)

    for i, row in subset.iterrows():
        allowable = row["Allowable trueness bias, %"]
        ax.hlines(
            i,
            -allowable,
            allowable,
            linewidth=5,
            alpha=0.35,
        )

    ax.errorbar(
        mean,
        y,
        xerr=np.vstack([mean - low, high - mean]),
        fmt="s",
        capsize=4,
        label="BA mean bias with 95% CI",
    )
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(subset["Analyte"])
    ax.invert_yaxis()
    ax.set_xlabel("BA percent bias with 95% CI")
    ax.set_title(f"Interchangeability summary — {method}")
    ax.grid(axis="x", alpha=0.25)

    for i, status in enumerate(subset["P/F"]):
        ax.text(
            1.01,
            i,
            status,
            transform=ax.get_yaxis_transform(),
            va="center",
            fontweight="bold",
        )

    ax.legend(loc="best")
    fig.tight_layout()
    return fig


def figure_to_png_bytes(fig: plt.Figure) -> bytes:
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=180, bbox_inches="tight")
    buffer.seek(0)
    return buffer.getvalue()


def dataframe_to_excel_bytes(sheets: dict[str, pd.DataFrame]) -> bytes:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
        for name, frame in sheets.items():
            safe_name = name[:31]
            frame.to_excel(writer, sheet_name=safe_name, index=False)
            workbook = writer.book
            worksheet = writer.sheets[safe_name]
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "bg_color": "#D9EAF7",
                    "border": 1,
                    "text_wrap": True,
                    "valign": "top",
                }
            )
            for col_num, value in enumerate(frame.columns.values):
                worksheet.write(0, col_num, value, header_format)
                width = min(max(len(str(value)) + 2, 12), 34)
                worksheet.set_column(col_num, col_num, width)
            worksheet.freeze_panes(1, 0)
    buffer.seek(0)
    return buffer.getvalue()


def build_results_zip(
    *,
    summary: pd.DataFrame,
    diagnostics: pd.DataFrame,
    pairs: pd.DataFrame,
    outliers: pd.DataFrame,
    sd_distance: pd.DataFrame,
    plots: dict[str, bytes],
    settings_text: str,
) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(
            "summary.csv",
            summary.to_csv(index=False, encoding="utf-8-sig"),
        )
        archive.writestr(
            "diagnostics.csv",
            diagnostics.to_csv(index=False, encoding="utf-8-sig"),
        )
        archive.writestr(
            "paired_points.csv",
            pairs.to_csv(index=False, encoding="utf-8-sig"),
        )
        archive.writestr(
            "outliers.csv",
            outliers.to_csv(index=False, encoding="utf-8-sig"),
        )
        archive.writestr(
            "sd_distance.csv",
            sd_distance.to_csv(index=False, encoding="utf-8-sig"),
        )
        archive.writestr("analysis_settings.txt", settings_text)

        excel_bytes = dataframe_to_excel_bytes(
            {
                "Summary": summary,
                "Diagnostics": diagnostics,
                "Paired points": pairs,
                "Outliers": outliers,
                "SD distance": sd_distance,
            }
        )
        archive.writestr("interchangeability_results.xlsx", excel_bytes)

        for filename, payload in plots.items():
            archive.writestr(f"plots/{filename}", payload)

    buffer.seek(0)
    return buffer.getvalue()
