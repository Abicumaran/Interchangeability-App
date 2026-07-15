
from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from analysis_core import (
    DEFAULT_COLUMN_MAP,
    INITIAL_EP09_BIAS,
    MHS_SCREENSHOT_BIAS,
    NORMAL_RANGES,
    analyze_pairs,
    build_results_zip,
    dataframe_to_excel_bytes,
    figure_to_png_bytes,
    flag_is_true,
    holm_adjust,
    plot_ba_percent,
    plot_forest,
    prepare_long_pairs,
    prepare_wide_pairs,
)


st.set_page_config(
    page_title="Interchangeability App",
    page_icon="🧪",
    layout="wide",
)

st.title("🧪 Interchangeability App")
st.caption(
    "Paired Bland–Altman percent-bias analysis, trueness P/F, normality-guided "
    "sensitivity analyses, optional outlier handling, Sysmex correction, "
    "paired tests, power, and SD-distance coverage."
)

with st.expander("Important use note", expanded=False):
    st.warning(
        "Use de-identified data only. The app provides statistical screening and "
        "does not replace a pre-specified validation protocol, clinical decision, "
        "or regulatory review. Outlier removal and calibration should be justified "
        "independently and reported transparently."
    )

uploaded = st.file_uploader(
    "1. Upload a CSV or Excel datasheet",
    type=["csv", "xlsx", "xls"],
)

if uploaded is None:
    st.info(
        "Upload a file to begin. Long format such as Donor + Type + analyte columns "
        "and wide format with explicit paired columns are both supported."
    )
    st.stop()


@st.cache_data(show_spinner=False)
def read_upload(file_name: str, file_bytes: bytes):
    suffix = Path(file_name).suffix.lower()
    if suffix == ".csv":
        return {"CSV": pd.read_csv(io.BytesIO(file_bytes))}
    excel = pd.ExcelFile(io.BytesIO(file_bytes))
    return {
        sheet: pd.read_excel(io.BytesIO(file_bytes), sheet_name=sheet)
        for sheet in excel.sheet_names
    }


file_bytes = uploaded.getvalue()
sheets = read_upload(uploaded.name, file_bytes)
sheet_name = st.selectbox("Sheet", list(sheets), index=0)
data = sheets[sheet_name].copy()

st.success(f"Loaded {len(data):,} rows and {len(data.columns):,} columns.")
with st.expander("Preview uploaded data", expanded=False):
    st.dataframe(data.head(30), use_container_width=True)

columns = list(data.columns)

# Dataset compatibility check for the EP09/APT long-format workbooks.
detected_value_columns = {
    analyte: candidate
    for analyte, (candidate, _) in DEFAULT_COLUMN_MAP.items()
    if candidate in columns
}
detected_reference_columns = {
    analyte: reference
    for analyte, (_, reference) in DEFAULT_COLUMN_MAP.items()
    if reference in columns
}
with st.expander("Automatic dataset compatibility check", expanded=True):
    check1, check2, check3, check4 = st.columns(4)
    check1.metric("Recognized analytes", f"{len(detected_value_columns)}/12")
    check2.metric("Recognized references", f"{len(detected_reference_columns)}/12")
    check3.metric("Has Donor + Type", "Yes" if {"Donor", "Type"}.issubset(columns) else "No")
    check4.metric("Has global_flag", "Yes" if "global_flag" in columns else "No")
    if {"Donor", "Type"}.issubset(columns):
        donor_labels = data["Donor"].dropna().astype(str).drop_duplicates().tolist()
        type_labels = data["Type"].dropna().astype(str).drop_duplicates().tolist()
        st.write(
            f"Detected **{len(donor_labels)} exact donor labels** and Type labels: "
            + ", ".join(type_labels)
        )
        if any(label.lower().endswith("b") for label in donor_labels):
            st.caption(
                "Exact donor strings are preserved, so labels such as D04 and D04b "
                "remain separate during donor averaging and pairing."
            )
    missing_values = sorted(set(DEFAULT_COLUMN_MAP) - set(detected_value_columns))
    missing_refs = sorted(set(DEFAULT_COLUMN_MAP) - set(detected_reference_columns))
    if missing_values:
        st.warning("Unrecognized default analyte columns: " + ", ".join(missing_values))
    if missing_refs:
        st.info(
            "Missing default reference columns: " + ", ".join(missing_refs)
            + ". Sysmex corrections remain optional and unavailable only for those analytes."
        )

layout = st.radio(
    "2. Data layout",
    ["Long format: Donor + Type", "Wide format: paired columns"],
    horizontal=True,
)

st.subheader("3. Pairing and analyte mapping")

if layout.startswith("Long"):
    col1, col2, col3 = st.columns(3)
    donor_default = columns.index("Donor") if "Donor" in columns else 0
    type_default = columns.index("Type") if "Type" in columns else min(1, len(columns) - 1)
    donor_col = col1.selectbox("Donor/subject column", columns, index=donor_default)
    type_col = col2.selectbox("Method/type column", columns, index=type_default)

    unique_types = (
        data[type_col].dropna().astype(str).drop_duplicates().tolist()
        if type_col in data.columns
        else []
    )
    if len(unique_types) < 2:
        st.error("The selected Type column needs at least two labels.")
        st.stop()

    method_a_default = unique_types.index("Mc-") if "Mc-" in unique_types else 0
    method_a_label = col3.selectbox(
        "Method A / capillary label",
        unique_types,
        index=method_a_default,
    )
    method_b_options = [value for value in unique_types if value != method_a_label]
    method_b_default = method_b_options.index("MK2") if "MK2" in method_b_options else 0
    method_b_label = col3.selectbox(
        "Method B / venous label",
        method_b_options,
        index=method_b_default,
    )

    optional_cols = ["None"] + columns
    sample_default = (
        optional_cols.index("bloodSampleId")
        if "bloodSampleId" in optional_cols
        else 0
    )
    flag_default = (
        optional_cols.index("global_flag")
        if "global_flag" in optional_cols
        else 0
    )
    col4, col5 = st.columns(2)
    sample_id_col_ui = col4.selectbox(
        "Optional sample/replicate ID",
        optional_cols,
        index=sample_default,
    )
    global_flag_col_ui = col5.selectbox(
        "Optional global flag column",
        optional_cols,
        index=flag_default,
    )
    sample_id_col = None if sample_id_col_ui == "None" else sample_id_col_ui
    global_flag_col = None if global_flag_col_ui == "None" else global_flag_col_ui
else:
    donor_options = ["Use row index"] + columns
    donor_col_ui = st.selectbox(
        "Optional donor/subject column",
        donor_options,
        index=(donor_options.index("Donor") if "Donor" in donor_options else 0),
    )
    donor_col = None if donor_col_ui == "Use row index" else donor_col_ui
    optional_cols = ["None"] + columns
    global_flag_col_ui = st.selectbox(
        "Optional global flag column",
        optional_cols,
        index=(
            optional_cols.index("global_flag")
            if "global_flag" in optional_cols
            else 0
        ),
    )
    global_flag_col = None if global_flag_col_ui == "None" else global_flag_col_ui
    type_col = None
    method_a_label = "Method A"
    method_b_label = "Method B"
    sample_id_col = None

detected_analytes = []
for analyte, (candidate, _) in DEFAULT_COLUMN_MAP.items():
    if layout.startswith("Long"):
        if candidate in columns:
            detected_analytes.append(analyte)
    else:
        if candidate in columns:
            detected_analytes.append(analyte)

default_selected = detected_analytes or list(NORMAL_RANGES)
selected_analytes = st.multiselect(
    "Analytes to analyze",
    options=list(NORMAL_RANGES),
    default=default_selected,
)
if not selected_analytes:
    st.warning("Select at least one analyte.")
    st.stop()

profile = st.radio(
    "Trueness-bias profile for P/F",
    [
        "MHS screenshot profile",
        "Initial EP09 profile",
        "Custom values below",
        "Upload a trueness-bias table",
    ],
    horizontal=True,
)
bias_profile = (
    INITIAL_EP09_BIAS.copy()
    if profile == "Initial EP09 profile"
    else MHS_SCREENSHOT_BIAS.copy()
)

bias_upload = None
if profile == "Upload a trueness-bias table":
    bias_upload = st.file_uploader(
        "Upload CSV/XLSX with columns such as Analyte and Bias %",
        type=["csv", "xlsx", "xls"],
        key="bias_upload",
    )
    if bias_upload is not None:
        try:
            if bias_upload.name.lower().endswith(".csv"):
                bias_table = pd.read_csv(bias_upload)
            else:
                bias_table = pd.read_excel(bias_upload)
            lower_names = {str(c).lower(): c for c in bias_table.columns}
            analyte_col = next(
                (lower_names[k] for k in lower_names if "analyte" in k or "measurand" in k),
                None,
            )
            bias_col = next(
                (lower_names[k] for k in lower_names if "bias" in k),
                None,
            )
            if analyte_col and bias_col:
                for _, row in bias_table.iterrows():
                    name = str(row[analyte_col]).strip().upper()
                    raw_value = str(row[bias_col]).replace("%", "").strip()
                    value = pd.to_numeric(pd.Series([raw_value]), errors="coerce").iloc[0]
                    if name in bias_profile and np.isfinite(value):
                        bias_profile[name] = float(value)
                st.success("Uploaded bias table applied.")
            else:
                st.error("Could not identify analyte and bias columns.")
        except Exception as exc:
            st.error(f"Could not read the bias table: {exc}")

config_rows = []
for analyte in selected_analytes:
    candidate, reference = DEFAULT_COLUMN_MAP[analyte]
    low, high = NORMAL_RANGES[analyte]
    row = {
        "Analyte": analyte,
        "Normal low": low,
        "Normal high": high,
        "Allowable bias %": bias_profile[analyte],
    }
    if layout.startswith("Long"):
        row["Value column"] = candidate if candidate in columns else columns[0]
        row["Reference column"] = reference if reference in columns else "None"
    else:
        row["Method A column"] = candidate if candidate in columns else columns[0]
        row["Method B column"] = candidate if candidate in columns else columns[0]
        row["Reference A column"] = reference if reference in columns else "None"
        row["Reference B column"] = reference if reference in columns else "None"
    config_rows.append(row)

config_df = pd.DataFrame(config_rows)

if profile == "Custom values below":
    st.info("Edit the Allowable bias % column directly.")

column_config = {
    "Analyte": st.column_config.TextColumn(disabled=True),
    "Normal low": st.column_config.NumberColumn(format="%.4g"),
    "Normal high": st.column_config.NumberColumn(format="%.4g"),
    "Allowable bias %": st.column_config.NumberColumn(format="%.4g"),
}
if layout.startswith("Long"):
    column_config["Value column"] = st.column_config.SelectboxColumn(
        options=columns,
        required=True,
    )
    column_config["Reference column"] = st.column_config.SelectboxColumn(
        options=["None"] + columns,
    )
else:
    for name in ["Method A column", "Method B column"]:
        column_config[name] = st.column_config.SelectboxColumn(
            options=columns,
            required=True,
        )
    for name in ["Reference A column", "Reference B column"]:
        column_config[name] = st.column_config.SelectboxColumn(
            options=["None"] + columns,
        )

edited_config = st.data_editor(
    config_df,
    use_container_width=True,
    hide_index=True,
    column_config=column_config,
    key="analyte_config",
)

st.subheader("4. Filters, outliers, corrections, and optional analyses")

with st.container(border=True):
    st.markdown("**Filtering**")
    filter_to_normal = st.checkbox(
        "Filter to the configured normal range before paired analysis",
        value=True,
    )
    require_reference_normal = st.checkbox(
        "Also require paired Sysmex/reference values to be within the normal range",
        value=False,
        disabled=not filter_to_normal,
    )

with st.container(border=True):
    st.markdown("**Outlier handling**")
    outlier_enabled = st.checkbox(
        "Run an outlier-removal sensitivity analysis",
        value=False,
    )
    if layout.startswith("Long"):
        outlier_stage = st.selectbox(
            "Outlier stage",
            [
                "Donor-pair BA bias",
                "Pre-averaging replicate residual",
            ],
            disabled=not outlier_enabled,
        )
    else:
        outlier_stage = "Donor-pair BA bias"

    outlier_method = st.selectbox(
        "Outlier method",
        ["Auto recommendation", "Grubbs", "MAD", "Generalized ESD"],
        disabled=not outlier_enabled,
    )
    outlier_alpha = st.number_input(
        "Outlier alpha",
        min_value=0.001,
        max_value=0.20,
        value=0.05,
        step=0.01,
        disabled=not outlier_enabled,
    )
    mad_threshold = st.number_input(
        "MAD threshold",
        min_value=2.0,
        max_value=8.0,
        value=3.5,
        step=0.1,
        disabled=not outlier_enabled,
    )
    max_outliers = st.number_input(
        "Maximum outliers removed per analyte/arm",
        min_value=1,
        max_value=5,
        value=1,
        step=1,
        disabled=not outlier_enabled,
    )
    st.caption(
        "Auto uses Grubbs when the relevant distribution is compatible with "
        "normality and n is adequate; otherwise it uses robust MAD. ESD is best "
        "reserved for a pre-specified expectation of multiple outliers."
    )

with st.container(border=True):
    st.markdown("**Optional correction/sensitivity methods**")
    include_direct_log = st.checkbox("Direct log-ratio sensitivity", value=True)
    include_median_log = st.checkbox(
        "Median log-ratio recalibration",
        value=True,
    )
    median_log_loo = st.checkbox(
        "Use leave-one-out median factor",
        value=True,
        disabled=not include_median_log,
    )
    include_sysmex_absolute = st.checkbox(
        "Donor-specific Sysmex absolute difference-in-differences",
        value=False,
    )
    include_sysmex_log = st.checkbox(
        "Donor-specific Sysmex log-ratio correction",
        value=False,
    )
    st.caption(
        "Sysmex corrections require paired reference columns. Absolute correction "
        "subtracts the donor-specific reference shift; the log correction subtracts "
        "the donor-specific log reference ratio."
    )

with st.container(border=True):
    st.markdown("**Tests, power, and SD distance**")
    power_enabled = st.checkbox(
        "Compute optional sample-size estimates",
        value=False,
    )
    power = st.select_slider(
        "Power target",
        options=[0.80, 0.85, 0.90],
        value=0.80,
        disabled=not power_enabled,
    )
    sd_distance_enabled = st.checkbox(
        "Compute optional SD-distance coverage table",
        value=False,
    )
    st.caption(
        "The primary output is BA mean percent bias with its 95% CI and P/F. "
        "Paired t/Wilcoxon results are secondary. Mean-bias equivalence is the "
        "default power concept; LoA precision is also reported."
    )

run_button = st.button("Run interchangeability analysis", type="primary", use_container_width=True)
if not run_button:
    st.stop()

all_summary = []
all_diagnostics = []
all_pairs = []
all_outliers = []
all_sd = []
plot_bytes = {}
run_messages = []

progress = st.progress(0.0)
status = st.empty()

for position, config in edited_config.reset_index(drop=True).iterrows():
    analyte = str(config["Analyte"]).strip().upper()
    status.write(f"Analyzing {analyte}…")

    normal_low = float(config["Normal low"])
    normal_high = float(config["Normal high"])
    allowable_bias = float(config["Allowable bias %"])

    try:
        if layout.startswith("Long"):
            value_col = str(config["Value column"])
            reference_col = str(config["Reference column"])
            reference_col = None if reference_col == "None" else reference_col

            pairs, pre_outliers, audit = prepare_long_pairs(
                data,
                donor_col=donor_col,
                type_col=type_col,
                method_a_label=method_a_label,
                method_b_label=method_b_label,
                value_col=value_col,
                reference_col=reference_col,
                normal_low=normal_low,
                normal_high=normal_high,
                filter_to_normal=filter_to_normal,
                require_reference_normal=require_reference_normal,
                global_flag_col=global_flag_col,
                sample_id_col=sample_id_col,
                outlier_enabled=outlier_enabled,
                outlier_stage=outlier_stage,
                outlier_method=outlier_method,
                outlier_alpha=float(outlier_alpha),
                mad_threshold=float(mad_threshold),
                max_outliers=int(max_outliers),
            )
        else:
            method_a_col = str(config["Method A column"])
            method_b_col = str(config["Method B column"])
            ref_a_col = str(config["Reference A column"])
            ref_b_col = str(config["Reference B column"])
            ref_a_col = None if ref_a_col == "None" else ref_a_col
            ref_b_col = None if ref_b_col == "None" else ref_b_col

            pairs, pre_outliers, audit = prepare_wide_pairs(
                data,
                donor_col=donor_col,
                method_a_col=method_a_col,
                method_b_col=method_b_col,
                ref_a_col=ref_a_col,
                ref_b_col=ref_b_col,
                normal_low=normal_low,
                normal_high=normal_high,
                filter_to_normal=filter_to_normal,
                require_reference_normal=require_reference_normal,
                global_flag_col=global_flag_col,
            )

        if len(pairs) < 2:
            run_messages.append(
                f"{analyte}: skipped because fewer than two paired observations remained."
            )
            continue

        result = analyze_pairs(
            pairs,
            analyte=analyte,
            allowable_bias=allowable_bias,
            normal_low=normal_low,
            normal_high=normal_high,
            outlier_enabled=outlier_enabled,
            outlier_stage=outlier_stage,
            outlier_method=outlier_method,
            outlier_alpha=float(outlier_alpha),
            mad_threshold=float(mad_threshold),
            max_outliers=int(max_outliers),
            include_direct_log=include_direct_log,
            include_median_log=include_median_log,
            median_log_loo=median_log_loo,
            include_sysmex_absolute=include_sysmex_absolute,
            include_sysmex_log=include_sysmex_log,
            power_enabled=power_enabled,
            power=float(power),
            sd_distance_enabled=sd_distance_enabled,
        )

        result["summary"]["Rows after filters"] = audit["rows_after_filters"]
        result["summary"]["Paired before outlier"] = audit["paired_donors_before_pair_outlier"]
        result["pairs"].insert(0, "Analyte", analyte)
        result["diagnostics"]["Rows initial"] = audit["rows_initial"]
        result["diagnostics"]["Rows after filters"] = audit["rows_after_filters"]

        if pre_outliers:
            pre_frame = pd.DataFrame(pre_outliers)
            pre_frame.insert(0, "Analyte", analyte)
            result_outliers = pd.concat(
                [pre_frame, result["outliers"]],
                ignore_index=True,
                sort=False,
            )
        else:
            result_outliers = result["outliers"].copy()

        if not result_outliers.empty and "Analyte" not in result_outliers.columns:
            result_outliers.insert(0, "Analyte", analyte)

        all_summary.append(result["summary"])
        all_diagnostics.append(result["diagnostics"])
        all_pairs.append(result["pairs"])
        if not result_outliers.empty:
            all_outliers.append(result_outliers)
        if not result["sd_distance"].empty:
            all_sd.append(result["sd_distance"])

        outlier_donors = (
            set(result_outliers["donor"].dropna().astype(str))
            if not result_outliers.empty and "donor" in result_outliers.columns
            else set()
        )
        for method in result["summary"]["Method"]:
            fig = plot_ba_percent(
                result["pairs"],
                analyte=analyte,
                method=method,
                allowable_bias=allowable_bias,
                outlier_donors=outlier_donors,
            )
            plot_bytes[f"{analyte}_{method.replace(' ', '_')}_BA.png"] = figure_to_png_bytes(fig)
            plt.close(fig)

    except Exception as exc:
        run_messages.append(f"{analyte}: {type(exc).__name__}: {exc}")

    progress.progress((position + 1) / len(edited_config))

status.empty()
progress.empty()

if not all_summary:
    st.error("No analyte produced a valid paired analysis.")
    if run_messages:
        st.code("\n".join(run_messages))
    st.stop()

summary_df = pd.concat(all_summary, ignore_index=True)
diagnostics_df = pd.concat(all_diagnostics, ignore_index=True)
pairs_df = pd.concat(all_pairs, ignore_index=True)
outliers_df = (
    pd.concat(all_outliers, ignore_index=True, sort=False)
    if all_outliers
    else pd.DataFrame(
        columns=["Analyte", "stage", "donor", "sample_id", "method"]
    )
)
sd_df = pd.concat(all_sd, ignore_index=True) if all_sd else pd.DataFrame()

# Holm correction across analytes for each primary paired test.
diagnostics_df["Recommended test Holm p"] = holm_adjust(
    diagnostics_df["recommended_test_p"].to_numpy()
)
diagnostics_df["Paired t Holm p"] = holm_adjust(
    diagnostics_df["paired_t_p"].to_numpy()
)
diagnostics_df["Wilcoxon Holm p"] = holm_adjust(
    diagnostics_df["wilcoxon_p"].to_numpy()
)

raw_summary = summary_df[summary_df["Method"] == "Raw percent BA"].copy()
pass_count = int((raw_summary["P/F"] == "PASS").sum())
fail_count = int((raw_summary["P/F"] == "FAIL").sum())

st.subheader("Results")
metric1, metric2, metric3, metric4 = st.columns(4)
metric1.metric("Analytes analyzed", raw_summary["Analyte"].nunique())
metric2.metric("Raw BA PASS", pass_count)
metric3.metric("Raw BA FAIL", fail_count)
metric4.metric("Paired observations", int(raw_summary["n paired"].sum()))

tab1, tab2, tab3, tab4, tab5 = st.tabs(
    [
        "BA and P/F",
        "Normality and tests",
        "Outliers",
        "Plots",
        "Power and SD distance",
    ]
)

with tab1:
    st.dataframe(
        summary_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "P/F": st.column_config.TextColumn(),
            "Mean BA bias, %": st.column_config.NumberColumn(format="%.3f"),
            "BA 95% CI low, %": st.column_config.NumberColumn(format="%.3f"),
            "BA 95% CI high, %": st.column_config.NumberColumn(format="%.3f"),
        },
    )
    st.caption(
        "Primary P/F: PASS only when the full BA mean-bias 95% CI lies inside "
        "the selected ± allowable trueness-bias interval."
    )

with tab2:
    st.dataframe(diagnostics_df, use_container_width=True, hide_index=True)
    for _, row in diagnostics_df.iterrows():
        st.info(f"**{row['Analyte']}** — {row['Recommendation']}")

with tab3:
    if outliers_df.empty:
        st.success("No outliers were removed or flagged under the selected settings.")
    else:
        st.dataframe(outliers_df, use_container_width=True, hide_index=True)
        st.warning(
            "Outlier removal is a sensitivity analysis unless the exclusion was "
            "pre-specified or independently supported by QC/preanalytical evidence."
        )

with tab4:
    analyte_choice = st.selectbox(
        "Analyte",
        sorted(summary_df["Analyte"].unique()),
        key="plot_analyte",
    )
    methods_for_analyte = summary_df.loc[
        summary_df["Analyte"] == analyte_choice,
        "Method",
    ].tolist()
    method_choice = st.selectbox(
        "Method",
        methods_for_analyte,
        key="plot_method",
    )
    plot_key = f"{analyte_choice}_{method_choice.replace(' ', '_')}_BA.png"
    st.image(plot_bytes[plot_key], use_container_width=True)

    forest_fig = plot_forest(summary_df, method_choice)
    forest_bytes = figure_to_png_bytes(forest_fig)
    plt.close(forest_fig)
    st.image(forest_bytes, use_container_width=True)
    plot_bytes[f"ALL_ANALYTES_{method_choice.replace(' ', '_')}_forest.png"] = forest_bytes

with tab5:
    power_columns = [
        "Analyte",
        "Method",
        "Mean-bias equivalence n",
        "Mean-bias CI precision n",
        "LoA precision n",
    ]
    if power_enabled:
        st.dataframe(
            summary_df[power_columns],
            use_container_width=True,
            hide_index=True,
        )
        st.caption(
            "Mean-bias equivalence n is the default power estimate. Infinite means "
            "the observed bias is already at or outside the selected allowable margin."
        )
    else:
        st.info("Power analysis was not selected.")

    if sd_distance_enabled and not sd_df.empty:
        st.dataframe(sd_df, use_container_width=True, hide_index=True)
    elif sd_distance_enabled:
        st.warning("SD-distance analysis could not be estimated.")
    else:
        st.info("SD-distance analysis was not selected.")

if run_messages:
    with st.expander("Warnings and skipped items"):
        st.code("\n".join(run_messages))

settings = {
    "uploaded_file": uploaded.name,
    "sheet": sheet_name,
    "layout": layout,
    "filter_to_normal": filter_to_normal,
    "require_reference_normal": require_reference_normal,
    "outlier_enabled": outlier_enabled,
    "outlier_stage": outlier_stage,
    "outlier_method": outlier_method,
    "outlier_alpha": float(outlier_alpha),
    "mad_threshold": float(mad_threshold),
    "max_outliers": int(max_outliers),
    "direct_log": include_direct_log,
    "median_log": include_median_log,
    "median_log_loo": median_log_loo,
    "sysmex_absolute": include_sysmex_absolute,
    "sysmex_log": include_sysmex_log,
    "power_enabled": power_enabled,
    "power": float(power),
    "sd_distance": sd_distance_enabled,
    "trueness_profile": profile,
}
settings_text = json.dumps(settings, indent=2)

excel_bytes = dataframe_to_excel_bytes(
    {
        "Summary": summary_df,
        "Diagnostics": diagnostics_df,
        "Paired points": pairs_df,
        "Outliers": outliers_df,
        "SD distance": sd_df,
        "Configuration": edited_config,
    }
)

results_zip = build_results_zip(
    summary=summary_df,
    diagnostics=diagnostics_df,
    pairs=pairs_df,
    outliers=outliers_df,
    sd_distance=sd_df,
    plots=plot_bytes,
    settings_text=settings_text,
)

download_col1, download_col2 = st.columns(2)
download_col1.download_button(
    "Download Excel report",
    data=excel_bytes,
    file_name="Interchangeability_App_results.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    use_container_width=True,
)
download_col2.download_button(
    "Download complete results ZIP",
    data=results_zip,
    file_name="Interchangeability_App_results.zip",
    mime="application/zip",
    use_container_width=True,
)
