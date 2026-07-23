from __future__ import annotations

import io
import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from analysis_core import (
    DEFAULT_ANALYTE_CONFIG,
    METHOD_DESCRIPTIONS,
    build_runtime_analytes,
    canonicalize_input_dataframe,
    collect_output_files,
    default_analyte_table,
    normalize_bool,
    parse_blood_sample_id,
    read_binary,
    run_complete_notebook_pipeline,
)


st.set_page_config(
    page_title="PROXIMA Trueness + Bland–Altman",
    page_icon="🧪",
    layout="wide",
)

st.title("🧪 PROXIMA Trueness + Bland–Altman App")
st.caption(
    "Streamlit implementation of the final validated notebook: replicate-level "
    "APT generalized-ESD trueness regression plus paired-specimen M02/M11/M05 "
    "Bland–Altman analysis, native-unit context, acceptance screens, full audits, "
    "Excel export, and all plot outputs."
)

with st.expander("Method lock and interpretation", expanded=False):
    st.markdown(
        """
- **Trueness branch:** `global_flag == FALSE` → validated generalized ESD on raw linked replicate rows → analyte-specific removal → donor means → Huber regression → Pearson *r* with Fisher-z 95% CI.
- **Bland–Altman branch:** `global_flag == FALSE` and analyte flags → residual normality by Shapiro–Wilk → manual Grubbs `Gcrit` if normal or MAD modified-Z if non-normal → maximum one linked replicate removed per analyte → donor means.
- **M02:** direct signed MHS specimen-A versus specimen-B comparison.
- **M11:** signed Sysmex reference-adjusted comparison and the preferred reference-adjusted option.
- **M05:** absolute-gap magnitude sensitivity analysis; it is not signed and does not replace M02/M11.

Use de-identified study data only. AC/CLIA entries are editable contextual criteria and must be verified before final or regulatory use.
"""
    )


@st.cache_data(show_spinner=False)
def read_upload(file_name: str, payload: bytes) -> dict[str, pd.DataFrame]:
    suffix = Path(file_name).suffix.lower()
    if suffix == ".csv":
        return {"CSV": pd.read_csv(io.BytesIO(payload))}
    excel = pd.ExcelFile(io.BytesIO(payload))
    return {
        sheet: pd.read_excel(io.BytesIO(payload), sheet_name=sheet)
        for sheet in excel.sheet_names
    }


def option_index(options: list[str], preferred: str | None, fallback: int = 0) -> int:
    if preferred is not None and preferred in options:
        return options.index(preferred)
    return min(fallback, max(0, len(options) - 1))


def metadata_preview(
    data: pd.DataFrame,
    *,
    sample_id_col: str | None,
    donor_col: str | None,
    specimen_col: str | None,
    replicate_col: str | None,
) -> pd.DataFrame:
    if sample_id_col is not None:
        parsed = pd.DataFrame([parse_blood_sample_id(v) for v in data[sample_id_col]])
        parsed.insert(0, "source_row", np.arange(1, len(parsed) + 1))
        return parsed
    preview = pd.DataFrame({
        "source_row": np.arange(1, len(data) + 1),
        "donor": data[donor_col].astype(str).str.strip(),
        "specimen_type": data[specimen_col].astype(str).str.strip(),
        "replicate_number": pd.to_numeric(data[replicate_col], errors="coerce"),
    })
    preview["parse_ok"] = preview[["donor", "specimen_type"]].notna().all(axis=1)
    return preview


def default_group_rows(specimen_types: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if any(value in specimen_types for value in ["Mc", "Mc2"]):
        selected = [value for value in ["Mc", "Mc2"] if value in specimen_types]
        rows.append({
            "Enabled": True,
            "Group name": "Capillary_Mc_Mc2",
            "Specimen types (comma-separated)": ", ".join(selected),
        })
    known = [
        ("MK2", "Venous_MK2", True),
        ("MK3", "Venous_MK3", False),
        ("Dir", "Direct_Capillary_Dir", False),
        ("Tas", "TASSO_Tas", False),
    ]
    already = {token.strip() for row in rows for token in str(row["Specimen types (comma-separated)"]).split(",")}
    for specimen, name, default_enabled in known:
        if specimen in specimen_types and specimen not in already:
            rows.append({
                "Enabled": default_enabled,
                "Group name": name,
                "Specimen types (comma-separated)": specimen,
            })
            already.add(specimen)
    for specimen in specimen_types:
        if specimen not in already:
            rows.append({
                "Enabled": False,
                "Group name": f"Type_{specimen}",
                "Specimen types (comma-separated)": specimen,
            })
    return pd.DataFrame(rows or [{
        "Enabled": True,
        "Group name": "Analysis_Group_1",
        "Specimen types (comma-separated)": "",
    }])


def parse_group_editor(editor: pd.DataFrame, available_types: set[str]) -> tuple[dict[str, list[str]], list[str]]:
    groups: dict[str, list[str]] = {}
    enabled: list[str] = []
    for _, row in editor.iterrows():
        if not bool(row.get("Enabled", False)):
            continue
        name = str(row.get("Group name", "")).strip()
        types = [token.strip() for token in str(row.get("Specimen types (comma-separated)", "")).split(",") if token.strip()]
        if not name:
            raise ValueError("Every enabled regression group needs a name.")
        if name in groups:
            raise ValueError(f"Duplicate regression group name: {name}")
        if not types:
            raise ValueError(f"{name}: select at least one specimen type.")
        missing = sorted(set(types) - available_types)
        if missing:
            raise ValueError(f"{name}: specimen types not present in the data: {', '.join(missing)}")
        groups[name] = types
        enabled.append(name)
    if not enabled:
        raise ValueError("Enable at least one trueness-regression group.")
    return groups, enabled


def pf_styler(frame: pd.DataFrame):
    def highlight(value):
        text = str(value).lower()
        if text == "pass":
            return "background-color: #d9ead3"
        if text == "fail":
            return "background-color: #f4cccc"
        if "unverified" in text:
            return "background-color: #fff2cc"
        return ""
    return frame.style.map(highlight)


uploaded = st.file_uploader("1. Upload CSV or Excel data", type=["csv", "xlsx", "xls"])
if uploaded is None:
    st.info("Upload the combined-results CSV/XLSX to configure and run the analysis.")
    st.stop()

sheets = read_upload(uploaded.name, uploaded.getvalue())
sheet_name = st.selectbox("Worksheet", list(sheets), index=0)
data = sheets[sheet_name].copy()
columns = list(map(str, data.columns))
data.columns = columns

st.success(f"Loaded {len(data):,} rows × {len(columns):,} columns from **{sheet_name}**.")
with st.expander("Uploaded-data preview", expanded=False):
    st.dataframe(data.head(50), use_container_width=True)

st.header("2. Metadata and specimen parsing")
meta1, meta2, meta3 = st.columns(3)
batch_options = ["<create from row number>"] + columns
batch_choice = meta1.selectbox(
    "Batch/replicate row ID",
    batch_options,
    index=option_index(batch_options, "batch_id"),
)
global_options = ["<no column>"] + columns
global_choice = meta2.selectbox(
    "Global flag column",
    global_options,
    index=option_index(global_options, "global_flag"),
)
treat_no_global_as_false = meta3.checkbox(
    "Treat all rows as global_flag = FALSE when no column is selected",
    value=False,
    disabled=global_choice != "<no column>",
)

id_mode = st.radio(
    "How should donor/specimen/replicate identity be obtained?",
    ["Parse a sample-ID column", "Use explicit donor/specimen/replicate columns"],
    horizontal=True,
)

sample_id_col: str | None = None
donor_col: str | None = None
specimen_col: str | None = None
replicate_col: str | None = None
if id_mode.startswith("Parse"):
    sample_id_col = st.selectbox(
        "Sample-ID column (expected ending such as D05-Mc-2)",
        columns,
        index=option_index(columns, "bloodSampleId"),
    )
else:
    id1, id2, id3 = st.columns(3)
    donor_col = id1.selectbox("Donor column", columns, index=option_index(columns, "Donor"))
    specimen_col = id2.selectbox("Specimen-type column", columns, index=option_index(columns, "Type"))
    replicate_col = id3.selectbox("Replicate-number column", columns, index=option_index(columns, "Replicate"))

try:
    parsed_preview = metadata_preview(
        data,
        sample_id_col=sample_id_col,
        donor_col=donor_col,
        specimen_col=specimen_col,
        replicate_col=replicate_col,
    )
    specimen_types = sorted(parsed_preview.loc[parsed_preview["parse_ok"], "specimen_type"].dropna().astype(str).unique().tolist())
    donor_count = int(parsed_preview.loc[parsed_preview["parse_ok"], "donor"].nunique())
    parse_failures = int((~parsed_preview["parse_ok"].fillna(False)).sum())
except Exception as exc:
    st.error(f"Could not inspect specimen IDs: {exc}")
    st.stop()

m1, m2, m3 = st.columns(3)
m1.metric("Parsed donors", donor_count)
m2.metric("Detected specimen types", len(specimen_types))
m3.metric("Unparsed rows", parse_failures)
st.write("Specimen labels detected: " + (", ".join(specimen_types) if specimen_types else "none"))
if parse_failures:
    st.warning("Unparsed rows are excluded by both notebook pipelines and are listed in the output audit.")

if global_choice != "<no column>":
    global_bool = data[global_choice].map(normalize_bool)
    false_n = int(global_bool.eq(False).sum())
    true_n = int(global_bool.eq(True).sum())
    unknown_n = int(global_bool.isna().sum())
    st.caption(f"Global-flag preview: FALSE={false_n}, TRUE={true_n}, missing/unparseable={unknown_n}.")

st.header("3. Analytes, mappings, ranges, and acceptance criteria")
defaults = default_analyte_table()
detected = []
for _, row in defaults.iterrows():
    if row["MHS column"] in columns and row["Reference column"] in columns:
        detected.append(str(row["Analyte"]))
selected_analytes = st.multiselect(
    "Analytes to run",
    options=list(DEFAULT_ANALYTE_CONFIG),
    default=detected or list(DEFAULT_ANALYTE_CONFIG),
)
if not selected_analytes:
    st.warning("Select at least one analyte.")
    st.stop()

config_df = defaults[defaults["Analyte"].isin(selected_analytes)].copy()
for index, row in config_df.iterrows():
    flag_name = str(row["Flag column"])
    if flag_name not in columns:
        config_df.at[index, "Flag column"] = "<none>"

edited_analytes = st.data_editor(
    config_df,
    hide_index=True,
    use_container_width=True,
    num_rows="fixed",
    column_config={
        "Analyte": st.column_config.TextColumn(disabled=True),
        "MHS column": st.column_config.SelectboxColumn(options=columns, required=True),
        "Reference column": st.column_config.SelectboxColumn(options=columns, required=True),
        "Flag column": st.column_config.SelectboxColumn(options=["<none>"] + columns),
        "Normal low": st.column_config.NumberColumn(format="%.6g", required=True),
        "Normal high": st.column_config.NumberColumn(format="%.6g", required=True),
        "Unit": st.column_config.TextColumn(required=True),
        "AC limit %": st.column_config.NumberColumn(format="%.4g", required=True),
        "CLIA limit %": st.column_config.NumberColumn(format="%.4g"),
    },
    key="analyte_editor_final",
)
st.caption(
    "AC/CLIA values are editable. A blank CLIA cell means that no CLIA percentage screen is applied for that analyte. "
    "Normal ranges also control the red vertical lines in regression plots and the 50%-midpoint context."
)

st.header("4. Trueness-regression configuration")
initial_groups = default_group_rows(specimen_types)
group_editor = st.data_editor(
    initial_groups,
    hide_index=True,
    use_container_width=True,
    num_rows="dynamic",
    column_config={
        "Enabled": st.column_config.CheckboxColumn(),
        "Group name": st.column_config.TextColumn(required=True),
        "Specimen types (comma-separated)": st.column_config.TextColumn(required=True),
    },
    key="group_editor_final",
)

reg1, reg2 = st.columns(2)
run_huber_bootstrap = reg1.checkbox(
    "Supplemental donor-level Huber bootstrap CIs",
    value=False,
    help="Disabled by default to match the notebook verification run. Huber point estimates and Fisher-z correlation CIs are always computed.",
)
n_huber_bootstraps = reg2.number_input(
    "Huber bootstrap resamples",
    min_value=200,
    max_value=20_000,
    value=1000,
    step=200,
    disabled=not run_huber_bootstrap,
)

with st.expander("Regression-plot label placement", expanded=False):
    st.caption(
        "Automatic placement searches for blank space and avoids observations, the fitted line, "
        "and the red normal-range boundaries. For reproducible exported PNGs, manual positions can "
        "also be selected independently for the equation and the R/95% CI block."
    )
    placement_options = {
        "Automatic — avoid points and red lines": "auto",
        "Top left": "top_left",
        "Top centre": "top_center",
        "Top right": "top_right",
        "Middle left": "middle_left",
        "Middle right": "middle_right",
        "Bottom left": "bottom_left",
        "Bottom centre": "bottom_center",
        "Bottom right": "bottom_right",
    }
    place1, place2 = st.columns(2)
    equation_position_label = place1.selectbox(
        "Huber equation position",
        list(placement_options),
        index=0,
        key="regression_equation_position",
    )
    metrics_position_label = place2.selectbox(
        "R, Fisher-z 95% CI, n and range position",
        list(placement_options),
        index=0,
        key="regression_metrics_position",
    )
    regression_equation_position = placement_options[equation_position_label]
    regression_metrics_position = placement_options[metrics_position_label]
    st.info(
        "The exported figures are static PNGs, so browser click-and-drag movement would not persist "
        "in the downloaded files. These controls provide the same adjustment reproducibly; automatic "
        "placement is recommended."
    )

try:
    preview_groups, preview_enabled_groups = parse_group_editor(group_editor, set(specimen_types))
except Exception:
    preview_groups, preview_enabled_groups = {}, []

cross_enabled = st.checkbox("Also compute a matched-donor cross-specimen trueness regression", value=len(preview_enabled_groups) >= 2)
cross_pair: tuple[str, str] | None = None
if cross_enabled and len(preview_enabled_groups) >= 2:
    cross1, cross2 = st.columns(2)
    group_a = cross1.selectbox("Cross-specimen group A", preview_enabled_groups, index=0)
    group_b_options = [name for name in preview_enabled_groups if name != group_a]
    group_b = cross2.selectbox("Cross-specimen group B", group_b_options, index=0)
    cross_pair = (group_a, group_b)
elif cross_enabled:
    st.warning("Enable at least two regression groups to run the cross-specimen comparison.")

st.header("5. Paired-specimen Bland–Altman configuration")
if len(specimen_types) < 2:
    st.error("At least two specimen labels are required for paired-specimen Bland–Altman analysis.")
    st.stop()

default_a = [value for value in ["Mc", "Mc2"] if value in specimen_types] or [specimen_types[0]]
default_b = ["MK2"] if "MK2" in specimen_types else [value for value in specimen_types if value not in default_a][:1]
ba1, ba2 = st.columns(2)
ba_a_types = ba1.multiselect("Specimen A types", specimen_types, default=default_a)
ba_b_types = ba2.multiselect("Specimen B types", specimen_types, default=default_b)
label1, label2 = st.columns(2)
ba_a_label = label1.text_input("Specimen A display label", value="Capillary")
ba_b_label = label2.text_input("Specimen B display label", value="MK2")
criteria_confirmed = st.checkbox(
    "I reviewed and verified the edited AC/CLIA table for this analysis",
    value=False,
)
if not criteria_confirmed:
    st.warning("Pass/Fail labels will be generated as **UNVERIFIED / PROVISIONAL CONTEXT ONLY**.")

plot_mode = st.radio(
    "Plot output",
    ["Full: all individual plots and combined panels", "Panels only: faster"],
    horizontal=True,
)
plot_mode_value = "full" if plot_mode.startswith("Full") else "panels_only"

with st.expander("Advanced notebook constants", expanded=False):
    adv1, adv2, adv3 = st.columns(3)
    shapiro_alpha = adv1.number_input("Shapiro branch alpha", min_value=0.001, max_value=0.20, value=0.05, step=0.01)
    outlier_alpha = adv2.number_input("Manual Grubbs alpha", min_value=0.001, max_value=0.20, value=0.05, step=0.01)
    mad_z = adv3.number_input("MAD modified-Z threshold", min_value=2.0, max_value=8.0, value=3.5, step=0.1)
    adv4, adv5, adv6 = st.columns(3)
    ba_bootstraps = adv4.number_input("BA donor bootstrap iterations", min_value=1000, max_value=100_000, value=10_000, step=1000)
    donor_profile_bootstraps = adv5.number_input("Within-donor profile bootstraps", min_value=500, max_value=50_000, value=5000, step=500)
    random_seed = adv6.number_input("Random seed", min_value=1, max_value=2_147_483_647, value=20260723, step=1)
    st.caption("The defaults exactly match the final notebook.")

run_clicked = st.button("Run final PROXIMA analysis", type="primary", use_container_width=True)
if run_clicked:
    try:
        runtime_analytes, ac_limits, clia_limits = build_runtime_analytes(edited_analytes.to_dict("records"))
        regression_groups, groups_to_run = parse_group_editor(group_editor, set(specimen_types))
        if not ba_a_types or not ba_b_types:
            raise ValueError("Select at least one specimen type in both Bland–Altman arms.")
        overlap = sorted(set(ba_a_types) & set(ba_b_types))
        if overlap:
            raise ValueError("The Bland–Altman specimen arms must not overlap: " + ", ".join(overlap))
        canonical = canonicalize_input_dataframe(
            data,
            runtime_analytes,
            batch_id_col=None if batch_choice == "<create from row number>" else batch_choice,
            global_flag_col=None if global_choice == "<no column>" else global_choice,
            sample_id_col=sample_id_col,
            donor_col=donor_col,
            specimen_col=specimen_col,
            replicate_col=replicate_col,
            treat_missing_global_flag_as_false=treat_no_global_as_false,
        )

        with st.status("Running the notebook pipeline from beginning to end…", expanded=True) as status:
            st.write("Applying global-flag and ID parsing rules…")
            st.write("Running replicate-level generalized-ESD trueness regressions…")
            st.write("Running M02/M11/M05 Bland–Altman analyses and bootstraps…")
            st.write("Saving all plots, CSV audits, Excel sheets, and ZIP outputs…")
            with tempfile.TemporaryDirectory(prefix="proxima_app_") as temp_dir:
                result = run_complete_notebook_pipeline(
                    canonical,
                    temp_dir,
                    runtime_analytes=runtime_analytes,
                    ac_limits=ac_limits,
                    clia_limits=clia_limits,
                    regression_groups=regression_groups,
                    groups_to_run=groups_to_run,
                    ba_specimen_a_types=ba_a_types,
                    ba_specimen_b_types=ba_b_types,
                    ba_specimen_a_label=ba_a_label.strip() or "Specimen A",
                    ba_specimen_b_label=ba_b_label.strip() or "Specimen B",
                    criteria_confirmed=criteria_confirmed,
                    run_huber_bootstrap=run_huber_bootstrap,
                    n_huber_bootstraps=int(n_huber_bootstraps),
                    plot_mode=plot_mode_value,
                    shapiro_alpha=float(shapiro_alpha),
                    outlier_alpha=float(outlier_alpha),
                    mad_z_threshold=float(mad_z),
                    ba_bootstrap_iterations=int(ba_bootstraps),
                    donor_profile_bootstraps=int(donor_profile_bootstraps),
                    random_seed=int(random_seed),
                    cross_specimen_groups=cross_pair,
                    regression_equation_position=regression_equation_position,
                    regression_metrics_position=regression_metrics_position,
                )

                regression_metrics = pd.concat(
                    [bundle["metrics"] for bundle in result["regression"]["groups"].values()],
                    ignore_index=True,
                ) if result["regression"]["groups"] else pd.DataFrame()
                regression_outliers = pd.concat(
                    [bundle["outliers"] for bundle in result["regression"]["groups"].values()],
                    ignore_index=True,
                ) if result["regression"]["groups"] else pd.DataFrame()
                regression_diagnostics = pd.concat(
                    [bundle["diagnostics"] for bundle in result["regression"]["groups"].values()],
                    ignore_index=True,
                ) if result["regression"]["groups"] else pd.DataFrame()
                plot_files = collect_output_files(result["output_dir"], suffixes=[".png"])
                plot_bytes = {name: read_binary(path) for name, path in plot_files.items()}
                inventory = list(collect_output_files(result["output_dir"]).keys())

                st.session_state["proxima_final_results"] = {
                    "zip": read_binary(result["zip_path"]),
                    "xlsx": read_binary(result["combined_xlsx"]),
                    "summary_csv": read_binary(result["summary_csv"]),
                    "canonical_csv": canonical.to_csv(index=False).encode("utf-8"),
                    "regression_metrics": regression_metrics,
                    "regression_outliers": regression_outliers,
                    "regression_diagnostics": regression_diagnostics,
                    "regression_groups": list(result["regression"]["groups"]),
                    "cross_metrics": result["regression"].get("cross_specimen_metrics", pd.DataFrame()).copy(),
                    "ba_context": result["ba"]["combined_context"].copy(),
                    "ba_percent": result["ba"]["percent_results"].copy(),
                    "ba_native": result["ba"]["native_results"].copy(),
                    "ba_donor_values": result["ba"]["donor_values"].copy(),
                    "ba_removed": result["ba"]["removed_outliers"].copy(),
                    "ba_audit": result["ba"]["outlier_audit"].copy(),
                    "ba_analyte_flags": result["ba"]["analyte_flag_exclusions"].copy(),
                    "global_exclusions": result["ba"]["global_exclusions"].copy(),
                    "criteria": result["ba"]["criteria"].copy(),
                    "plots": plot_bytes,
                    "inventory": inventory,
                    "input_rows": len(canonical),
                    "global_false_rows": int(canonical["global_flag"].map(normalize_bool).eq(False).sum()),
                    "selected_analytes": list(runtime_analytes),
                    "criteria_confirmed": criteria_confirmed,
                }
            status.update(label="Analysis complete", state="complete", expanded=False)
        st.success("The final notebook-equivalent app analysis completed successfully.")
    except Exception as exc:
        st.exception(exc)


results = st.session_state.get("proxima_final_results")
if results is None:
    st.stop()

st.divider()
st.header("Results")
summary1, summary2, summary3, summary4 = st.columns(4)
summary1.metric("Input rows", f"{results['input_rows']:,}")
summary2.metric("global_flag = FALSE", f"{results['global_false_rows']:,}")
summary3.metric("Analytes", len(results["selected_analytes"]))
summary4.metric("BA outliers removed", len(results["ba_removed"]))

if not results["criteria_confirmed"]:
    st.warning("Acceptance results are provisional because the criteria-verification checkbox was not confirmed.")

reg_tab, ba_tab, audit_tab, plot_tab, download_tab = st.tabs([
    "Trueness regression",
    "Bland–Altman",
    "Exclusions and outliers",
    "Plots",
    "Downloads",
])

with reg_tab:
    st.subheader("Huber regression and Fisher-z correlation")
    st.dataframe(results["regression_metrics"], use_container_width=True, hide_index=True)
    if not results["cross_metrics"].empty:
        st.subheader("Matched-donor cross-specimen regression")
        st.dataframe(results["cross_metrics"], use_container_width=True, hide_index=True)
    regression_plot_names = [name for name in results["plots"] if "regression" in name.lower()]
    if regression_plot_names:
        chosen = st.selectbox("Regression plot", regression_plot_names, key="reg_plot_selector")
        st.image(results["plots"][chosen], caption=chosen, use_container_width=True)

with ba_tab:
    st.info(
        "M11 is the preferred signed reference-adjusted option; M02 is the direct signed comparison; "
        "M05 is magnitude-only sensitivity. Exact native values are calculated directly, while midpoint-scaled "
        "unit values are contextual approximations."
    )
    method = st.radio("Method", ["M11", "M02", "M05"], horizontal=True)
    method_context = results["ba_context"][results["ba_context"]["method"] == method].copy()
    st.dataframe(pf_styler(method_context), use_container_width=True, hide_index=True)
    with st.expander("Full percentage results"):
        st.dataframe(pf_styler(results["ba_percent"]), use_container_width=True, hide_index=True)
    with st.expander("Full exact native-unit results"):
        st.dataframe(results["ba_native"], use_container_width=True, hide_index=True)
    with st.expander("Acceptance criteria used"):
        st.dataframe(results["criteria"], use_container_width=True, hide_index=True)

with audit_tab:
    st.subheader("Bland–Altman replicate removals")
    st.dataframe(results["ba_removed"], use_container_width=True, hide_index=True)
    with st.expander("Full Shapiro / Grubbs / MAD audit"):
        st.dataframe(results["ba_audit"], use_container_width=True, hide_index=True)
    with st.expander("Trueness generalized-ESD removals"):
        st.dataframe(results["regression_outliers"], use_container_width=True, hide_index=True)
    with st.expander("Trueness ESD diagnostics"):
        st.dataframe(results["regression_diagnostics"], use_container_width=True, hide_index=True)
    with st.expander("Global-flag and unparsed exclusions"):
        st.dataframe(results["global_exclusions"], use_container_width=True, hide_index=True)
    with st.expander("Analyte-specific flag exclusions"):
        st.dataframe(results["ba_analyte_flags"], use_container_width=True, hide_index=True)

with plot_tab:
    plot_names = list(results["plots"])
    if not plot_names:
        st.info("No plot files were generated.")
    else:
        category = st.selectbox(
            "Plot category",
            ["All", "Regression", "M02", "M11", "M05", "Donor profiles"],
        )
        filtered_names = plot_names
        if category == "Regression":
            filtered_names = [name for name in plot_names if "regression" in name.lower()]
        elif category in {"M02", "M11", "M05"}:
            filtered_names = [name for name in plot_names if category.lower() in name.lower() and "bland_altman" in name.lower()]
        elif category == "Donor profiles":
            filtered_names = [name for name in plot_names if "donor_profiles" in name.lower()]
        chosen_plot = st.selectbox("Saved PNG", filtered_names, key="all_plot_selector")
        st.image(results["plots"][chosen_plot], caption=chosen_plot, use_container_width=True)
        st.download_button(
            "Download selected PNG",
            data=results["plots"][chosen_plot],
            file_name=Path(chosen_plot).name,
            mime="image/png",
        )

with download_tab:
    st.download_button(
        "Download complete results ZIP",
        data=results["zip"],
        file_name="PROXIMA_COMBINED_RESULTS.zip",
        mime="application/zip",
        use_container_width=True,
    )
    st.download_button(
        "Download combined Excel workbook",
        data=results["xlsx"],
        file_name="PROXIMA_Trueness_and_BlandAltman_Output.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        use_container_width=True,
    )
    st.download_button(
        "Download main percent/native summary CSV",
        data=results["summary_csv"],
        file_name="PROXIMA_BA_MAIN_PERCENT_NATIVE_SUMMARY.csv",
        mime="text/csv",
        use_container_width=True,
    )
    st.download_button(
        "Download canonicalized input used by the analysis",
        data=results["canonical_csv"],
        file_name="canonical_uploaded_input.csv",
        mime="text/csv",
        use_container_width=True,
    )
    with st.expander("Output inventory"):
        st.code("\n".join(results["inventory"]))

if st.button("Clear stored results"):
    st.session_state.pop("proxima_final_results", None)
    st.rerun()
