from __future__ import annotations

import io
import inspect
import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

APP_BUILD = "PROXIMA V6.1 — matched core import fix"
REQUIRED_CORE_API_VERSION = "2026-07-24-proxima-v6"

st.set_page_config(
    page_title="PROXIMA Trueness + Bland–Altman",
    page_icon="🧪",
    layout="wide",
)

# Import the core as a module first instead of importing many names directly.
# This prevents the redacted Streamlit Cloud ImportError shown when app.py and
# analysis_core.py come from different releases, and gives an actionable error.
try:
    import analysis_core as _analysis_core
except Exception as exc:
    st.error(
        "The statistical engine could not be imported. Upload app.py and "
        "analysis_core.py from the same fixed package, then reboot the app."
    )
    st.code(f"{type(exc).__name__}: {exc}")
    st.stop()

_REQUIRED_CORE_SYMBOLS = [
    "DEFAULT_ANALYTE_CONFIG",
    "build_runtime_analytes",
    "canonicalize_input_dataframe",
    "collect_output_files",
    "default_analyte_table",
    "normalize_bool",
    "normalize_specimen_label",
    "parse_blood_sample_id",
    "read_binary",
    "run_complete_notebook_pipeline",
]
_missing_core_symbols = [
    name for name in _REQUIRED_CORE_SYMBOLS if not hasattr(_analysis_core, name)
]
_loaded_core_version = getattr(_analysis_core, "CORE_API_VERSION", None)
if _missing_core_symbols or _loaded_core_version != REQUIRED_CORE_API_VERSION:
    st.error(
        "app.py and analysis_core.py are from different releases. Replace both "
        "files together and reboot the Streamlit app."
    )
    st.code(
        "Expected core version: " + REQUIRED_CORE_API_VERSION
        + "\nLoaded core version: " + str(_loaded_core_version)
        + "\nMissing symbols: "
        + (", ".join(_missing_core_symbols) if _missing_core_symbols else "none")
    )
    st.stop()

DEFAULT_ANALYTE_CONFIG = _analysis_core.DEFAULT_ANALYTE_CONFIG
build_runtime_analytes = _analysis_core.build_runtime_analytes
canonicalize_input_dataframe = _analysis_core.canonicalize_input_dataframe
collect_output_files = _analysis_core.collect_output_files
default_analyte_table = _analysis_core.default_analyte_table
normalize_bool = _analysis_core.normalize_bool
normalize_specimen_label = _analysis_core.normalize_specimen_label
parse_blood_sample_id = _analysis_core.parse_blood_sample_id
read_binary = _analysis_core.read_binary
run_complete_notebook_pipeline = _analysis_core.run_complete_notebook_pipeline

st.title("🧪 PROXIMA Trueness + Bland–Altman App")
st.caption(f"Build: {APP_BUILD}")
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
- **Optional sensitivity run:** the user may disable both replicate-level outlier branches; all otherwise eligible rows are then retained before donor averaging, with this choice recorded in every audit/configuration output.
- **Donor parsing:** the app extracts the D-token, collection prefix, and optional explicit Donor column, then asks the user which identity rule to use. This prevents accidental merging when the same parsed token occurs under different prefixes or maps to labels such as `D03` and `D03b`.
- **Bland–Altman is optional:** trueness regression runs independently. With one specimen type, the optional `REF` mode compares MHS directly with its matched Sysmex `_ref` value.
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


def donor_sort_key(value: object):
    text = str(value)
    match = __import__("re").fullmatch(r"([A-Za-z]+)(\d+)([A-Za-z0-9]*)", text)
    if match:
        prefix, number, suffix = match.groups()
        return (prefix.lower(), int(number), suffix.lower())
    return (text.lower(), 0, "")


def metadata_preview(
    data: pd.DataFrame,
    *,
    sample_id_col: str | None,
    donor_col: str | None,
    specimen_col: str | None,
    replicate_col: str | None,
    donor_identity_mode: str = "parsed_token",
    verification_donor_col: str | None = None,
) -> pd.DataFrame:
    if sample_id_col is not None:
        parsed = pd.DataFrame([parse_blood_sample_id(v) for v in data[sample_id_col]])
        parsed.insert(0, "source_row", np.arange(1, len(parsed) + 1))
        parsed.insert(1, "source_sample_id", data[sample_id_col].astype(str).values)
        explicit = None
        if verification_donor_col is not None and verification_donor_col in data.columns:
            explicit = data[verification_donor_col].astype(str).str.strip().replace({"": np.nan, "nan": np.nan})
            parsed["explicit_donor_value"] = explicit.values
        mode = str(donor_identity_mode).strip().lower()
        if mode == "explicit_column" and explicit is not None:
            parsed["donor"] = explicit.where(explicit.notna(), parsed["parsed_donor_token"]).values
        elif mode == "prefix_plus_token":
            parsed["donor"] = parsed["prefix_donor_key"]
        else:
            parsed["donor"] = parsed["parsed_donor_token"]
        parsed["specimen_type"] = parsed["specimen_type"].map(normalize_specimen_label)
        parsed["parse_ok"] = parsed["parse_ok"].fillna(False) & parsed["donor"].notna()
        parsed["donor_identity_mode"] = mode
        return parsed
    preview = pd.DataFrame({
        "source_row": np.arange(1, len(data) + 1),
        "donor": data[donor_col].astype(str).str.strip(),
        "specimen_type": data[specimen_col].map(normalize_specimen_label),
        "replicate_number": pd.to_numeric(data[replicate_col], errors="coerce"),
    })
    preview["parsed_donor_token"] = preview["donor"]
    preview["sample_prefix"] = ""
    preview["prefix_donor_key"] = preview["donor"]
    preview["donor_identity_mode"] = "explicit_columns"
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
verification_donor_col: str | None = None
donor_identity_mode = "parsed_token"
donor_mapping_reviewed = True
prefix_collision_count = 0
explicit_conflict_count = 0

if id_mode.startswith("Parse"):
    sample_id_col = st.selectbox(
        "Sample-ID column (expected ending such as D05-Mc-2)",
        columns,
        index=option_index(columns, "bloodSampleId"),
    )
    preliminary = pd.DataFrame([parse_blood_sample_id(v) for v in data[sample_id_col]])
    donor_candidates = ["<none>"] + columns
    default_donor_column = "Donor" if "Donor" in columns else next(
        (c for c in columns if str(c).strip().lower() in {"donor", "donor_id", "donorid"}),
        "<none>",
    )
    verification_choice = st.selectbox(
        "Optional donor-verification column",
        donor_candidates,
        index=option_index(donor_candidates, default_donor_column),
        help=(
            "When present, this column can preserve distinctions such as D03 versus D03b even when the "
            "sample-ID text contains the same parsed D03 token under different collection prefixes."
        ),
    )
    verification_donor_col = None if verification_choice == "<none>" else verification_choice

    prefix_counts = (
        preliminary.dropna(subset=["parsed_donor_token", "prefix_donor_key"])
        .groupby("parsed_donor_token")["prefix_donor_key"].nunique()
    )
    prefix_collision_count = int((prefix_counts > 1).sum())
    identity_labels = {
        "Use explicit Donor column (recommended when available)": "explicit_column",
        "Use parsed D-token only (for example D03)": "parsed_token",
        "Use collection prefix + D-token (for example 130726D03)": "prefix_plus_token",
    }
    identity_options = list(identity_labels)
    if verification_donor_col is None:
        identity_options = identity_options[1:]
        default_identity_index = 1 if prefix_collision_count else 0
    else:
        default_identity_index = 0
    identity_label = st.radio(
        "Donor identity rule",
        identity_options,
        index=default_identity_index,
        help=(
            "The app always parses the D-token and collection prefix. Choose the explicit Donor column when it "
            "contains reviewed labels such as D03b; otherwise prefix + D-token prevents accidental merging "
            "when the same D-token occurs under different date/collection prefixes."
        ),
    )
    donor_identity_mode = identity_labels[identity_label]

    if verification_donor_col is not None:
        explicit = data[verification_donor_col].astype(str).str.strip().replace({"": np.nan, "nan": np.nan})
        comparison = pd.DataFrame({
            "parsed": preliminary["parsed_donor_token"],
            "prefix_key": preliminary["prefix_donor_key"],
            "explicit": explicit,
        })
        explicit_map_counts = comparison.dropna(subset=["parsed", "explicit"]).groupby("parsed")["explicit"].nunique()
        explicit_conflict_count = int((explicit_map_counts > 1).sum())

    try:
        parsed_preview = metadata_preview(
            data,
            sample_id_col=sample_id_col,
            donor_col=None,
            specimen_col=None,
            replicate_col=None,
            donor_identity_mode=donor_identity_mode,
            verification_donor_col=verification_donor_col,
        )
    except Exception as exc:
        st.error(f"Could not inspect specimen IDs: {exc}")
        st.stop()

    review_needed = prefix_collision_count > 0 or explicit_conflict_count > 0
    if review_needed:
        st.warning(
            f"Donor-identity review required: {prefix_collision_count} parsed D-token(s) occur under multiple "
            f"prefix keys and {explicit_conflict_count} parsed D-token(s) map to multiple explicit donor labels."
        )
        donor_mapping_reviewed = st.checkbox(
            "I reviewed the donor identity mapping below and confirm the selected donor rule",
            value=False,
        )
    with st.expander("Verify donor identity mapping", expanded=review_needed):
        mapping_cols = [
            c for c in ["source_sample_id", "sample_prefix", "parsed_donor_token", "prefix_donor_key",
                        "explicit_donor_value", "donor", "specimen_type", "replicate_number", "parse_ok"]
            if c in parsed_preview.columns
        ]
        mapping = parsed_preview[mapping_cols].drop_duplicates().sort_values(
            [c for c in ["donor", "specimen_type", "replicate_number"] if c in mapping_cols],
            kind="stable",
        )
        st.dataframe(mapping, use_container_width=True, hide_index=True)
else:
    id1, id2, id3 = st.columns(3)
    donor_col = id1.selectbox("Donor column", columns, index=option_index(columns, "Donor"))
    specimen_col = id2.selectbox("Specimen-type column", columns, index=option_index(columns, "Type"))
    replicate_col = id3.selectbox("Replicate-number column", columns, index=option_index(columns, "Replicate"))
    try:
        parsed_preview = metadata_preview(
            data,
            sample_id_col=None,
            donor_col=donor_col,
            specimen_col=specimen_col,
            replicate_col=replicate_col,
            donor_identity_mode="explicit_columns",
        )
    except Exception as exc:
        st.error(f"Could not inspect specimen IDs: {exc}")
        st.stop()

try:
    detected_specimen_types = sorted(
        parsed_preview.loc[parsed_preview["parse_ok"], "specimen_type"].dropna().astype(str).unique().tolist()
    )
    detected_donors = sorted(
        parsed_preview.loc[parsed_preview["parse_ok"], "donor"].dropna().astype(str).unique().tolist(),
        key=donor_sort_key,
    )
    donor_count = len(detected_donors)
    parse_failures = int((~parsed_preview["parse_ok"].fillna(False)).sum())
except Exception as exc:
    st.error(f"Could not inspect specimen IDs: {exc}")
    st.stop()

m1, m2, m3 = st.columns(3)
m1.metric("Detected unique donors", donor_count)
m2.metric("Detected specimen types", len(detected_specimen_types))
m3.metric("Unparsed rows", parse_failures)
donor_rule_text = {
    "explicit_column": "the selected explicit Donor column, with parsed IDs retained for verification",
    "prefix_plus_token": "collection prefix + parsed D-token (for example 130726D03)",
    "parsed_token": "the parsed D-token up to the first hyphen (for example D03 or D03b)",
    "explicit_columns": "the selected explicit donor/specimen/replicate columns",
}.get(donor_identity_mode, donor_identity_mode)
st.caption(f"Final donor identity uses {donor_rule_text}. All detected specimen labels are normalized and shown for user selection.")
if parse_failures:
    st.warning("Unparsed rows are excluded and listed in the selection audit.")

select1, select2 = st.columns(2)
selected_donors = select1.multiselect(
    "Donors to include",
    options=detected_donors,
    default=detected_donors,
    help="All detected donors are selected by default. Unselect only donors you intentionally want to exclude.",
)
selected_specimen_types = select2.multiselect(
    "Specimen/sample types to include",
    options=detected_specimen_types,
    default=detected_specimen_types,
    help="Automatically detected labels such as Mc, Mc2, MK2, MK3, Dir, Tas or Tasso can be selected independently.",
)
if not selected_donors:
    st.warning("Select at least one donor.")
    st.stop()
if not selected_specimen_types:
    st.warning("Select at least one specimen/sample type.")
    st.stop()

specimen_types = list(selected_specimen_types)
selected_preview_mask = (
    parsed_preview["parse_ok"].fillna(False)
    & parsed_preview["donor"].astype(str).isin(selected_donors)
    & parsed_preview["specimen_type"].astype(str).isin(selected_specimen_types)
)
selected_preview = parsed_preview.loc[selected_preview_mask].copy()

with st.expander("Detected donor × specimen audit", expanded=False):
    donor_specimen_table = (
        parsed_preview.loc[parsed_preview["parse_ok"].fillna(False)]
        .groupby(["donor", "specimen_type"], dropna=False)
        .size().rename("row_count").reset_index()
        .pivot(index="donor", columns="specimen_type", values="row_count").fillna(0).astype(int)
        .reindex(detected_donors)
    )
    donor_specimen_table.insert(0, "Selected donor", donor_specimen_table.index.astype(str).isin(selected_donors))
    st.dataframe(donor_specimen_table, use_container_width=True)

if global_choice != "<no column>":
    global_bool = data[global_choice].map(normalize_bool)
    false_n = int(global_bool.eq(False).sum())
    true_n = int(global_bool.eq(True).sum())
    unknown_n = int(global_bool.isna().sum())
    selected_global_false_donors = int(
        selected_preview.loc[global_bool.loc[selected_preview.index].eq(False), "donor"].nunique()
    ) if not selected_preview.empty else 0
    st.caption(
        f"Global-flag preview: FALSE={false_n}, TRUE={true_n}, missing/unparseable={unknown_n}. "
        f"Among the current donor/specimen selection, {selected_global_false_donors} donors have at least one global_flag=FALSE row."
    )

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

ignore_outlier_removal = st.checkbox(
    "Ignore replicate-level outlier removal (sensitivity run)",
    value=False,
    help=(
        "Default OFF preserves the validated automatic ESD and Grubbs/MAD procedures. "
        "Turn this ON only to retain every otherwise eligible replicate before donor averaging. "
        "Global and analyte-specific flag exclusions still apply."
    ),
)
apply_outlier_removal = not ignore_outlier_removal
if ignore_outlier_removal:
    st.warning("Sensitivity mode: replicate-level outlier removal is disabled; all otherwise eligible replicates will be retained.")

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

st.header("5. Optional Bland–Altman analysis")
if len(specimen_types) == 1:
    st.info(
        "One specimen type was detected. Leave Bland–Altman unchecked to run trueness regression only, "
        "or check it to run single-specimen MHS-versus-matched-Sysmex (_ref) Bland–Altman analysis."
    )
else:
    st.caption(
        "Bland–Altman is optional. Leave it unchecked for regression-only analysis; check it only when BA outputs are required."
    )

run_bland_altman = st.checkbox(
    "Also run Bland–Altman analysis",
    value=False,
    key="run_bland_altman_optional_v6",
    help=(
        "Optional and off by default. Trueness regression runs independently. With two or more specimen types, choose paired "
        "M02/M11/M05; with one specimen type, compare MHS directly with its matched Sysmex _ref column."
    ),
)

ba_mode = "disabled"
ba_a_types: list[str] = []
ba_b_types: list[str] = []
ba_reference_types: list[str] = []
ba_a_label = "Specimen A"
ba_b_label = "Specimen B"
criteria_confirmed = False

if run_bland_altman:
    mode_options = ["MHS versus matched Sysmex reference within selected specimen"]
    if len(specimen_types) >= 2:
        mode_options.insert(0, "Paired specimens: M02, M11 and M05")
    ba_mode_label = st.radio("Bland–Altman mode", mode_options, horizontal=True)

    if ba_mode_label.startswith("Paired"):
        ba_mode = "paired_specimens"
        default_a = [value for value in ["Mc", "Mc2"] if value in specimen_types] or [specimen_types[0]]
        default_b = ["MK2"] if "MK2" in specimen_types else [value for value in specimen_types if value not in default_a][:1]
        ba1, ba2 = st.columns(2)
        ba_a_types = ba1.multiselect("Specimen A types", specimen_types, default=default_a)
        ba_b_types = ba2.multiselect("Specimen B types", specimen_types, default=default_b)
        if ba_a_types and ba_b_types:
            selected_pair_preview = selected_preview.copy()
            donor_type_sets = selected_pair_preview.groupby("donor")["specimen_type"].agg(lambda x: set(map(str, x)))
            pair_eligible = donor_type_sets.map(lambda values: bool(values & set(ba_a_types)) and bool(values & set(ba_b_types)))
            pair1, pair2, pair3 = st.columns(3)
            pair1.metric("Selected donors", len(selected_donors))
            pair2.metric("BA pair-eligible before flags", int(pair_eligible.sum()))
            pair3.metric("Not pair-eligible before flags", int((~pair_eligible).sum()))
            if int((~pair_eligible).sum()) > 0:
                missing_names = ", ".join(pair_eligible.index[~pair_eligible].astype(str).tolist())
                st.info(
                    "Some donors lack at least one selected BA arm and cannot enter paired M02/M11/M05 analysis: "
                    + missing_names
                )
        label1, label2 = st.columns(2)
        ba_a_label = label1.text_input("Specimen A display label", value="Capillary")
        ba_b_label = label2.text_input("Specimen B display label", value="MK2")
    else:
        ba_mode = "single_specimen_reference"
        default_ref = ["Mc"] if "Mc" in specimen_types else [specimen_types[0]]
        ba_reference_types = st.multiselect(
            "Specimen type(s) for MHS-versus-Sysmex reference BA",
            specimen_types,
            default=default_ref,
            help=(
                "For each donor and analyte, retained replicates are averaged and MHS is compared with the matched "
                "_ref/Sysmex value. Percentage bias = 100×(MHS−Sysmex)/Sysmex; native bias = MHS−Sysmex."
            ),
        )
        ba_a_label = st.text_input("Selected specimen display label", value=", ".join(ba_reference_types) or "Selected specimen")
        donor_type_sets = selected_preview.groupby("donor")["specimen_type"].agg(lambda x: set(map(str, x)))
        ref_eligible = donor_type_sets.map(lambda values: bool(values & set(ba_reference_types))) if ba_reference_types else pd.Series(False, index=donor_type_sets.index)
        ref1, ref2, ref3 = st.columns(3)
        ref1.metric("Selected donors", len(selected_donors))
        ref2.metric("Reference-BA eligible before flags", int(ref_eligible.sum()))
        ref3.metric("Not eligible before flags", int((~ref_eligible).sum()))

    criteria_confirmed = st.checkbox(
        "I reviewed and verified the edited AC/CLIA table for this analysis",
        value=False,
    )
    if not criteria_confirmed:
        st.warning("Pass/Fail labels will be generated as **UNVERIFIED / PROVISIONAL CONTEXT ONLY**.")
else:
    st.info("Bland–Altman is disabled for this run. Trueness regression and its plots/exports will still run normally.")

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
        if not donor_mapping_reviewed:
            raise ValueError("Review and confirm the donor identity mapping before running the analysis.")
        if run_bland_altman and ba_mode == "paired_specimens":
            if not ba_a_types or not ba_b_types:
                raise ValueError("Select at least one specimen type in both paired Bland–Altman arms.")
            overlap = sorted(set(ba_a_types) & set(ba_b_types))
            if overlap:
                raise ValueError("The paired Bland–Altman specimen arms must not overlap: " + ", ".join(overlap))
        if run_bland_altman and ba_mode == "single_specimen_reference" and not ba_reference_types:
            raise ValueError("Select at least one specimen type for MHS-versus-reference Bland–Altman analysis.")
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
            donor_identity_mode=donor_identity_mode,
            verification_donor_col=verification_donor_col,
        )

        with st.status("Running the notebook pipeline from beginning to end…", expanded=True) as status:
            st.write("Applying global-flag and ID parsing rules…")
            st.write(
                "Running replicate-level generalized-ESD trueness regressions…"
                if apply_outlier_removal else
                "Running trueness regressions with outlier removal disabled…"
            )
            if run_bland_altman:
                if ba_mode == "paired_specimens":
                    ba_status_text = "Running paired M02/M11/M05 Bland–Altman analyses"
                else:
                    ba_status_text = "Running MHS-versus-Sysmex reference Bland–Altman analysis"
                st.write(
                    ba_status_text + ", outlier screening and bootstraps…"
                    if apply_outlier_removal else
                    ba_status_text + " and bootstraps with outlier removal disabled…"
                )
            else:
                st.write("Bland–Altman disabled; continuing with trueness regression only…")
            st.write("Saving all plots, CSV audits, Excel sheets, and ZIP outputs…")
            with tempfile.TemporaryDirectory(prefix="proxima_app_") as temp_dir:
                pipeline_kwargs = {
                    "runtime_analytes": runtime_analytes,
                    "ac_limits": ac_limits,
                    "clia_limits": clia_limits,
                    "regression_groups": regression_groups,
                    "groups_to_run": groups_to_run,
                    "ba_specimen_a_types": ba_a_types,
                    "ba_specimen_b_types": ba_b_types,
                    "ba_specimen_a_label": ba_a_label.strip() or "Specimen A",
                    "ba_specimen_b_label": ba_b_label.strip() or "Specimen B",
                    "criteria_confirmed": criteria_confirmed,
                    "run_bland_altman": run_bland_altman,
                    "ba_mode": ba_mode,
                    "ba_reference_specimen_types": ba_reference_types,
                    "run_huber_bootstrap": run_huber_bootstrap,
                    "n_huber_bootstraps": int(n_huber_bootstraps),
                    "plot_mode": plot_mode_value,
                    "shapiro_alpha": float(shapiro_alpha),
                    "outlier_alpha": float(outlier_alpha),
                    "mad_z_threshold": float(mad_z),
                    "ba_bootstrap_iterations": int(ba_bootstraps),
                    "donor_profile_bootstraps": int(donor_profile_bootstraps),
                    "random_seed": int(random_seed),
                    "cross_specimen_groups": cross_pair,
                    "regression_equation_position": regression_equation_position,
                    "regression_metrics_position": regression_metrics_position,
                    "apply_outlier_removal": apply_outlier_removal,
                    "selected_donors": selected_donors,
                    "selected_specimen_types": selected_specimen_types,
                }

                # Version-safe call: the current analysis_core.py supports both
                # annotation-position arguments. Filtering prevents a stale
                # cached/deployed core from crashing with an unexpected-keyword
                # TypeError and reports exactly which file needs replacement.
                core_signature = inspect.signature(run_complete_notebook_pipeline)
                accepts_extra_kwargs = any(
                    parameter.kind == inspect.Parameter.VAR_KEYWORD
                    for parameter in core_signature.parameters.values()
                )
                unsupported_kwargs = [
                    key for key in pipeline_kwargs
                    if key not in core_signature.parameters and not accepts_extra_kwargs
                ]
                if unsupported_kwargs:
                    st.warning(
                        "The loaded analysis_core.py is older than app.py and does not "
                        "support: " + ", ".join(unsupported_kwargs) + ". "
                        "The analysis will continue with those display-only options "
                        "omitted; replace both app.py and analysis_core.py together."
                    )
                    pipeline_kwargs = {
                        key: value for key, value in pipeline_kwargs.items()
                        if key in core_signature.parameters or accepts_extra_kwargs
                    }

                result = run_complete_notebook_pipeline(
                    canonical,
                    temp_dir,
                    **pipeline_kwargs,
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
                    "canonical_csv": result["canonical_selected"].to_csv(index=False).encode("utf-8"),
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
                    "ba_enabled": bool(result["ba"].get("ba_enabled", run_bland_altman)),
                    "ba_mode": str(result["ba"].get("ba_mode", ba_mode)),
                    "donor_identity_audit": result["donor_identity_audit"].copy(),
                    "donor_specimen_audit": result["donor_specimen_audit"].copy(),
                    "selection_exclusions": result["selection_exclusions"].copy(),
                    "ba_pairing_audit": result["ba_pairing_audit"].copy(),
                    "plots": plot_bytes,
                    "inventory": inventory,
                    "input_rows": len(result["canonical_selected"]),
                    "global_false_rows": int(result["canonical_selected"]["global_flag"].map(normalize_bool).eq(False).sum()),
                    "detected_donor_count": len(detected_donors),
                    "selected_donor_count": len(selected_donors),
                    "pair_eligible_donor_count": int(result["ba_pairing_audit"]["BA_pair_eligible_before_flags"].sum()) if not result["ba_pairing_audit"].empty else 0,
                    "selected_analytes": list(runtime_analytes),
                    "criteria_confirmed": criteria_confirmed,
                    "apply_outlier_removal": apply_outlier_removal,
                    "donor_identity_mode": donor_identity_mode,
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
summary1, summary2, summary3, summary4, summary5, summary6 = st.columns(6)
summary1.metric("Selected input rows", f"{results['input_rows']:,}")
summary2.metric("global_flag = FALSE", f"{results['global_false_rows']:,}")
summary3.metric("Detected donors", results["detected_donor_count"])
summary4.metric("Selected donors", results["selected_donor_count"])
if results.get("ba_enabled", False):
    eligible_label = "BA eligible" if results.get("ba_mode") == "single_specimen_reference" else "BA pair-eligible"
    summary5.metric(eligible_label, results["pair_eligible_donor_count"])
    summary6.metric("BA outliers removed", len(results["ba_removed"]))
else:
    summary5.metric("Bland–Altman", "Not run")
    summary6.metric("BA outliers removed", "—")
st.caption(
    "Outlier removal: " + ("enabled (validated default)" if results["apply_outlier_removal"] else "disabled by user for sensitivity analysis")
)

if results.get("ba_enabled", False) and not results["criteria_confirmed"]:
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
    if not results.get("ba_enabled", False):
        st.info("Bland–Altman was not requested. Trueness regression results and exports remain complete.")
    elif results.get("ba_mode") == "single_specimen_reference":
        st.info(
            "REF compares donor-mean MHS with its matched Sysmex _ref value within the selected specimen type. "
            "Percentage bias is 100×(MHS−Sysmex)/Sysmex and native bias is MHS−Sysmex."
        )
        st.dataframe(pf_styler(results["ba_context"]), use_container_width=True, hide_index=True)
        with st.expander("Full percentage results"):
            st.dataframe(pf_styler(results["ba_percent"]), use_container_width=True, hide_index=True)
        with st.expander("Full exact native-unit results"):
            st.dataframe(results["ba_native"], use_container_width=True, hide_index=True)
        with st.expander("Acceptance criteria used"):
            st.dataframe(results["criteria"], use_container_width=True, hide_index=True)
    else:
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
    st.subheader("Donor identity, specimen detection, and retention")
    with st.expander("Donor identity mapping audit", expanded=True):
        st.dataframe(results["donor_identity_audit"], use_container_width=True, hide_index=True)
    st.dataframe(results["ba_pairing_audit"], use_container_width=True, hide_index=True)
    with st.expander("Full donor × specimen detection audit"):
        st.dataframe(results["donor_specimen_audit"], use_container_width=True, hide_index=True)
    with st.expander("Rows excluded by donor/specimen selection or parsing"):
        st.dataframe(results["selection_exclusions"], use_container_width=True, hide_index=True)

    if results.get("ba_enabled", False):
        st.subheader("Bland–Altman replicate removals")
        st.dataframe(results["ba_removed"], use_container_width=True, hide_index=True)
        with st.expander("Full Shapiro / Grubbs / MAD audit"):
            st.dataframe(results["ba_audit"], use_container_width=True, hide_index=True)
        with st.expander("Global-flag and unparsed exclusions"):
            st.dataframe(results["global_exclusions"], use_container_width=True, hide_index=True)
        with st.expander("Analyte-specific flag exclusions"):
            st.dataframe(results["ba_analyte_flags"], use_container_width=True, hide_index=True)
    else:
        st.info("Bland–Altman was disabled, so no BA-specific outlier or analyte-flag audit was generated.")

    with st.expander("Trueness generalized-ESD removals"):
        st.dataframe(results["regression_outliers"], use_container_width=True, hide_index=True)
    with st.expander("Trueness ESD diagnostics"):
        st.dataframe(results["regression_diagnostics"], use_container_width=True, hide_index=True)

with plot_tab:
    plot_names = list(results["plots"])
    if not plot_names:
        st.info("No plot files were generated.")
    else:
        plot_categories = ["All", "Regression"]
        if results.get("ba_enabled", False):
            if results.get("ba_mode") == "single_specimen_reference":
                plot_categories.append("REF")
            else:
                plot_categories.extend(["M02", "M11", "M05", "Donor profiles"])
        category = st.selectbox("Plot category", plot_categories)
        filtered_names = plot_names
        if category == "Regression":
            filtered_names = [name for name in plot_names if "regression" in name.lower()]
        elif category in {"M02", "M11", "M05", "REF"}:
            filtered_names = [name for name in plot_names if category.lower() in name.lower() and "bland_altman" in name.lower()]
        elif category == "Donor profiles":
            filtered_names = [name for name in plot_names if "donor_profiles" in name.lower()]
        if not filtered_names:
            st.info("No plots are available in this category.")
        else:
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
