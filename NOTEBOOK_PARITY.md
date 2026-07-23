# Notebook parity checklist

The app core was built directly from the code cells of
`PROXIMA_Trueness_and_BlandAltman_v2.ipynb`.

## Preserved functions and outputs

- validated generalized ESD transformation and branching;
- analyte-specific replicate removal before donor averaging;
- Huber regression;
- Pearson r and Fisher-z 95% CI;
- normal-range red lines and the same regression annotations;
- M02, M11, M05, and S02 equations;
- Shapiro-guided manual Grubbs/MAD outlier branch;
- 10,000-iteration BA donor bootstrap by default;
- 5,000-iteration within-donor profile bootstrap by default;
- percentage and exact native-unit summaries;
- normal midpoint and 50%-midpoint context;
- mean-bias, LoA point, and full endpoint-uncertainty AC/CLIA screens;
- individual regression plots and panel;
- individual M02/M11/M05 percentage/native plots and panels;
- percentage/native donor-profile plots and Appendix B panels;
- CSV audits and combined Excel workbook;
- complete results ZIP.

The app adds only input/configuration controls and download/display components.
