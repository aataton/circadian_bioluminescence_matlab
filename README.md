# MATLAB Circadian Bioluminescence Analysis Pipeline

A configurable MATLAB pipeline for analysing circadian bioluminescence time-series data. The workflow reads plate-based luminescence measurements, groups replicate wells defined in plain-text parameter files, interpolates missing values, optionally detrends the data, fits chirped damped cosine models to user-defined time segments, and exports per-sample results and publication-ready figures.

The repository also includes configurable workflows to:

- assemble mean ± SD luminescence traces into overlay or tiled summary figures;
- quantify period-versus-temperature relationships and calculate \(Q_{10}\); and
- generate phase-response curves (PRCs) from reference and pulse-condition acrophase data.

This code is configured and tested with **MATLAB R2025a**.

> **Important:** This repository should contain scripts, MATLAB functions, shell launchers, raw input data, and parameter files needed for reproducibility. Generated figures, fitted-result CSVs, MATLAB `.fig` files, and other analysis products are intentionally not required in the repository: they can be regenerated from the included inputs by following the workflow below.

## Contents

| File | Role |
|---|---|
| `Plot_rawLuminescence.m` | Main single-condition analysis script. Reads one raw CSV file plus one parameter file, calculates mean and SD, runs fits, and exports outputs. |
| `processSamples.m` | Reads and validates the 13-line (optional 14th line) parameter-file format. |
| `analyzeAndPlotSegments.m` | Performs per-segment chirped damped-cosine fitting and writes period, phase, damping, chirp, peak, and drift outputs. |
| `fitCosineCurve.m` | FFT-initialized, nonlinear least-squares fitting engine for a chirped damped cosine. |
| `_RunMatlabFronTerm_*.sh` | Example Bash launchers that run `Plot_rawLuminescence.m` in MATLAB batch mode for multiple conditions. |
| `Plot_summary.m` | User-editable configuration script for generating multi-series overlays and tiled summary plots. |
| `run_luminescence_plots.m` | Summary-plot engine called by `Plot_summary.m`. |
| `Q10_config.m` | User-editable configuration script for period-versus-temperature and \(Q_{10}\) figures. |
| `Q10_plot.m` | Plotting and calculation engine called by `Q10_config.m`. |
| `PRC_config.m` | User-editable configuration script for phase-response-curve figures. |
| `PRC_plot.m` | PRC calculation and plotting engine called by `PRC_config.m`. |

The exact example folder and file names in the repository are dataset-specific. For a new experiment, retain the relative-path convention or update the relevant shell launcher/configuration file accordingly.

## Requirements

### Required software

- MATLAB **R2025a**.
- A Unix-like shell environment for the provided `.sh` batch-launch scripts (Linux, macOS, or Windows through WSL/Git Bash, as appropriate).
- MATLAB available on your shell `PATH` as `matlab` when using the `.sh` scripts.

### Required MATLAB products

The fitting workflow relies on functions including `fit`, `fittype`, `fitoptions`, `findpeaks`, `detrend`, and `tcdf`. Install and license the following products:

- MATLAB
- Curve Fitting Toolbox
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

A quick environment check from the MATLAB Command Window is:

```matlab
ver
which fit -all
which findpeaks -all
which tcdf -all
```

The listed functions should resolve to installed MATLAB products rather than returning an empty result.

### No external MATLAB packages

The repository does not require third-party MATLAB packages, Python, R, or internet access to execute the pipeline. All project `.m` files must remain together in the repository root, or that root must be added to the MATLAB path before execution.

## Installation and setup

1. Clone or download the repository.

   ```bash
   git clone <REPOSITORY-URL>
   cd <REPOSITORY-DIRECTORY>
   ```

2. Confirm that MATLAB is callable from the terminal.

   ```bash
   matlab -batch "disp(version)"
   ```

3. Retain the project directory layout. In particular, keep raw CSV data and its parameter files in the relative locations referenced by the shell scripts, `Plot_summary.m`, `Q10_config.m`, and `PRC_config.m`.

4. Make the launcher scripts executable if needed.

   ```bash
   chmod +x _RunMatlabFronTerm_*.sh
   ```

5. For interactive MATLAB use, open MATLAB from the repository root or run:

   ```matlab
   addpath(pwd)
   ```

## Input data

### Raw luminescence CSV

`Plot_rawLuminescence.m` expects a comma-separated file with:

- a first column named exactly `Time`, in hours; and
- one or more subsequent numeric columns, each representing an individual sample, well, or replicate.

The names in the first line of the parameter file must match the relevant CSV variable names exactly. MATLAB-valid variable names are strongly recommended; avoid spaces and punctuation in raw column headers.

Example:

```csv
Time,Sample_00001,Sample_00002,Sample_00003
0,11250,10830,11705
1,10640,10325,10915
2,9870,9610,10020
```

The main workflow sorts data by `Time`. Missing values are filled by linear interpolation with extrapolation when a sample contains at least two non-missing observations. Rows that still contain a missing value in any selected sample after this procedure are removed. Values used for the mean/SD and fitting workflow are rounded to whole numbers; optional linear detrending is applied after interpolation.

### Parameter file

Each analysis group requires a plain-text parameter file. The parameter file selects sample columns, defines output names and fitting windows, and controls plotting behavior. It contains **13 required lines** and an **optional 14th line**.

| Line | Setting | Format and meaning |
|---:|---|---|
| 1 | `samplesList` | Comma-separated raw-CSV column names to include in the group. |
| 2 | `sampleName` | Output prefix, optionally including a relative folder path. |
| 3 | `segments` | Semicolon-separated fitting windows, each given as `start end` in hours. Example: `96 192; 192 288`. |
| 4 | `detrenddata` | `1` to linearly detrend each trace before analysis; `0` to retain the original trend. |
| 5 | `resettingHours` | Four polygon x-coordinates defining a dark/resetting patch: `x1 x2 x2 x1`; use `0 0 0 0` to omit it. |
| 6 | `numShadedRegions` | Number of standard 12 h night-shading regions. |
| 7 | `plotOriginalTimesamples` | `1` to plot individual raw traces; `0` to hide them. |
| 8 | `plotTimesamplesMean` | `1` to plot the group mean trace; `0` to hide it. |
| 9 | `plotTimesamplesStDeviation` | `1` to show the ±SD band; `0` to hide it. |
| 10 | `plotIndividualFittedCurves` | `1` to overlay individual fitted curves; `0` to hide them. |
| 11 | `plotAverageFittedCurve` | `1` to overlay the mean fitted curve and associated segment summary; `0` to hide it. |
| 12 | `figureHeight` | Main figure height in inches. |
| 13 | `colorChoice` | `shadesOfRed`, `shadesOfGreen`, `shadesOfBlue`, `shadesOfYellow`, `shadesOfOrange`, `shadesOfBrown`, `shadesOfGrey`, `shadesOfPurple`, `multicolor`, or `blackAndWhite`. |
| 14 (optional) | `phaseOffsetHours` | Reference time in hours for reported acrophase. If omitted or blank, it defaults to `0`. |

Example parameter file:

```text
Sample_00001,Sample_00002,Sample_00003
Exp01Temperature/33C/example_strain_33C
96 192; 192 288
1
0 0 0 0
8
1
1
1
0
1
12
shadesOfOrange
96
```

In this example, acrophase is reported relative to 96 h, rather than relative to the experimental time origin.

### Path convention

Line 2 of the parameter file controls the output prefix. A recommended organization is:

```text
repository-root/
├── Plot_rawLuminescence.m
├── processSamples.m
├── analyzeAndPlotSegments.m
├── fitCosineCurve.m
├── Plot_summary.m
├── run_luminescence_plots.m
├── Q10_config.m
├── Q10_plot.m
├── PRC_config.m
├── PRC_plot.m
├── _RunMatlabFronTerm_MyExperiment.sh
└── MyExperiment/
    ├── plate_A.csv
    ├── condition_A/
    │   └── strain_A.txt
    └── condition_B/
        └── strain_A.txt
```

Use forward slashes in paths for maximum portability. Keep parameter files and their derived outputs in the same condition-specific subdirectory whenever possible.

## Workflow

Run the pipeline in the following order.

### 1. Analyse individual groups

The supplied `.sh` scripts are batch launchers. Each command assigns `fileName` and `paramfileName` in MATLAB's base workspace, then executes `Plot_rawLuminescence.m` without opening the MATLAB desktop.

Run an existing launcher from the repository root:

```bash
bash _RunMatlabFronTerm_Exp01_Q10s.sh
```

or, if executable:

```bash
./_RunMatlabFronTerm_Exp01_Q10s.sh
```

A minimal launcher command has this form:

```bash
matlab -nodisplay -batch "assignin('base','fileName','MyExperiment/plate_A.csv'); assignin('base','paramfileName','MyExperiment/condition_A/strain_A.txt'); run('Plot_rawLuminescence.m')"
```

To create a launcher for a new experiment, copy an existing `.sh` file and change only the raw-CSV and parameter-file paths. One MATLAB invocation corresponds to one raw data file/parameter-file analysis group.

Alternatively, run a single group interactively from the repository root:

```matlab
fileName = 'MyExperiment/plate_A.csv';
paramfileName = 'MyExperiment/condition_A/strain_A.txt';
run('Plot_rawLuminescence.m')
```

### 2. Create summary luminescence figures

After the single-group analyses have generated `_MeanStd.csv` files, edit and run:

```matlab
run('Plot_summary.m')
```

`Plot_summary.m` is a configuration file. It specifies:

- the `_MeanStd.csv` files to combine;
- optional per-series colors as `'#RRGGBB'` hex values;
- overlay and/or tiled output PDF paths;
- which series appear in each tile;
- tile layout and ordering;
- standard night shading and custom shading intervals;
- y-axis scaling mode (`adaptive`, `fixed`, or `normalized`); and
- annotation and legend settings.

For annotation, `run_luminescence_plots.m` looks next to each input `_MeanStd.csv` for its matching `_individual_periods.csv` and `_individual_phases.csv` files. If either file is absent or has no valid numeric values, the corresponding annotation component is omitted.

### 3. Create period/Q10 figures

After individual group analyses have produced `_individual_periods.csv` files across temperatures, edit the strain registry and `plotDefs` blocks in:

```matlab
run('Q10_config.m')
```

For every temperature condition, `Q10_config.m` constructs a path to:

```text
<rootDir>/<temperature-folder>/<strain-key>_<temperature-folder>_individual_periods.csv
```

The Q10 workflow displays per-temperature boxplots and jittered individual values, fits a selected period-versus-temperature model, calculates regression-based and adjacent-pair \(Q_{10}\) values, and produces vector SVG plus 300 dpi PNG output.

The regression-based calculation is:

\[
Q_{10} = \left(\frac{\tau_{T_{\min}}}{\tau_{T_{\max}}}\right)^{10/(T_{\max}-T_{\min})},
\]

where \(\tau\) values are predictions from the selected polynomial model. Configurable model choices are `linear`, `quadratic`, `cubic`, and `linearsubset`. `blankTemps` preserves an x-axis position while omitting data from display and regression; `excludeTemps` can exclude specified temperatures from a `linearsubset` regression.

### 4. Create phase-response curves

After reference and pulsed-condition analyses have produced `_individual_phases.csv` files, edit the color map, strain-style registry, and `plotDefs` blocks in:

```matlab
run('PRC_config.m')
```

`PRC_config.m` uses expected folder and file naming patterns for reference and pulsed measurements. Verify that the generated paths match your dataset before running. The PRC engine reads individual phases and calculates:

\[
\Delta\Phi = \overline{\Phi}_{\mathrm{reference}} - \overline{\Phi}_{\mathrm{pulsed}}.
\]

Under this convention, positive values indicate a phase advance and negative values indicate a phase delay. Error bars are propagated as:

\[
\mathrm{SEM}_{\Delta\Phi} = \sqrt{\mathrm{SEM}_{\mathrm{reference}}^2 + \mathrm{SEM}_{\mathrm{pulsed}}^2}.
\]

The PRC workflow writes vector PDF figures.

## Analysis model

### Chirped damped cosine

Each selected sample is fitted independently within each requested segment using an FFT-initialized nonlinear least-squares model:

\[
y(t) = A e^{-\lambda u}\cos\left(\omega_0 u + \frac{1}{2}\kappa u^2 + \phi_0\right) + B,
\]

where \(u=t-t_0\), with \(t_0\) defined as the segment start. The model parameters are:

- \(A\): amplitude at the start of the segment;
- \(\lambda\): damping coefficient in h\(^{-1}\); positive values indicate declining amplitude and negative values indicate increasing amplitude;
- \(\omega_0\): angular frequency at the segment start;
- \(\kappa\): chirp rate in rad h\(^{-2}\); positive values correspond to frequency increase/period shortening and negative values to frequency decrease/period lengthening;
- \(\phi_0\): initial phase; and
- \(B\): fitted vertical offset (mesor).

The instantaneous period is:

\[
T(t) = \frac{2\pi}{\omega_0 + \kappa(t-t_0)}.
\]

The reported primary period for a segment is the instantaneous period at the segment midpoint. The code bounds periods to the 18–30 h fitting interval and applies multi-start optimization across three chirp-rate initializations (zero, positive, and negative) to reduce sensitivity to local minima. Acrophase is calculated numerically from the fitted curve on a dense time grid and is wrapped relative to `phaseOffsetHours`.

### Segment adjustment

Before fitting, the analysis may adjust a requested segment end time using landmarks detected in the group-average trace. The procedure estimates an FFT period, identifies the last peak or trough near the requested endpoint, and targets approximately 66.7% of the corresponding half-cycle before snapping to the nearest sampled time. The console reports whether each segment endpoint was adjusted. This behavior is intentional and should be considered when specifying fitting windows and comparing results across analyses.

## Output files

For a parameter-file output prefix `<prefix>`, the main analysis writes the following files alongside that prefix.

| Output | Contents |
|---|---|
| `<prefix>_MeanStd.csv` | Time, group mean, and group standard deviation calculated from selected, interpolated traces. |
| `<prefix>_plot.pdf` | Main luminescence figure exported as vector graphics. |
| `<prefix>_plot.fig` | MATLAB figure file for interactive editing. |
| `<prefix>_individual_periods.csv` | Individual fitted periods at each segment midpoint. |
| `<prefix>_individual_periods_mean.csv` | Per-sample time-averaged periods within each segment. |
| `<prefix>_individual_deltaT.csv` | Per-sample period change from the start to end of each segment. |
| `<prefix>_individual_phases.csv` | Individual fitted acrophases for each segment. |
| `<prefix>_individual_chirprates.csv` | Chirp-rate estimates \(\kappa\) for each sample and segment. |
| `<prefix>_individual_lambdas.csv` | Per-fit damping estimates \(\lambda\). |
| `<prefix>_individual_peak_times.csv` | First two raw peak times detected for each sample and segment. |
| `<prefix>_damping_coefficients.csv` | Cross-segment damping estimates for individual samples. |
| `<prefix>_phase_coherence.csv` | Phase-coherence summaries for each segment. |
| `<prefix>_phase_drift.csv` | Per-sample phase-drift slopes across segments. |

Generated results are reproducible products and can be excluded from version control if raw CSV files, parameter files, code, and launcher/configuration scripts are retained.

## Customization guide

### Add a condition or strain

1. Place the raw CSV within the relevant experiment directory.
2. Create one parameter file per biological group/strain/condition.
3. Add a command for that parameter file to an appropriate `.sh` launcher.
4. Run the launcher to generate the group-level outputs.
5. Add the relevant `_MeanStd.csv` paths to `Plot_summary.m` for summary luminescence figures.
6. Add temperature-specific `_individual_periods.csv` paths in `Q10_config.m` for Q10 analysis, if applicable.
7. Add reference/pulsed phase-file relationships in `PRC_config.m` for PRC analysis, if applicable.

### Change the fitting window

Modify line 3 of the relevant parameter file. For example:

```text
96 192; 192 288
```

specifies two independent fitting windows. The code requires an \(N\times2\) matrix after parsing; use semicolons between segments.

### Change phase reference

Set line 14 in the parameter file to the time to be treated as phase zero. For example, `0` reports phases relative to the experiment origin; `96` reports phases relative to 96 h. Omit the line or leave it blank to use the default of `0`.

### Change plot appearance

- Use parameter-file lines 7–13 to control individual traces, mean/SD, fitted curves, figure height, and color palette for raw-analysis plots.
- Use `Plot_summary.m` for colors, figure dimensions, overlays, tile order, shading, annotations, legends, and y-axis scaling in multi-condition summaries.
- Use `Q10_config.m` for strain labels, fit model, period limits, temperature colors, figure layout, and omitted/blank temperature positions.
- Use `PRC_config.m` for PRC colors, strain labels, marker/line styles, figure layout, tile definitions, and data-path conventions.

## Reproducibility and version control

Recommended practices:

- Commit the MATLAB source files, shell launchers, raw input data when sharing permissions allow, and every parameter/configuration file required to recreate a figure.
- Record the MATLAB release and toolbox versions used for each release or manuscript analysis.
- Keep a stable tag or archived commit associated with each submitted figure set.
- Do not overwrite original raw CSV files. Use parameter files and output prefixes to define distinct analyses.
- Treat output files as derived products. If they are not committed, document the exact commands used to recreate them.
- Review the terminal output after each batch run for missing files, invalid sample-column names, fit warnings, and segment-end adjustments.

A useful `.gitignore` starting point is:

```gitignore
# MATLAB autosave and temporary files
*.asv
*.autosave
*.m~

# Regenerable MATLAB figures and analysis products
*_plot.fig
*_plot.pdf
*_MeanStd.csv
*_individual_*.csv
*_damping_coefficients.csv
*_phase_coherence.csv
*_phase_drift.csv
*.png
*.svg

# OS/editor files
.DS_Store
Thumbs.db
```

Do **not** ignore parameter `.txt` files, raw input `.csv` files, `.m` source files, or `.sh` launchers if they are required to reproduce analyses. Adjust the ignore list if you intentionally want selected final manuscript figures or curated output tables tracked in Git.

## Troubleshooting

| Symptom | Likely cause and resolution |
|---|---|
| `Undefined function 'fit'` or `fittype` | Curve Fitting Toolbox is unavailable. Install/license it, then restart MATLAB. |
| `Undefined function 'findpeaks'` | Signal Processing Toolbox is unavailable. Install/license it. |
| `Undefined function 'tcdf'` or mixed-effects-model errors in Q10 analysis | Statistics and Machine Learning Toolbox is unavailable. Install/license it. |
| `matlab: command not found` | Add the MATLAB executable directory to your shell `PATH`, use MATLAB's full executable path, or run the workflow interactively. |
| `Could not open parameter file` | Run from the repository root or correct the relative parameter-file path in the shell launcher. |
| `Sample names not found in data` | Confirm that line 1 of the parameter file exactly matches the raw CSV column headers. |
| Q10/PRC input file not found | Check that output prefixes, folder names, strain keys, and temperature/ZT naming in the configuration file match the actual generated output names. |
| Empty/missing annotations in summary plots | Verify that matching `_individual_periods.csv` and `_individual_phases.csv` files were generated next to the referenced `_MeanStd.csv` file. |
| Fits return `NaN` | Check that the selected segment has at least five time points, adequate temporal sampling, and a signal with a plausible 18–30 h periodic component. Inspect raw data and fitting-window choices. |
| Batch run opens windows or fails under a display-less session | Use the provided `-nodisplay -batch` form and ensure the MATLAB release/environment supports headless graphics export. |

## Limitations and interpretation

- The pipeline estimates a model-based period and acrophase; results depend on sampling density, fit window, detrending choice, biological replicate quality, and the suitability of the chirped damped-cosine model.
- Missing observations are interpolated/extrapolated before group summary and fitting. Extensive missing data can therefore materially affect estimates and should be reviewed rather than treated as neutral.
- The period search is constrained to 18–30 h. Oscillations outside this range are not appropriate for the default fit settings without modifying `Plot_rawLuminescence.m`/`analyzeAndPlotSegments.m`.
- Q10 and PRC calculations are descriptive analysis outputs. They do not replace experimental-design decisions, independent replicate structure, or statistical inference tailored to the underlying biological question.
- Generated filenames are defined by parameter-file and configuration conventions. Renaming outputs manually can break downstream summary, Q10, or PRC path resolution.

## License and copyright

Copyright (c) 2026 **[Arnaud Taton / Susan S. Golden Laboratory / University of California San Diego]**.

Unless a different license file is included in this repository, all rights are reserved. You may inspect and use this code for internal research purposes, but redistribution, modification, publication of derivative software, or commercial use is not granted by this notice alone.

Before making a public GitHub repository, choose and add a dedicated `LICENSE` file. A common permissive option is the MIT License; it permits reuse, modification, and redistribution provided that the copyright and license notice are retained. If you select MIT, replace the copyright section above with the license notice and include the complete MIT text in `LICENSE`.

Suggested source-file header for project-created `.m` and `.sh` files:

```text
Copyright (c) 2026 [Copyright holder / laboratory / institution name]

This file is part of the MATLAB Circadian Bioluminescence Analysis Pipeline.
See the LICENSE file in the repository root for license terms.
```

### Third-party and platform notices

- MATLAB and associated toolbox names are trademarks of The MathWorks, Inc. This repository is not affiliated with, endorsed by, or sponsored by The MathWorks, Inc.
- Input datasets may be subject to separate ownership, consent, collaboration, publication, or repository-sharing restrictions. Confirm that raw data and metadata may legally and ethically be distributed before placing them in a public repository.
- If any code or data in the repository originated from another person, laboratory, or project, preserve its original attribution and license terms and obtain permission before relicensing or redistributing it.

## Citation

If you use this pipeline in a manuscript, cite the associated publication when available and identify the repository version or commit used for analysis. Until a formal archival release is created, a suggested acknowledgment is:

```text
Circadian bioluminescence data were analyzed using the MATLAB Circadian
Bioluminescence Analysis Pipeline (MATLAB R2025a; repository version/commit:
[insert identifier]).
```
## Acknowledgments
The code for this pipeline was developed with assistance from Perplexity AI, and this README was drafted by Perplexity AI and reviewed by the repository maintainer.

## Contact

For questions, bug reports, or requests, open a GitHub issue in the repository or contact **[Arnaud Taton and ataton@ucsd.edu]**.
