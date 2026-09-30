%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENTS CONFIGURATION FILE

% HOW TO RUN ONE EXPERIMENT:
% 1. Run the SHARED SETTINGS section once (Ctrl+Enter on that block)
% 2. Click inside any experiment block and press Ctrl+Enter

% HOW TO RUN ALL EXPERIMENTS: press F5

% ── Per-experiment fields ─────────────────────────────────────────────
% fileNames cell array of CSV paths (or {path, '#hexcolor'} pairs)
%           If a hex color is supplied next to a filename it overrides
%           the palette color for that series.  Examples:
%             {'myfile_MeanStd.csv'}           → use palette color
%             {'myfile_MeanStd.csv', '#E86A9D'} → use this pink
% tileSeriesIdx {[1 2 3], [4 5 6], ...} — series per tile
% tileColumns number of tile columns (1 = stacked, 2+ = grid)
% overlayPDF output PDF for overlay ('' = skip)
% tilesPDF output PDF for tiles ('' = skip)
% numShadedRegions number of standard night light-grey bands
% customShading cell array, one cell per tile — Nx3 [tStart tEnd grey]
%               grey: 0=black … 1=white (0.35 = dark grey); [] = none

% ── tileOrder (set via cfg.tileOrder before each call) ───────────────
% cfg.tileOrder = 'rows' → fill left→right then down (default)
% cfg.tileOrder = 'cols' → fill top→down then right
% Set it right before the run_luminescence_plots() call so it applies
% only to that experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHARED DISPLAY SETTINGS (run once)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.figureWidth = 12;
cfg.overlayHeight = 6;
cfg.tileHeightFixed = 4;
cfg.marginBottom = 0.55;
cfg.marginTop = 0.15;
cfg.marginLeft = 0.85;
cfg.marginRight = 0.15;
cfg.tileGap = 1/2.54; % 1 cm between tiles
cfg.colGap = 0.10;
cfg.axLabelSize = 18;
cfg.xLabelOffset = 0.85;
cfg.yLabelOffset = 0.55;
cfg.legendFontSize = 10;
cfg.legendLocation = 'northeast';
cfg.legendMode = 'none';

% cfg.yMode options:
% 'adaptive' — each tile auto-scales to its own data range (default)
% 'fixed' — all tiles share the same Y limits (set cfg.yLimFixed)
% 'normalized' — each tile scaled to [-1, +1] based on its own
% global min/max across ALL series (mean ± std).
% Y-axis label changes to "Normalized Luminescence".
cfg.yMode = 'normalized';

cfg.yLimFixed = [-15000 25000];
cfg.yPadFrac = 0.10;
cfg.yTickSpacing = 5000;
cfg.overlayYTick = 2000;
cfg.tileOrder = 'rows'; % default — override per experiment below

% ── Annotation box settings ──────────────────────────────────────────
% cfg.showAnnotation — true (default): show period/phase box on every tile
%                      false : suppress all annotation boxes
% cfg.annotFontSize — font size for the annotation text (default: 9)

% The box displays the grand mean period and grand mean phase computed
% from the *individual* values of all series visible in each tile.
% Data are read from sibling CSV files in the same folder as each
% _MeanStd.csv:
%   _individual_periods_mean.csv → average period (h)
%   _individual_phases.csv       → average phase (h)
% If either file is absent the corresponding line is omitted silently.

cfg.showAnnotation = true;
cfg.annotFontSize = 8;
cfg.annotRight = true; % annotation box in upper-right corner
cfg.legendRowSpacing = 100; % extra vertical space between legend items
cfg.pal_seq = { ...
'#6A57B4','#E86A9D','#E5B832','#C95A46','#76B8E6','#6CA665'
};
% % purple, pink, yellow, red, blue, green  (palette; per-series colors set inline above)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exp01Temperature — dkaiA strain series, 33C, noDTR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileNames = { ...
    {'Exp01Temperature/33CnoDTR/dkaiA_wt_33C_MeanStd.csv',       '#6B55CC'}, ... % 1
    {'Exp01Temperature/33CnoDTR/dkaiA_ctrl_33C_MeanStd.csv',     '#B0B0B0'}, ... % 2
    {'Exp01Temperature/33CnoDTR/dkaiA_148_33C_MeanStd.csv',      '#E8B400'}, ... % 3
    {'Exp01Temperature/33CnoDTR/dkaiA_148L156A_33C_MeanStd.csv', '#D63B25'}, ... % 4
    {'Exp01Temperature/33CnoDTR/dkaiA_L156A_33C_MeanStd.csv',    '#28A046'}, ... % 5
    {'Exp01Temperature/33CnoDTR/dkaiA_181_33C_MeanStd.csv',      '#1F7DC4'}, ... % 6
};
tileSeriesIdx    = { 1:6 };
tileColumns      = 1;
numShadedRegions = 1;
customShading    = {};
overlayPDF       = 'Exp01Temperature/luminescence_33CnoDTR_dkaiAkaiA.pdf';
tilesPDF         = '';
cfg.tileOrder = 'cols'; % 1 column → order irrelevant
run_luminescence_plots(fileNames, overlayPDF, tilesPDF, ...
    tileSeriesIdx, tileColumns, numShadedRegions, customShading, cfg);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exp01Temperature — kaiA CRISPR strain series, 33C, noDTR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileNames = { ...
    {'Exp01Temperature/33CnoDTR/kaiA_33C_MeanStd.csv',       '#6B55CC'}, ... % 1 green
    {'Exp01Temperature/33CnoDTR/kaiA148_33C_MeanStd.csv',      '#E8B400'}, ... % 3 yellow
    {'Exp01Temperature/33CnoDTR/kaiA148L156A_33C_MeanStd.csv', '#D63B25'}, ... % 4 red
};
tileSeriesIdx    = { 1:3 };
tileColumns      = 1;
numShadedRegions = 1;
customShading    = {};
overlayPDF       = 'Exp01Temperature/luminescence_33CnoDTR_kaiACRISPR.pdf';
tilesPDF         = '';
cfg.tileOrder = 'cols'; % 1 column → order irrelevant
run_luminescence_plots(fileNames, overlayPDF, tilesPDF, ...
    tileSeriesIdx, tileColumns, numShadedRegions, customShading, cfg);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exp01Temperature — dkaiA str. WT and dNL156A, 35-25C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileNames = { ...
    {'Exp01Temperature/35C/dkaiA_wt_35C_MeanStd.csv',        '#7B2D8B'}, ... % 1  wt       35°C deep purple
    {'Exp01Temperature/35C/dkaiA_148L156A_35C_MeanStd.csv',  '#9B1C1C'}, ... % 2  148L156A 35°C dark red
    {'Exp01Temperature/33C/dkaiA_wt_33C_MeanStd.csv',        '#C44B9A'}, ... % 3  wt       33°C magenta-pink
    {'Exp01Temperature/33C/dkaiA_148L156A_33C_MeanStd.csv',  '#C95A20'}, ... % 4  148L156A 33°C burnt orange
    {'Exp01Temperature/30C/dkaiA_wt_30C_MeanStd.csv',        '#E8728C'}, ... % 5  wt       30°C warm pink
    {'Exp01Temperature/30C/dkaiA_148L156A_30C_MeanStd.csv',  '#E5B832'}, ... % 6  148L156A 30°C yellow
    {'Exp01Temperature/27C/dkaiA_wt_27C_MeanStd.csv',        '#E070A8'}, ... % 7  wt       27°C soft pink
    {'Exp01Temperature/27C/dkaiA_148L156A_27C_MeanStd.csv',  '#4A9B5E'}, ... % 8  148L156A 27°C green
    {'Exp01Temperature/25C/dkaiA_wt_25C_MeanStd.csv',      '#E898C0'}, ... % 9  wt       25°C light blush
    {'Exp01Temperature/25C/dkaiA_148L156A_25C_MeanStd.csv','#1A6BAF'}, ... % 10 148L156A 25°C blue
};
tileSeriesIdx    = {1:2, 3:4, 5:6, 7:8, 9:10};
tileColumns      = 1;
numShadedRegions = 1;
customShading    = {[], [], [], [], []};
overlayPDF       = '';
tilesPDF         = 'Exp01Temperature/luminescence_dkaiAkaiA_35-25C_tiles.svg';
cfg.tileOrder = 'cols'; % 1 column → order irrelevant
cfg.yMode = 'normalized';     % 'adaptive' | 'fixed' | 'normalized'
run_luminescence_plots(fileNames, overlayPDF, tilesPDF, ...
    tileSeriesIdx, tileColumns, numShadedRegions, customShading, cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exp01Temperature — kaiA str. WT and dNL156A, 35-25C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileNames = { ...
    {'Exp01Temperature/35C/kaiA_35C_MeanStd.csv',        '#7B2D8B'}, ... % 1  wt       35°C deep purple
    {'Exp01Temperature/35C/kaiA148L156A_35C_MeanStd.csv',  '#9B1C1C'}, ... % 2  148L156A 35°C dark red
    {'Exp01Temperature/33C/kaiA_33C_MeanStd.csv',        '#C44B9A'}, ... % 3  wt       33°C magenta-pink
    {'Exp01Temperature/33C/kaiA148L156A_33C_MeanStd.csv',  '#C95A20'}, ... % 4  148L156A 33°C burnt orange
    {'Exp01Temperature/30C/kaiA_30C_MeanStd.csv',        '#E8728C'}, ... % 5  wt       30°C warm pink
    {'Exp01Temperature/30C/kaiA148L156A_30C_MeanStd.csv',  '#E5B832'}, ... % 6  148L156A 30°C yellow
    {'Exp01Temperature/27C/kaiA_27C_MeanStd.csv',        '#E070A8'}, ... % 7  wt       27°C soft pink
    {'Exp01Temperature/27C/kaiA148L156A_27C_MeanStd.csv',  '#4A9B5E'}, ... % 8  148L156A 27°C green
    {'Exp01Temperature/25C/kaiA_25C_MeanStd.csv',      '#E898C0'}, ... % 9  wt       25°C light blush
    {'Exp01Temperature/25C/kaiA148L156A_25C_MeanStd.csv','#1A6BAF'}, ... % 10 148L156A 25°C blue
};
tileSeriesIdx    = {1:2, 3:4, 5:6, 7:8, 9:10};
tileColumns      = 1;
numShadedRegions = 1;
customShading    = {[], [], [], [], []};
overlayPDF       = '';
tilesPDF         = 'Exp01Temperature/luminescence_kaiACRISPR_35-25C_tiles.svg';
cfg.tileOrder = 'cols'; % 1 column → order irrelevant
cfg.yMode = 'normalized';     % 'adaptive' | 'fixed' | 'normalized'
run_luminescence_plots(fileNames, overlayPDF, tilesPDF, ...
    tileSeriesIdx, tileColumns, numShadedRegions, customShading, cfg);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exp02Anabaena4LDcycles — dkaiA str. incl. Anabaena, 30°C, EVO, noDTR
% 10 strains, 2 tiles in 1 column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileNames = { ...
    {'Exp02Anabaena4LDcycles/dkaiA_CTRL_MeanStd.csv',    '#B0B0B0'}, ... % 1   '#B0B0B0'
    {'Exp02Anabaena4LDcycles/dkaiA_WT_MeanStd.csv',      '#6B55CC'}, ... % 2   '#6B55CC'
    {'Exp02Anabaena4LDcycles/dkaiA_148_MeanStd.csv',     '#E8B400'}, ... % 3   '#E8B400'
    {'Exp02Anabaena4LDcycles/dkaiA_148L156A_MeanStd.csv','#D63B25'}, ... % 4   '#D63B25'
    {'Exp02Anabaena4LDcycles/dkaiA_L156A_MeanStd.csv',   '#28A046'}, ... % 5   '#28A046'
    {'Exp02Anabaena4LDcycles/dkaiA_181_MeanStd.csv',     '#1F7DC4'}, ... % 6   '#1F7DC4'
    {'Exp02Anabaena4LDcycles/dkaiA_CTRLgm_MeanStd.csv',  '#707070'}, ... % 7   '#707070'
    {'Exp02Anabaena4LDcycles/dkaiA_WTgm_MeanStd.csv',    '#7B52CC'}, ... % 8   '#7B52CC'
    {'Exp02Anabaena4LDcycles/dkaiA_AnPan_MeanStd.csv',   '#F5AE55'}, ... % 9   '#F5AE55'
    {'Exp02Anabaena4LDcycles/dkaiA_AnPse_MeanStd.csv',   '#D4306A'}, ... % 10  '#D4306A'
}; 
tileSeriesIdx    = { [1:6], [7:10] };
tileColumns      = 1;
numShadedRegions = 4;
customShading    = { [], [] };
overlayPDF       = 'Exp02Anabaena4LDcycles/luminescence_30C_inclAnabaena_4LDc_overlay.svg';
tilesPDF         = 'Exp02Anabaena4LDcycles/luminescence_30C_inclAnabaena_4LDc_tiles.svg';
cfg.tileOrder = 'rows'; % 1 column → order irrelevant
cfg.yMode = 'adaptive';     % 'adaptive' | 'fixed' | 'normalized'
run_luminescence_plots(fileNames, overlayPDF, tilesPDF, ...
    tileSeriesIdx, tileColumns, numShadedRegions, customShading, cfg);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exp03Resetting — ZT phi reset at 29/34°C deltaKaiA
% 20 conditions, 5 columns × 4 rows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileNames = { ...
    {'Exp03Resetting/ZT_29/dkaiAWT_29_ZT_MeanStd.csv',             '#2B1A6B'}, ... % 1  WT  ZT    deep indigo
    {'Exp03Resetting/ZT5_29/dkaiAWT_29_ZT5_MeanStd.csv',           '#6B55CC'}, ... % 2  WT  ZT5   indigo
    {'Exp03Resetting/ZT8_29/dkaiAWT_29_ZT8_MeanStd.csv',           '#A8389A'}, ... % 3  WT  ZT8   rich pink-purple
    {'Exp03Resetting/ZT11_29/dkaiAWT_29_ZT11_MeanStd.csv',         '#D4306A'}, ... % 4  WT  ZT11  crimson rose
    {'Exp03Resetting/ZT_29/dkaiA148L156A_29_ZT_MeanStd.csv',       '#8C1A0A'}, ... % 5  148 ZT    dark red
    {'Exp03Resetting/ZT5_29/dkaiA148L156A_29_ZT5_MeanStd.csv',     '#D63B25'}, ... % 6  148 ZT5   red
    {'Exp03Resetting/ZT8_29/dkaiA148L156A_29_ZT8_MeanStd.csv',     '#E8820C'}, ... % 7  148 ZT8   orange
    {'Exp03Resetting/ZT11_29/dkaiA148L156A_29_ZT11_MeanStd.csv',   '#E8B400'}, ... % 8  148 ZT11  amber-gold
    {'Exp03Resetting/ZT_34/dkaiAWT_34_ZT_MeanStd.csv',             '#2B1A6B'}, ... % 9  WT  ZT    deep indigo
    {'Exp03Resetting/ZT5_34/dkaiAWT_34_ZT5_MeanStd.csv',           '#6B55CC'}, ... % 10 WT  ZT5   indigo
    {'Exp03Resetting/ZT8_34/dkaiAWT_34_ZT8_MeanStd.csv',           '#A8389A'}, ... % 11 WT  ZT8   rich pink-purple
    {'Exp03Resetting/ZT11_34/dkaiAWT_34_ZT11_MeanStd.csv',         '#D4306A'}, ... % 12 WT  ZT11  crimson rose
    {'Exp03Resetting/ZT_34/dkaiA148L156A_34_ZT_MeanStd.csv',       '#8C1A0A'}, ... % 13 148 ZT    dark red
    {'Exp03Resetting/ZT5_34/dkaiA148L156A_34_ZT5_MeanStd.csv',     '#D63B25'}, ... % 14 148 ZT5   red
    {'Exp03Resetting/ZT8_34/dkaiA148L156A_34_ZT8_MeanStd.csv',     '#E8820C'}, ... % 15 148 ZT8   orange
    {'Exp03Resetting/ZT11_34/dkaiA148L156A_34_ZT11_MeanStd.csv',   '#E8B400'}, ... % 16 148 ZT11  amber-gold
};
tileSeriesIdx = { [1 2], [1 3], [1 4], [5 6], [5 7], [5 8], ...
                  [9 10], [9 11], [9 12], [13 14], [13 15], [13 16] };
tileColumns      = 2;
numShadedRegions = 3;
customShading    = { [77 81 0.025], [80 84 0.025], [83 87 0.025], ...
                     [77 81 0.025], [80 84 0.025], [83 87 0.025], ...
                     [77 81 0.025], [80 84 0.025], [83 87 0.025], ...
                     [77 81 0.025], [80 84 0.025], [83 87 0.025] };
overlayPDF       = '';
tilesPDF         = 'Exp03Resetting/luminescence_dkaiAkaiA_ZTreset_tiles.svg';
cfg.tileOrder = 'cols';
run_luminescence_plots(fileNames, overlayPDF, tilesPDF, ...
    tileSeriesIdx, tileColumns, numShadedRegions, customShading, cfg);






