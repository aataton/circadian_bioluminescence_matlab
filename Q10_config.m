%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q10_config -- Data and layout configuration for Q10_plot
%
% This is the ONLY file you need to edit.
% Run this script to generate all period/Q10 figures.
%
% ADDING A NEW TILE / STRAIN
%   Copy a makeTile() call, change strainKey, update file paths.
%
% ADDING A NEW FIGURE
%   Increment p, set outputName/title/nCols/tileTitles, add tiles.
%
% FILE PATH CONVENTION
%   <rootDir>/<tempFolder>/<strainKey>_<tempFolder>_individual_periods.csv
%   e.g. Exp01Temperature/25C/dkaiA_wt_25C_individual_periods.csv
%
% CSV FORMAT (2 columns, 1 header row)
%   Sample,<any_column_name>
%   <sample_id>,<period_value>
%
% MODEL CHOICES
%   'linear'       -- Period ~ Temperature          (all temps)
%   'quadratic'    -- Period ~ Temperature + Temp^2
%   'cubic'        -- Period ~ Temperature + Temp^2 + Temp^3
%   'linearsubset' -- linear but excluding excludeTemps
%
% Q10 FORMULA
%   Q10 = (tau_Tmin / tau_Tmax)^(10 / (Tmax - Tmin))
%   where tau values come from the polyfit regression model.
%   Uncertainty via delta-method. Pairwise Q10 also computed.
%
% COLORS
%   colors: [] = automatic jet gradient
%   colors: Nx3 RGB matrix = one row per temperature, sorted ascending
%           N must equal the total number of temps including blankTemps.
%   Use hex2rgb() helper at the bottom to convert hex strings to [r g b].
%
% BLANK TEMPS
%   blankTemps: [] = show all temps (default)
%   blankTemps: e.g. [25, 27] = show x-axis tick but no data box/scatter
%               These temps are also excluded from regression automatically.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% -- SETTINGS -----------------------------------------------------------

rootDir = 'Exp01Temperature';   % root data folder

cfg.tileWidth      = 5.5;
cfg.tileHeight     = 4.5;
cfg.hGap           = 0.7;
cfg.vGap           = 0.9;
cfg.marginL        = 1.0;
cfg.marginR        = 0.4;
cfg.marginT        = 0.7;
cfg.marginB        = 0.9;
cfg.tickFontSize   = 14;
cfg.axLabelSize    = 16;
cfg.annotFontSize  = 15;
cfg.boxLineWidth   = 2.0;
cfg.markerSize     = 30;
cfg.jitterWidth    = 0.28;

%% -- STRAIN STYLE REGISTRY ----------------------------------------------
% colors rows: one per temperature, ASCENDING order (25->27->30->33->35)
% N rows must equal total temps including blankTemps.

strainStyles = struct();

% dkaiA_wt
strainStyles.dkaiA_wt.displayName  = '\DeltakaiA WT';
strainStyles.dkaiA_wt.periodMin    = 22;
strainStyles.dkaiA_wt.periodMax    = 28;
strainStyles.dkaiA_wt.model        = 'linearsubset';
strainStyles.dkaiA_wt.excludeTemps = [];
strainStyles.dkaiA_wt.blankTemps   = [25];      % tick shown, no data
strainStyles.dkaiA_wt.colors       = [ ...
    hex2rgb('#E898C0'); ...   % 25C  cornflower blue  (blank — tick only)
    hex2rgb('#E070A8'); ...   % 27C  sky blue
    hex2rgb('#E8728C'); ...   % 30C  sage green
    hex2rgb('#C44B9A'); ...   % 33C  coral-orange
    hex2rgb('#7B2D8B')  ];    % 35C  crimson

% dkaiA_148L156A
strainStyles.dkaiA_148L156A.displayName  = '\DeltakaiA^{148L156A}';
strainStyles.dkaiA_148L156A.periodMin    = 20;
strainStyles.dkaiA_148L156A.periodMax    = 30;
strainStyles.dkaiA_148L156A.model        = 'linearsubset';
strainStyles.dkaiA_148L156A.excludeTemps = [];
strainStyles.dkaiA_148L156A.blankTemps   = [25];   % tick shown, no data
strainStyles.dkaiA_148L156A.colors       = [ ...
    hex2rgb('#2166AC'); ...   % 25C  steel blue       (blank — tick only)
    hex2rgb('#4DAC26'); ...   % 27C  teal-green
    hex2rgb('#F7D027'); ...   % 30C  golden yellow
    hex2rgb('#E8751A'); ...   % 33C  burnt orange
    hex2rgb('#A50026')  ];    % 35C  deep red

% kaiA
strainStyles.kaiA.displayName  = 'kaiA';
strainStyles.kaiA.periodMin    = 22;
strainStyles.kaiA.periodMax    = 28;
strainStyles.kaiA.model        = 'linear';
strainStyles.kaiA.excludeTemps = [];
strainStyles.kaiA.blankTemps   = [];
strainStyles.kaiA.colors       = [ ...
    hex2rgb('#E898C0'); ...   % 25C
    hex2rgb('#E070A8'); ...   % 27C
    hex2rgb('#E8728C'); ...   % 30C
    hex2rgb('#C44B9A'); ...   % 33C
    hex2rgb('#7B2D8B')  ];    % 35C

% kaiA148L156A
strainStyles.kaiA148L156A.displayName  = 'kaiA^{148L156A}';
strainStyles.kaiA148L156A.periodMin    = 20;
strainStyles.kaiA148L156A.periodMax    = 30;
strainStyles.kaiA148L156A.model        = 'linearsubset';
strainStyles.kaiA148L156A.excludeTemps = [25 27];
strainStyles.kaiA148L156A.blankTemps   = [25 27];  % ticks shown, no data
strainStyles.kaiA148L156A.colors       = [ ...
    hex2rgb('#2166AC'); ...   % 25C  steel blue       (blank — tick only)
    hex2rgb('#4DAC26'); ...   % 27C  teal-green       (blank — tick only)
    hex2rgb('#F7D027'); ...   % 30C  golden yellow
    hex2rgb('#E8751A'); ...   % 33C  burnt orange
    hex2rgb('#A50026')  ];    % 35C  deep red

% --- Add more strains here ---
% strainStyles.myStrain.displayName  = 'My Strain';
% strainStyles.myStrain.periodMin    = 20;
% strainStyles.myStrain.periodMax    = 30;
% strainStyles.myStrain.model        = 'linear';
% strainStyles.myStrain.excludeTemps = [];
% strainStyles.myStrain.blankTemps   = [];
% strainStyles.myStrain.colors       = [];   % [] = auto jet, or Nx3 matrix

%% -- PLOT DEFINITIONS ---------------------------------------------------

p = 0;

% =========================================================================
% Plot 1: all four strains, 2x2 layout
% =========================================================================
p = p+1;
plotDefs(p).outputName = 'Exp01Temperature/Periods_Q10_dkaiA_wt_148L156A';
plotDefs(p).title      = 'Period vs Temperature & Q_{10}';
plotDefs(p).nCols      = 2;
plotDefs(p).tileTitles = { ...
    'KaiA WT', ...
    'dN-KaiA L156A' ...
    };

plotDefs(p).tiles = { ...
    makeTile(rootDir, 'dkaiA_wt', ...
        {25,'25C'},{27,'27C'},{30,'30C'},{33,'33C'},{35,'35C'}), ...
    makeTile(rootDir, 'dkaiA_148L156A', ...
        {25,'25C'},{27,'27C'},{30,'30C'},{33,'33C'},{35,'35C'}), ...
};


% =========================================================================
% Plot 1: all four strains, 2x2 layout
% =========================================================================
p = p+1;
plotDefs(p).outputName = 'Exp01Temperature/Periods_Q10_CRISPRkaiA_wt_148L156A';
plotDefs(p).title      = 'Period vs Temperature & Q_{10}';
plotDefs(p).nCols      = 2;
plotDefs(p).tileTitles = { ...
    'KaiA WT*', ...
    'dN-KaiA L156A*'...
    };

plotDefs(p).tiles = { ...
    makeTile(rootDir, 'kaiA', ...
        {25,'25C'},{27,'27C'},{30,'30C'},{33,'33C'},{35,'35C'}), ...
    makeTile(rootDir, 'kaiA148L156A', ...
        {25,'25C'},{27,'27C'},{30,'30C'},{33,'33C'},{35,'35C'}) ...
};


% =========================================================================
% Add more plots here -- copy the block above and increment p
% =========================================================================
% p = p+1;
% plotDefs(p).outputName = 'Q10_dkaiA_only';
% plotDefs(p).title      = 'Period vs Temperature -- dkaiA strains';
% plotDefs(p).nCols      = 2;
% plotDefs(p).tileTitles = {'\DeltakaiA WT', '\DeltakaiA^{148L156A}'};
% plotDefs(p).tiles = { ...
%     makeTile(rootDir, 'dkaiA_wt', ...
%         {25,'25C'},{30,'30C'},{35,'35C'}), ...
%     makeTile(rootDir, 'dkaiA_148L156A', ...
%         {25,'25C'},{30,'30C'},{35,'35C'}) ...
% };

%% -- RUN ----------------------------------------------------------------
Q10_plot(plotDefs, strainStyles, cfg);
fprintf('\n=== All Q10 plots complete ===\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- makeTile ------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tile = makeTile(rootDir, strainKey, varargin)
% makeTile  Build one tile definition.
%
%   tile = makeTile(rootDir, strainKey, {temp1,folder1}, {temp2,folder2}, ...)
%
%   Each varargin entry: {numericTemp, folderName}
%   Builds path: <rootDir>/<folderName>/<strainKey>_<folderName>_individual_periods.csv

tile.strainKey = strainKey;
tile.temps     = zeros(1, numel(varargin));
tile.files     = cell(1,  numel(varargin));

for k = 1:numel(varargin)
    tempNum       = varargin{k}{1};
    tempFolder    = varargin{k}{2};
    tile.temps(k) = tempNum;
    tile.files{k} = fullfile(rootDir, tempFolder, ...
        sprintf('%s_%s_individual_periods.csv', strainKey, tempFolder));
end

end % makeTile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- hex2rgb -------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = hex2rgb(hex)
% Convert '#RRGGBB' hex string to [r g b] in 0-1 range.
hex = strtrim(hex);
if hex(1) == '#', hex = hex(2:end); end
rgb = double([hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))]) / 255;
end
