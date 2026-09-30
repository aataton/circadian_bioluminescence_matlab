%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRC_config — Experiment configuration for PRC_plot
%
% This is the ONLY file you need to edit.
% Run this script to generate all PRC figures.
%
% HOW TO ADD A NEW PLOT
% ---------------------
%   1. Copy an existing plotDefs(N) block and increment p.
%   2. Set outputName, title, nCols, tileTitles.
%   3. Build each tile with makeTile() — pass multiple {strain,temp} pairs
%      to overlay several strains on the same axes.
%
% HOW TO ADD A NEW STRAIN / COLOR
% --------------------------------
%   Add entries in the STRAIN STYLE REGISTRY (marker, lineStyle, displayName).
%   Add per-ZT colors in the COLOR MAP section using the key:
%     '<strain>_<temp>_ZT'      for the reference (unused in plot)
%     '<strain>_<temp>_ZT<zt>'  for each pulsed condition
%
% PHASE-SHIFT CONVENTION
% ----------------------
%   ΔΦ = mean(reference) − mean(pulsed)
%   Positive → phase advance;   Negative → phase delay.
%
% makeTile SYNTAX
% ---------------
%   makeTile(fileDir, strainList, ZTs, colorMap)
%
%   strainList — cell array of {strain_key, temp_string} pairs, e.g.:
%                  { {'kaiA','29'}, {'kaiA148L156A','29'} }
%   ZTs        — vector of pulsed ZT times, e.g. [5 8 11]
%   colorMap   — the colorMap struct defined below
%
%   File naming convention (Exp03Resetting folder):
%     reference : Exp03Resetting/ZT_<temp>/<strain>_<temp>_ZT_individual_phases.csv
%     pulsed    : Exp03Resetting/ZT<zt>_<temp>/<strain>_<temp>_ZT<zt>_individual_phases.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% -- SETTINGS -----------------------------------------------------------

fileDir = '';   % root folder; leave '' if running from Exp03Resetting parent dir

cfg.tileWidth      = 5.5;   % inches per tile
cfg.tileHeight     = 4.0;   % inches per tile
cfg.hGap           = 0.6;   % horizontal gap between tiles (inches)
cfg.vGap           = 0.8;   % vertical gap between tiles (inches)
cfg.marginL        = 0.9;   % left margin (inches)
cfg.marginR        = 0.3;   % right margin (inches)
cfg.marginT        = 0.6;   % top margin for overall title (inches)
cfg.marginB        = 0.8;   % bottom margin (inches)
cfg.tickFontSize   = 18;
cfg.legendFontSize = 15;
cfg.axLabelSize    = 20;

% --- Line / marker / error-bar appearance (edit here to adjust) ---------
cfg.lineWidth   = 2.5;  % connecting line between ZT points
cfg.markerSize  = 9;    % marker diameter
cfg.ebLineWidth = 2.0;  % error bar stem thickness
cfg.capSize     = 8;    % error bar cap width
% ------------------------------------------------------------------------

%% -- COLOR MAP ----------------------------------------------------------
% Key format:  '<strain>_<temp>_ZT'      → reference condition (not plotted)
%              '<strain>_<temp>_ZT<zt>'  → color for that pulsed ZT point
% Values are hex strings, e.g. '#2B1A6B'
% Add a new row for each new strain/temp/ZT combination.

colorMap = struct();

% colorMap.dkaiAWT_29_ZT      = '#2B1A6B';
colorMap.dkaiAWT_29_ZT5     = '#6B55CC';
colorMap.dkaiAWT_29_ZT8     = '#6B55CC';
colorMap.dkaiAWT_29_ZT11    = '#6B55CC';

% colorMap.dkaiA148L156A_29_ZT    = '#8C1A0A';
colorMap.dkaiA148L156A_29_ZT5   = '#D63B25';
colorMap.dkaiA148L156A_29_ZT8   = '#D63B25';
colorMap.dkaiA148L156A_29_ZT11  = '#D63B25';

% colorMap.dkaiAWT_34_ZT      = '#2B1A6B';
colorMap.dkaiAWT_34_ZT5     = '#6B55CC';
colorMap.dkaiAWT_34_ZT8     = '#6B55CC';
colorMap.dkaiAWT_34_ZT11    = '#6B55CC';

% colorMap.dkaiA148L156A_34_ZT    = '#8C1A0A';
colorMap.dkaiA148L156A_34_ZT5   = '#D63B25';
colorMap.dkaiA148L156A_34_ZT8   = '#D63B25';
colorMap.dkaiA148L156A_34_ZT11  = '#D63B25';


%% -- STRAIN STYLE REGISTRY ----------------------------------------------
% Defines marker shape and line style per strain.
% Colors come from colorMap above (per ZT), not from here.

strainStyles = struct();


strainStyles.dkaiA148L156A.displayName = '\DeltakaiA^{148L156A}';
strainStyles.dkaiA148L156A.marker      = '.';
strainStyles.dkaiA148L156A.lineStyle   = '-';

strainStyles.dkaiAWT.displayName      = '\DeltakaiA WT';
strainStyles.dkaiAWT.marker           = '.';
strainStyles.dkaiAWT.lineStyle        = '-';

% Add more strains here:
% strainStyles.myStrain.displayName = 'My Strain';
% strainStyles.myStrain.marker      = 'p';
% strainStyles.myStrain.lineStyle   = '-.';

%% -- PLOT DEFINITIONS ---------------------------------------------------

p = 0;

% =========================================================================
% Plot 1 — dkaiA and kaiA CRISPR strains, 2×2 tile layout
%   Tile 1: dkaiA148L156A + dkaiAWT at 29°C
%   Tile 2: dkaiA148L156A + dkaiAWT at 34°C
% =========================================================================
p = p+1;
plotDefs(p).outputName = 'Exp03Resetting/PRC_dkaiA.pdf';
plotDefs(p).title      = 'PRC — \DeltakaiA Strains (4 h dark pulse)';
plotDefs(p).nCols      = 2;
% plotDefs(p).tileTitles = { ...
%     '\DeltakaiA strains  29°C', ...
%     '\DeltakaiA strains  34°C' ...
%     };

plotDefs(p).tiles = { ...
    makeTile(fileDir, {{'dkaiAWT','29'}, {'dkaiA148L156A','29'}}, [5 8 11], colorMap), ...
    makeTile(fileDir, {{'dkaiAWT','34'}, {'dkaiA148L156A','34'}}, [5 8 11], colorMap) ...
};


% =========================================================================
% Add more plots here — copy the block above and increment p
% =========================================================================
% p = p+1;
% plotDefs(p).outputName = 'PRC_myComparison';
% plotDefs(p).title      = 'PRC — My Comparison';
% plotDefs(p).nCols      = 2;
% plotDefs(p).tileTitles = {'Tile A', 'Tile B'};
% plotDefs(p).tiles = { ...
%     makeTile(fileDir, {{'kaiA','25'}, {'kaiA148L156A','25'}}, [5 8 11 16], colorMap), ...
%     makeTile(fileDir, {{'kaiA','30'}, {'kaiA148L156A','30'}}, [5 8 11 16], colorMap)  ...
% };

%% -- RUN ----------------------------------------------------------------
PRC_plot(plotDefs, strainStyles, cfg);
fprintf('\n=== All plots complete ===\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- makeTile ------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conditions = makeTile(fileDir, strainList, ZTs, colorMap)
% makeTile   Build the condition struct array for one PRC tile.
%
%   Each condition gets a .color field resolved from colorMap using the
%   key '<strain>_<temp>_ZT<zt>'. Falls back to black if key is missing.

conditions = struct('strain',{}, 'temp',{}, 'zt_pulse',{}, ...
                    'ref_file',{}, 'pulse_file',{}, 'color',{});

for si = 1:numel(strainList)
    strain = strainList{si}{1};
    temp   = strainList{si}{2};

    refFile = fullfile(fileDir, 'Exp03Resetting', ['ZT_' temp], ...
        sprintf('%s_%s_ZT_individual_phases.csv', strain, temp));

    for i = 1:numel(ZTs)
        zt = ZTs(i);
        pulsedFolder = sprintf('ZT%d_%s', zt, temp);
        pulsedFile   = fullfile(fileDir, 'Exp03Resetting', pulsedFolder, ...
            sprintf('%s_%s_ZT%d_individual_phases.csv', strain, temp, zt));

        % Resolve color from colorMap
        colorKey = sprintf('%s_%s_ZT%d', strain, temp, zt);
        colorKey = matlab.lang.makeValidName(colorKey);
        if isfield(colorMap, colorKey)
            col = hex2rgb(colorMap.(colorKey));
        else
            warning('PRC_config: no color defined for key "%s", using black.', colorKey);
            col = [0 0 0];
        end

        c.strain     = strain;
        c.temp       = [temp '°C'];
        c.zt_pulse   = zt;
        c.ref_file   = refFile;
        c.pulse_file = pulsedFile;
        c.color      = col;
        conditions(end+1) = c; %#ok
    end
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
