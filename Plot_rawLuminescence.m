% =========================================================================
% Main script -- Circadian luminescence analysis
%
% Dependencies (must be on MATLAB path):
%   fitCosineCurve.m          -- FFT-NLLS cosine fitting
%   analyzeAndPlotSegments.m  -- per-segment analysis and plotting
%   processSamples.m          -- reads parameter file; returns a params struct
%
% Parameter file format (13 lines, optional 14th -- see processSamples.m):
%   Line  1 : Sample column names, comma-separated
%   Line  2 : Output name/prefix
%   Line  3 : Segments  e.g.  96 192; 192 288
%   Line  4 : Detrend  (1/0)
%   Line  5 : resettingHours  e.g.  120 124 124 120  (or  0 0 0 0)
%   Line  6 : numShadedRegions
%   Line  7 : plotOriginalTimesamples  (1/0)
%   Line  8 : plotTimesamplesMean      (1/0)
%   Line  9 : plotTimesamplesStDeviation (1/0)
%   Line 10 : plotIndividualFittedCurves (1/0)
%   Line 11 : plotAverageFittedCurve     (1/0)
%   Line 12 : figureHeight (inches)
%   Line 13 : colorChoice  (e.g. shadesOfOrange, multicolor)
%   Line 14 : phaseOffsetHours  [OPTIONAL -- defaults to 0 when absent]
%             Set to the reference time (h) for acrophase reporting.
%             e.g. 0 = ZT0/CT0;  96 = relative to t=96 h
% =========================================================================

% -------------------------------------------------------------------------
% Load data
% -------------------------------------------------------------------------
% filePrompt    = 'Enter the name of the CSV file with your data: ';
% fileName      = input(filePrompt, 's');
% paramfileName = input('Enter the name of the parameter file: ', 's');

data         = readtable(fileName);
time         = data.Time;
samplesNames = data.Properties.VariableNames(2:end);
fprintf('Available sample names: %s\n', strjoin(samplesNames, ', '));

% -------------------------------------------------------------------------
% Load parameters from file
% -------------------------------------------------------------------------
params = processSamples(paramfileName, samplesNames);

samplesList                = params.samplesList;
sampleName                 = params.sampleName;
segments                   = params.segments;
detrenddata                = params.detrenddata;
phaseOffsetHours           = params.phaseOffsetHours;
plotOriginalTimesamples    = params.plotOriginalTimesamples;
plotTimesamplesMean        = params.plotTimesamplesMean;
plotTimesamplesStDeviation = params.plotTimesamplesStDeviation;
plotIndividualFittedCurves = params.plotIndividualFittedCurves;
plotAverageFittedCurve     = params.plotAverageFittedCurve;
figureHeight               = params.figureHeight;
colorChoice                = params.colorChoice;
numShadedRegions           = params.numShadedRegions;
resettingHours             = params.resettingHours;
OutFileName                = params.OutFileName;
FigureName                 = params.FigureName;

% -------------------------------------------------------------------------
% Interpolate / fill missing values, then optionally detrend
% -------------------------------------------------------------------------
dataInterp = sortrows(data, 'Time');
xtime      = dataInterp.Time;

for col = 2:width(dataInterp)
    y         = dataInterp{:, col};
    validMask = ~isnan(y);
    if sum(validMask) >= 2
        dataInterp{:, col} = interp1(xtime(validMask), y(validMask), ...
                                      xtime, 'linear', 'extrap');
    end
end

% Remove rows where any sample column is still NaN
dataInterp = dataInterp(all(~isnan(dataInterp{:, 2:end}), 2), :);

for col = 2:width(dataInterp)
    dataInterp{:, col} = round(dataInterp{:, col});
    if detrenddata
        dataInterp{:, col} = detrend(dataInterp{:, col}, 'linear');
    end
end

dataInterp = sortrows(dataInterp, 'Time');

% -------------------------------------------------------------------------
% Colour palette
% -------------------------------------------------------------------------
numShades = max(40, length(samplesList));

shadesOfRed     = [ones(numShades,1), linspace(0.1,0.7,numShades)', linspace(0.1,0.7,numShades)'];
shadesOfGreen   = [linspace(0.2,0.5,numShades)', linspace(0.5,0.9,numShades)', linspace(0.1,0.4,numShades)'];
shadesOfBlue    = [zeros(numShades,1), linspace(0.3,1,numShades)', ones(numShades,1)];
shadesOfYellow  = [ones(numShades,1), linspace(0.8,1,numShades)', linspace(0,0.4,numShades)'];
shadesOfOrange  = [ones(numShades,1), linspace(0.3,1,numShades)', linspace(0.1,0.4,numShades)'];
shadesOfBrown   = [linspace(0.5,0.8,numShades)', linspace(0.3,0.5,numShades)', linspace(0.1,0.2,numShades)'];
shadesOfGrey    = repmat(linspace(0.2,0.8,numShades)', 1, 3);
shadesOfPurple  = [linspace(0.4,1,numShades)', zeros(numShades,1), linspace(0.4,1,numShades)'];
multicolorPalette    = lines(numShades);
blackAndWhitePalette = repmat(linspace(0.2,0.8,numShades)', 1, 3);

switch colorChoice
    case 'shadesOfRed';    selectedPalette = shadesOfRed;
    case 'shadesOfGreen';  selectedPalette = shadesOfGreen;
    case 'shadesOfBlue';   selectedPalette = shadesOfBlue;
    case 'shadesOfYellow'; selectedPalette = shadesOfYellow;
    case 'shadesOfOrange'; selectedPalette = shadesOfOrange;
    case 'shadesOfBrown';  selectedPalette = shadesOfBrown;
    case 'shadesOfGrey';   selectedPalette = shadesOfGrey;
    case 'shadesOfPurple'; selectedPalette = shadesOfPurple;
    case 'multicolor';     selectedPalette = multicolorPalette;
    otherwise
        if ~isempty(colorChoice)
            warning('Unknown colorChoice "%s". Defaulting to black and white.', colorChoice);
        end
        selectedPalette = blackAndWhitePalette;
end

% Same palette for both raw traces and fitted curves (visual consistency)
colors = selectedPalette(1:length(samplesList), :);

% -------------------------------------------------------------------------
% Create figure and plot raw data
% -------------------------------------------------------------------------
hMainFig = figure('Units', 'inches', 'Position', [0.75, 2.5, 10, (figureHeight/2) + 5]);
hold on;

if plotOriginalTimesamples
    for j = 1:length(samplesList)
        samplesName  = samplesList{j};
        samples      = data.(samplesName);
        labelStr = strrep(samplesName, '_', ' ');
        validIdx     = ~isnan(samples);
        validTime    = time(validIdx);
        validSamples = samples(validIdx);
        if detrenddata
            validSamples = detrend(validSamples, 'linear');
        end
        plot(validTime, validSamples, ':', ...
             'Marker', '.', 'MarkerSize', 10, 'LineWidth', 0.5, ...
             'Color', colors(j,:), 'DisplayName', labelStr);
    end
end

% -------------------------------------------------------------------------
% Mean and standard deviation (on interpolated data)
% -------------------------------------------------------------------------
overallMean = mean(dataInterp{:, samplesList}, 2);
overallStd  = std( dataInterp{:, samplesList}, 0, 2);

fileID = fopen(OutFileName, 'w');
fprintf(fileID, 'Time,Mean,StandardDeviation\n');
for i = 1:height(dataInterp)
    fprintf(fileID, '%f,%f,%f\n', dataInterp.Time(i), overallMean(i), overallStd(i));
end
fclose(fileID);
fprintf('Mean and SD exported to %s\n', OutFileName);

if plotTimesamplesMean
    plot(dataInterp.Time, overallMean, '-', 'LineWidth', 1, ...
         'Color', colors(1,:), 'DisplayName', 'Overall Mean');
end

if plotTimesamplesStDeviation
    fill([dataInterp.Time; flipud(dataInterp.Time)], ...
         [overallMean + overallStd; flipud(overallMean - overallStd)], ...
         colors(1,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'SD');
end

% -------------------------------------------------------------------------
% Cosine fitting -- analyzeAndPlotSegments
% -------------------------------------------------------------------------
[allPeriods, allPhases, allAdjustedPhases, allFittedCurves, segmentResults, ...
 dampingCoeffs, meanDamping, stdDamping, ...
 phaseDriftSlopes, meanPhaseDrift, stdPhaseDrift] = ...
    analyzeAndPlotSegments(dataInterp, dataInterp.Time, segments, samplesList, ...
                           colors, phaseOffsetHours, sampleName, ...
                           plotIndividualFittedCurves, plotAverageFittedCurve);

% Restore main figure before formatting - the floating window (if created)
% becomes gcf inside analyzeAndPlotSegments; without this line all
% formatting, night shading, and exportgraphics target the wrong figure.
figure(hMainFig);

% -------------------------------------------------------------------------
% Figure formatting
% -------------------------------------------------------------------------
ax = gca;

ax.XTick                 = 0:24:max(time);
ax.XMinorTick            = 'on';
ax.XAxis.MinorTickValues = 12:24:max(time)-12;

yTickMin = round(min(overallMean - overallStd) - 500, -3);
yTickMax = round(max(overallMean + overallStd) + 500, -3);
ax.YTick                 = yTickMin:5000:yTickMax;
ax.YAxis.TickLabelFormat = '%.0f';
ax.YAxis.Exponent        = 3;

set(ax, 'Units',    'inches');
set(ax, 'Position', [0.75, 0.5, 8.5, (figureHeight/2) + 2]);

% Night shading
yLimits = ylim;
for i = 1:numShadedRegions
    startTime = 12 + (i-1)*24;
    endTime   = 24 + (i-1)*24;
    fill([startTime endTime endTime startTime], ...
         [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Optional dark/resetting pulse patch
if ~isequal(resettingHours, [0 0 0 0])
    fill(ax, resettingHours, [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

title(extractBaseName(sampleName), 'Interpreter', 'none');
xlabel('Time (hours)');
ylabel('Luminescence (cps)');

grid on;  box on;
ax.XGrid              = 'on';
ax.YGrid              = 'off';
ax.XMinorGrid         = 'on';
ax.YMinorGrid         = 'off';
ax.GridColor          = [0 0 0];
ax.GridLineStyle      = ':';                
ax.GridAlpha          = 0.5;
ax.MinorGridColor     = [0.5 0.5 0.5];
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha     = 0.2;
ax.FontSize           = 15;
set(gca, 'LineWidth', 2);

lgd             = legend('show');
lgd.NumColumns  = 4;
lgd.Location    = 'southoutside';
lgd.Box         = 'off';
lgd.FontSize    = 6;
lgd.Interpreter = 'tex';

hold off;

% -------------------------------------------------------------------------
% Export figure
% -------------------------------------------------------------------------
exportgraphics(hMainFig, FigureName, 'ContentType', 'vector');
savefig(hMainFig, strcat(sampleName, '_plot.fig'));
fprintf('Figure exported to %s\n', FigureName);

% =========================================================================
% Local helper
% =========================================================================
function baseName = extractBaseName(fullPath)
    [~, baseName, ~] = fileparts(fullPath);
end
