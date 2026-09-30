function params = processSamples(inputFileName, samplesNames)
% processSamples  Read an analysis parameter file and return all settings
%                 as a struct.
%
%   params = processSamples(inputFileName, samplesNames)
%
%   Inputs
%   ------
%   inputFileName - Path to the plain-text parameter file (see format below).
%   samplesNames  - Cell array of column names present in the data file
%                   (used to validate the requested samplesList).
%
%   Output
%   ------
%   params  - Struct with fields:
%               samplesList, sampleName, OutFileName, FigureName,
%               segments, detrenddata, phaseOffsetHours, resettingHours,
%               numShadedRegions, plotOriginalTimesamples,
%               plotTimesamplesMean, plotTimesamplesStDeviation,
%               plotIndividualFittedCurves, plotAverageFittedCurve,
%               figureHeight, colorChoice
%
%   Parameter file format -- UNCHANGED from original (13 lines)
%   -----------------------------------------------------------
%   Line  1 : Sample column names, comma-separated
%               e.g.  Sample_00001, Sample_00002, Sample_00003
%   Line  2 : Output name/prefix (used for CSV and figure file names)
%               e.g.  MyExperiment_WT
%   Line  3 : Segments to fit, semicolon-separated [start end] pairs
%               e.g.  96 192; 192 288
%   Line  4 : Detrend data before fitting  (1 = yes, 0 = no)
%   Line  5 : Dark-pulse patch coordinates [x1 x2 x2 x1]
%               e.g.  120 124 124 120   (or  0 0 0 0  to skip)
%   Line  6 : Number of shaded night regions
%   Line  7 : Plot original time samples  (1/0)
%   Line  8 : Plot mean time series       (1/0)
%   Line  9 : Plot SD shaded region       (1/0)
%   Line 10 : Plot individual fitted curves (1/0)
%   Line 11 : Plot average fitted curve    (1/0)
%   Line 12 : Figure height (inches)
%   Line 13 : Color palette name
%               shadesOfRed | shadesOfGreen | shadesOfBlue | shadesOfYellow |
%               shadesOfOrange | shadesOfBrown | shadesOfGrey | shadesOfPurple |
%               multicolor | blackAndWhite
%
%   Line 14 : [OPTIONAL] Phase offset reference time in hours.
%               Acrophase = mod(t_peak - phaseOffsetHours, period).
%               Defaults to 0 (ZT0/CT0) when line is absent or blank.
%               e.g.  96  means acrophase is relative to t = 96 h.
%
%   All existing 13-line parameter files work without modification.

% -------------------------------------------------------------------------
% Validate inputs
% -------------------------------------------------------------------------
if ~iscell(samplesNames)
    error('processSamples: samplesNames must be a cell array of strings.');
end

fileID = fopen(inputFileName, 'r');
if fileID == -1
    error('processSamples: Could not open parameter file: %s', inputFileName);
end

try

    % --- Line 1: sample names to analyse ---
    samplesList  = readStringList(fileID, 'samplesList');

    invalidNames = setdiff(samplesList, samplesNames);
    if ~isempty(invalidNames)
        error('processSamples: Sample names not found in data: %s', ...
              strjoin(invalidNames, ', '));
    end

    % --- Line 2: output prefix ---
    sampleName  = readString(fileID, 'sampleName');
    OutFileName = [sampleName, '_MeanStd.csv'];
    FigureName  = [sampleName, '_plot.pdf'];

    % --- Line 3: segments ---
    segLine  = readString(fileID, 'segments');
    segments = str2num(segLine); %#ok<ST2NM>
    if isempty(segments) || size(segments, 2) ~= 2
        error('processSamples: segments must be an Nx2 matrix, e.g.  96 192; 192 288');
    end

    % --- Line 4: detrend flag ---
    detrenddata = logical(readScalar(fileID, 'detrenddata', 0));

    % --- Line 5: resettingHours ---
    resetLine      = readString(fileID, 'resettingHours');
    resettingHours = str2num(resetLine); %#ok<ST2NM>
    if isempty(resettingHours)
        resettingHours = [0 0 0 0];
    end

    % --- Line 6: number of shaded night regions ---
    numShadedRegions = round(readScalar(fileID, 'numShadedRegions', 1));

    % --- Lines 7-11: plot flags ---
    plotOriginalTimesamples    = logical(readScalar(fileID, 'plotOriginalTimesamples',    1));
    plotTimesamplesMean        = logical(readScalar(fileID, 'plotTimesamplesMean',        0));
    plotTimesamplesStDeviation = logical(readScalar(fileID, 'plotTimesamplesStDeviation', 1));
    plotIndividualFittedCurves = logical(readScalar(fileID, 'plotIndividualFittedCurves', 0));
    plotAverageFittedCurve     = logical(readScalar(fileID, 'plotAverageFittedCurve',     1));

    % --- Line 12: figure height ---
    figureHeight = readScalar(fileID, 'figureHeight', 15);

    % --- Line 13: colour palette ---
    colorChoice   = readString(fileID, 'colorChoice');
    validPalettes = {'shadesOfRed', 'shadesOfGreen',  'shadesOfBlue',   'shadesOfYellow', ...
                     'shadesOfOrange', 'shadesOfBrown', 'shadesOfGrey', 'shadesOfPurple', ...
                     'multicolor', 'blackAndWhite'};
    if isempty(colorChoice) || ~any(strcmp(colorChoice, validPalettes))
        if ~isempty(colorChoice)
            warning('processSamples: Unknown colorChoice "%s". Using blackAndWhite.', colorChoice);
        end
        colorChoice = 'blackAndWhite';
    end

    % --- Line 14 (OPTIONAL): phaseOffsetHours ---
    % If the file has only 13 lines this silently defaults to 0.
    rawLine = fgetl(fileID);
    if isequal(rawLine, -1) || isempty(strtrim(rawLine))
        phaseOffsetHours = 0;
    else
        phaseOffsetHours = str2double(strtrim(rawLine));
        if isnan(phaseOffsetHours)
            warning('processSamples: Could not parse line 14 as a number. Using 0.');
            phaseOffsetHours = 0;
        end
    end

catch ME
    fclose(fileID);
    rethrow(ME);
end

fclose(fileID);

% -------------------------------------------------------------------------
% Pack into output struct
% -------------------------------------------------------------------------
params.samplesList                = samplesList;
params.sampleName                 = sampleName;
params.OutFileName                = OutFileName;
params.FigureName                 = FigureName;
params.segments                   = segments;
params.detrenddata                = detrenddata;
params.phaseOffsetHours           = phaseOffsetHours;
params.resettingHours             = resettingHours;
params.numShadedRegions           = numShadedRegions;
params.plotOriginalTimesamples    = plotOriginalTimesamples;
params.plotTimesamplesMean        = plotTimesamplesMean;
params.plotTimesamplesStDeviation = plotTimesamplesStDeviation;
params.plotIndividualFittedCurves = plotIndividualFittedCurves;
params.plotAverageFittedCurve     = plotAverageFittedCurve;
params.figureHeight               = figureHeight;
params.colorChoice                = colorChoice;

fprintf('\n=== Parameters loaded from: %s ===\n', inputFileName);
fprintf('  Samples       : %s\n',     strjoin(params.samplesList, ', '));
fprintf('  Output prefix : %s\n',     params.sampleName);
fprintf('  Segments      : ');
for k = 1:size(params.segments, 1)
    fprintf('%d-%d h  ', params.segments(k,1), params.segments(k,2));
end
fprintf('\n');
fprintf('  Detrend       : %d\n',     params.detrenddata);
fprintf('  Phase ref     : %.1f h\n', params.phaseOffsetHours);
fprintf('  Color palette : %s\n',     params.colorChoice);
fprintf('  Figure height : %.1f in\n',params.figureHeight);

end % processSamples

% =========================================================================
% Private helpers
% =========================================================================

function line = readRawLine(fileID, fieldName)
    line = fgetl(fileID);
    if isequal(line, -1)
        error('processSamples: End of file reached while reading field "%s". Check the parameter file has 13 lines.', fieldName);
    end
    line = strtrim(line);
end

function str = readString(fileID, fieldName)
    str = readRawLine(fileID, fieldName);
end

function cellArr = readStringList(fileID, fieldName)
    line    = readRawLine(fileID, fieldName);
    cellArr = strtrim(strsplit(line, ','));
    cellArr = cellArr(~cellfun(@isempty, cellArr));
end

function val = readScalar(fileID, fieldName, defaultVal)
    line = readRawLine(fileID, fieldName);
    val  = str2double(line);
    if isnan(val)
        warning('processSamples: Cannot parse "%s" for field "%s". Using default: %g.', line, fieldName, defaultVal);
        val = defaultVal;
    end
end
