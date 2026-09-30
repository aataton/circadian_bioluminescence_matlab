function run_luminescence_plots(fileNames, overlayPDF, tilesPDF, ...
    tileSeriesIdx, tileColumns, ...
    numShadedRegions, customShading, cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_luminescence_plots
% Called from experiments_config.m — do not edit unless changing behaviour.
%
% fileNames format (two styles are accepted):
%   Old style — plain strings:
%       fileNames = {'path/to/file_MeanStd.csv', ...}
%   New style — {path, hexcolor} pairs:
%       fileNames = {{'path/to/file_MeanStd.csv','#E86A9D'}, ...}
%   Styles can be mixed freely within the same cell array.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -- Safety-net defaults ------------------------------------------------
if ~isfield(cfg,'colGap'),           cfg.colGap = 0.10;        end
if ~isfield(cfg,'tileGap'),          cfg.tileGap = 0.08;       end
if ~isfield(cfg,'xLabelOffset'),     cfg.xLabelOffset = 0.30;  end
if ~isfield(cfg,'yLabelOffset'),     cfg.yLabelOffset = 0.55;  end
if ~isfield(cfg,'axLabelSize'),      cfg.axLabelSize = 18;     end
if ~isfield(cfg,'overlayYTick'),     cfg.overlayYTick = 2000;  end
if ~isfield(cfg,'tileOrder'),        cfg.tileOrder = 'rows';   end
if ~isfield(cfg,'yMode'),            cfg.yMode = 'adaptive';   end
if ~isfield(cfg,'showAnnotation'),   cfg.showAnnotation = true; end
if ~isfield(cfg,'annotFontSize'),    cfg.annotFontSize = 10;   end
if ~isfield(cfg,'annotRight'),       cfg.annotRight = true;   end
if ~isfield(cfg,'legendRowSpacing'), cfg.legendRowSpacing = 2; end

%% -- 1. PARSE fileNames — extract paths and per-series colors -----------
% Each element of fileNames may be:
%   (a) a plain string:          'path/file_MeanStd.csv'
%   (b) a 1-element cell:        {'path/file_MeanStd.csv'}
%   (c) a 2-element cell pair:   {'path/file_MeanStd.csv', '#RRGGBB'}

nSeries       = numel(fileNames);
filePaths     = cell(1, nSeries);   % resolved file paths
fileHexColors = cell(1, nSeries);   % hex strings or '' (use palette)

for i = 1:nSeries
    entry = fileNames{i};
    if ischar(entry) || isstring(entry)
        % plain string
        filePaths{i}     = char(entry);
        fileHexColors{i} = '';
    elseif iscell(entry)
        filePaths{i}     = char(entry{1});
        if numel(entry) >= 2 && ~isempty(entry{2})
            fileHexColors{i} = char(entry{2});
        else
            fileHexColors{i} = '';
        end
    else
        error('fileNames{%d} must be a string or a cell {path} / {path,hex}.', i);
    end
end

%% -- 1a. LOAD AND ASSEMBLE DATA -----------------------------------------

seriesNames = cell(1, nSeries);
for i = 1:nSeries
    [~, name, ~]  = fileparts(filePaths{i});
    seriesNames{i} = strrep(name, '_MeanStd', '');
end

combinedTable = table();
for i = 1:nSeries
    currentTable = readtable(filePaths{i});
    uniqueName   = seriesNames{i};
    varsToRename = currentTable.Properties.VariableNames( ...
        ~strcmp(currentTable.Properties.VariableNames, 'Time'));
    currentTable.Properties.VariableNames(varsToRename) = ...
        strcat(uniqueName, "_", varsToRename);
    if isempty(combinedTable)
        combinedTable = currentTable;
    else
        combinedTable = outerjoin(combinedTable, currentTable, ...
            'Keys','Time','MergeKeys',true,'Type','full');
    end
end

dataInterp = sortrows(combinedTable, 'Time');
xtime = dataInterp.Time;
for col = 2:width(dataInterp)
    col_data = dataInterp{:,col};
    valid    = ~isnan(col_data);
    if sum(valid) > 1
        dataInterp{:,col} = interp1(xtime(valid), col_data(valid), ...
            xtime, 'linear', 'extrap');
    end
end
dataInterp    = dataInterp(all(~isnan(dataInterp{:,2:end}), 2), :);
for col = 2:width(dataInterp)
    dataInterp{:,col} = round(dataInterp{:,col});
end
combinedTable = sortrows(dataInterp, 'Time');
time          = combinedTable.Time;

%% -- 1b. BUILD COLOUR ARRAY --------------------------------------------
% Start from the interpolated palette, then override any series that
% supplied an explicit hex colour in fileNames.

fullcolors = hex2rgb(cfg.pal_seq);
numColors  = length(seriesNames);
colors     = interp1(linspace(0,1,size(fullcolors,1)), fullcolors, ...
                     linspace(0,1,numColors));

for i = 1:nSeries
    if ~isempty(fileHexColors{i})
        colors(i,:) = hex2rgb({fileHexColors{i}});
    end
end

%% -- 1c. PRE-LOAD period / phase data for each series -------------------

seriesPeriodMean = nan(1, nSeries);
seriesPeriodStd  = nan(1, nSeries);
seriesPhaseMean  = nan(1, nSeries);
seriesPhaseStd   = nan(1, nSeries);

fprintf('\n--- Annotation data loading ---\n');
for i = 1:nSeries
    basePath   = strrep(filePaths{i}, '_MeanStd.csv', '');
    periodFile = [basePath '_individual_periods.csv'];
    phaseFile  = [basePath '_individual_phases.csv'];

    fprintf('  [%d] Trying period file: %s\n', i, periodFile);
    vals = readCol2_textscan(periodFile);
    if ~isempty(vals)
        seriesPeriodMean(i) = mean(vals);
        seriesPeriodStd(i)  = std(vals);
        fprintf('      -> n=%d, mean=%.4f, std=%.4f\n', ...
            numel(vals), seriesPeriodMean(i), seriesPeriodStd(i));
    else
        fprintf('      -> no valid numeric values found\n');
    end

    fprintf('  [%d] Trying phase file: %s\n', i, phaseFile);
    vals = readCol2_textscan(phaseFile);
    if ~isempty(vals)
        seriesPhaseMean(i) = mean(vals);
        seriesPhaseStd(i)  = std(vals);
        fprintf('      -> n=%d, mean=%.4f, std=%.4f\n', ...
            numel(vals), seriesPhaseMean(i), seriesPhaseStd(i));
    else
        fprintf('      -> no valid numeric values found\n');
    end
end
fprintf('--- End annotation loading ---\n\n');

%% -- 1d. Helper: annotation coverage check ------------------------------
function tf = annotationCoversAll(idx)
    tf = cfg.showAnnotation && all( ...
        isfinite(seriesPeriodMean(idx)) | isfinite(seriesPhaseMean(idx)) );
end

%% -- 2. OVERLAY FIGURE --------------------------------------------------

if ~isempty(overlayPDF)

    fig = figure('Units','inches','Position', ...
        [1.5, 1.5, cfg.figureWidth, cfg.overlayHeight]);
    hold on;

    for i = 1:length(seriesNames)
        meanData = combinedTable.([seriesNames{i} '_Mean']);
        stdData  = combinedTable.([seriesNames{i} '_StandardDeviation']);
        upper    = meanData + stdData;
        lower    = meanData - stdData;
        fill([time; flipud(time)], [upper; flipud(lower)], colors(i,:), ...
            'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
        plot(time, meanData, '-', 'Color',colors(i,:), 'LineWidth',2.5, ...
            'DisplayName', seriesNames{i});
    end

    yLimits = ylim;
    for i = 1:numShadedRegions
        fill([12+(i-1)*24, 24+(i-1)*24, 24+(i-1)*24, 12+(i-1)*24], ...
            [yLimits(1), yLimits(1), yLimits(2), yLimits(2)], ...
            'k','FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
    end

    ax = gca;
    grid on; box off;
    ax.XGrid = 'on'; ax.YGrid = 'off'; ax.XMinorGrid = 'off';
    ax.XTick = 0:24:max(time);
    xlim(ax, [0, max(time)]);
    ax.XMinorTick = 'off';
    ax.XAxis.MinorTickValues = 12:24:max(time)-12;
    ax.YTick = min(ylim):cfg.overlayYTick:max(ylim);
    ax.YAxis.TickLabelFormat = '%.0f';
    ax.YAxis.Exponent = 3;
    ax.GridColor = [0 0 0]; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
    ax.MinorGridColor = [0.5 0.5 0.5];
    ax.MinorGridLineStyle = '-'; ax.MinorGridAlpha = 0.2;
    ax.FontSize = 15;
    set(ax,'LineWidth',2);
    xlabel('Time (hours)','FontSize',20);
    ylabel('Luminescence (cps)','FontSize',20);

    allIdx = 1:length(seriesNames);
    if ~annotationCoversAll(allIdx)
        lgd = legend('show');
        lgd.NumColumns = 1; lgd.Location = 'best';
        lgd.FontSize = 13; lgd.Interpreter = 'none';
        lgd.ItemTokenSize(2) = cfg.legendRowSpacing;
    else
        legend(ax, 'off');
    end

    if cfg.showAnnotation
        addAnnotationBox(ax, allIdx, seriesNames, colors, ...
            seriesPeriodMean, seriesPeriodStd, ...
            seriesPhaseMean,  seriesPhaseStd, ...
            cfg.annotFontSize, cfg.annotRight);
    end

    hold(ax, 'off');

    fig.PaperUnits    = 'inches';
    fig.PaperSize     = [cfg.figureWidth, cfg.overlayHeight];
    fig.PaperPosition = [0, 0, cfg.figureWidth, cfg.overlayHeight];
    exportgraphics(fig, overlayPDF, 'ContentType','vector');

    fprintf('Overlay saved: %s\n', overlayPDF);
end

%% -- 3. TILE FIGURE -----------------------------------------------------

if ~isempty(tilesPDF)

    numTiles = numel(tileSeriesIdx);
    nCols    = max(1, round(tileColumns));
    nRows    = ceil(numTiles / nCols);

    if numel(customShading) ~= numTiles
        error(['customShading must have exactly %d cells (one per tile). ' ...
            'Got %d.'], numTiles, numel(customShading));
    end

    tileRow = zeros(numTiles,1);
    tileCol = zeros(numTiles,1);
    for tIdx = 1:numTiles
        if strcmp(cfg.tileOrder, 'cols')
            tileRow(tIdx) = mod(tIdx-1, nRows) + 1;
            tileCol(tIdx) = floor((tIdx-1) / nRows) + 1;
        else
            tileRow(tIdx) = floor((tIdx-1) / nCols) + 1;
            tileCol(tIdx) = mod(tIdx-1, nCols) + 1;
        end
    end

    yLims       = zeros(numTiles, 2);
    normParams  = zeros(numTiles, 2);

    for tIdx = 1:numTiles
        serIdx  = tileSeriesIdx{tIdx};
        allVals = [];
        for i = serIdx
            mn      = combinedTable.([seriesNames{i} '_Mean']);
            sd      = combinedTable.([seriesNames{i} '_StandardDeviation']);
            allVals = [allVals; mn+sd; mn-sd]; %#ok
        end
        if strcmp(cfg.yMode,'fixed')
            yLims(tIdx,:) = cfg.yLimFixed;
        elseif strcmp(cfg.yMode,'normalized')
            dMin = min(allVals); dMax = max(allVals);
            normParams(tIdx,:) = [dMin, dMax - dMin];
            yLims(tIdx,:)      = [-1, 1];
        else
            dMin = min(allVals); dMax = max(allVals);
            pad  = (dMax-dMin)*cfg.yPadFrac;
            yLims(tIdx,:) = [dMin-pad, dMax+pad];
        end
    end

    tileHeightsAll = zeros(numTiles,1);
    if strcmp(cfg.yMode,'adaptive')
        yRanges        = yLims(:,2) - yLims(:,1);
        tileHeightsAll = cfg.tileHeightFixed * yRanges / max(yRanges);
    else
        tileHeightsAll(:) = cfg.tileHeightFixed;
    end

    rowHeights = zeros(nRows,1);
    for row = 1:nRows
        tilesInRow      = find(tileRow == row);
        rowHeights(row) = max(tileHeightsAll(tilesInRow));
    end

    tW   = (cfg.figureWidth - cfg.marginLeft - cfg.marginRight ...
            - cfg.colGap*(nCols-1)) / nCols;
    figW = cfg.figureWidth;
    figH = sum(rowHeights) + cfg.marginBottom + cfg.marginTop ...
           + cfg.tileGap*(nRows-1);

    rowBottoms        = zeros(nRows,1);
    rowBottoms(nRows) = cfg.marginBottom;
    for row = nRows-1:-1:1
        rowBottoms(row) = rowBottoms(row+1) + rowHeights(row+1) + cfg.tileGap;
    end

    fig2 = figure('Units','inches','Position',[1.5,1.5,figW,figH]);

    axHandles    = gobjects(numTiles,1);
    sharedHandles = [];
    sharedLabels  = {};
    seenLabels    = {};
    tileAnnotFull = false(numTiles,1);

    for tIdx = 1:numTiles
        row = tileRow(tIdx);
        col = tileCol(tIdx);

        xLeft = cfg.marginLeft + (col-1)*(tW + cfg.colGap);
        xN    = xLeft / figW;
        yN    = rowBottoms(row) / figH;
        wN    = tW / figW;
        hN    = rowHeights(row) / figH;

        ax           = axes('Units','normalized','Position',[xN,yN,wN,hN]); %#ok
        axHandles(tIdx) = ax;
        hold(ax,'on');

        serIdx       = tileSeriesIdx{tIdx};
        tileHandles  = gobjects(1,numel(serIdx));
        tileLabels   = {};

        doNorm   = strcmp(cfg.yMode,'normalized');
        nrmMin   = normParams(tIdx,1);
        nrmScale = normParams(tIdx,2);
        if nrmScale == 0, nrmScale = 1; end

        for k = 1:numel(serIdx)
            i        = serIdx(k);
            meanData = combinedTable.([seriesNames{i} '_Mean']);
            stdData  = combinedTable.([seriesNames{i} '_StandardDeviation']);

            if doNorm
                meanData = 2*(meanData - nrmMin)/nrmScale - 1;
                stdData  = 2*stdData / nrmScale;
            end

            upper = meanData + stdData;
            lower = meanData - stdData;
            fill(ax,[time;flipud(time)],[upper;flipud(lower)],colors(i,:), ...
                'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
            h     = plot(ax,time,meanData,'-','Color',colors(i,:),'LineWidth',2.5);
            label = strrep(seriesNames{i},'_','\_');
            tileHandles(k) = h;
            tileLabels{k}  = label;
            if ~ismember(label,seenLabels)
                sharedHandles(end+1) = h; %#ok
                sharedLabels{end+1}  = label;
                seenLabels{end+1}    = label;
            end
        end

        ylim(ax, yLims(tIdx,:));
        tileShadingMatrix = customShading{tIdx};
        formatAxes(ax, time, numShadedRegions, tileShadingMatrix, cfg.yTickSpacing);

        if doNorm
            ax.YTick = -1:0.5:1;
            ax.YAxis.Exponent = 0;
            ax.YAxis.TickLabelFormat = '%.1f';
            ax.YTickLabel = arrayfun(@(v) sprintf('%.1f',v), ...
                ax.YTick, 'UniformOutput', false);
        else
            ax.YAxis.Exponent = 0;
            ax.YAxis.TickLabelFormat = '%.0f';
            ax.YTickLabel = arrayfun(@(v) sprintf('%.0f',v/1000), ...
                ax.YTick,'UniformOutput',false);
        end

        tileAnnotFull(tIdx) = annotationCoversAll(serIdx);
        if strcmp(cfg.legendMode,'individual')
            if ~tileAnnotFull(tIdx)
                lgd = legend(ax, tileHandles, tileLabels, ...
                    'Location', cfg.legendLocation, ...
                    'FontSize',  cfg.legendFontSize, ...
                    'Interpreter','tex');
                lgd.ItemTokenSize(2) = cfg.legendRowSpacing;
            else
                legend(ax, 'off');
            end
        end

        if cfg.showAnnotation
            addAnnotationBox(ax, serIdx, seriesNames, colors, ...
                seriesPeriodMean, seriesPeriodStd, ...
                seriesPhaseMean,  seriesPhaseStd, ...
                cfg.annotFontSize, cfg.annotRight);
        end

        hold(ax, 'off');

        if row < nRows
            set(ax,'XTickLabel',[],'XTickLabelMode','manual');
        end
        if col > 1
            set(ax,'YTickLabel',[],'YTickLabelMode','manual');
        end
    end

    if strcmp(cfg.legendMode,'shared')
        if ~all(tileAnnotFull)
            lgd = legend(axHandles(1), sharedHandles, sharedLabels, ...
                'Location', cfg.legendLocation, ...
                'FontSize',  cfg.legendFontSize, ...
                'Interpreter','tex');
            lgd.ItemTokenSize(2) = cfg.legendRowSpacing;
        else
            legend(axHandles(1), 'off');
        end
    end

    hFrame = axes('Units','normalized','Position',[0 0 1 1], ...
        'Visible','off','HitTest','off','HandleVisibility','off');
    uistack(hFrame,'bottom');

    plotAreaWidth = nCols*tW + (nCols-1)*cfg.colGap;
    xLabelX = (cfg.marginLeft + plotAreaWidth/2) / figW;
    xLabelY = (cfg.marginBottom - cfg.xLabelOffset) / figH;
    text(hFrame, xLabelX, xLabelY, 'Time (hours)', ...
        'Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontSize',cfg.axLabelSize, ...
        'Interpreter','tex');

    if strcmp(cfg.yMode,'normalized')
        yAxisLabel = 'Normalized Luminescence';
    else
        yAxisLabel = 'Luminescence (cps, \times10^{3})';
    end
    yLabelX = (cfg.marginLeft - cfg.yLabelOffset) / figW;
    yLabelY = cfg.marginBottom/figH + ...
              (sum(rowHeights) + cfg.tileGap*(nRows-1)) / figH / 2;
    text(hFrame, yLabelX, yLabelY, yAxisLabel, ...
        'Units','normalized','Rotation',90, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',cfg.axLabelSize,'Interpreter','tex');

    fig2.PaperUnits    = 'inches';
    fig2.PaperSize     = [figW, figH];
    fig2.PaperPosition = [0, 0, figW, figH];
    exportgraphics(fig2, tilesPDF, 'ContentType','vector');

    fprintf('Tiles saved: %s\n', tilesPDF);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- HELPER: readCol2_textscan --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vals = readCol2_textscan(filepath)
    vals = [];
    fid  = fopen(filepath, 'r');
    if fid == -1
        fprintf('  -> could not open file\n');
        return;
    end
    try
        fgetl(fid);
        C = textscan(fid, '%s %f', 'Delimiter', ',', 'EmptyValue', NaN);
        fclose(fid);
        if numel(C) >= 2 && ~isempty(C{2})
            v    = C{2};
            vals = v(isfinite(v));
        end
    catch ME
        fclose(fid);
        fprintf('  -> FAILED: %s\n', ME.message);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- HELPER: formatAxes --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function formatAxes(ax, time, numShadedRegions, tileShadingMatrix, yTickSpacing)
    yL       = ylim(ax);
    ax.YTick = yL(1):yTickSpacing:yL(2);
    ax.XTick = 0:24:max(time);
    ax.XMinorTick = 'off';
    ax.XAxis.MinorTickValues = 12:24:max(time)-12;
    grid(ax,'on'); box(ax,'off');
    ax.XGrid = 'on'; ax.YGrid = 'off'; ax.XMinorGrid = 'off';
    ax.GridColor       = [0 0 0]; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
    ax.MinorGridColor  = [0.5 0.5 0.5];
    ax.MinorGridLineStyle = '-'; ax.MinorGridAlpha = 0.2;
    ax.FontSize = 15;
    set(ax,'LineWidth',1.5);

    yLimits = ylim(ax);
    for i = 1:numShadedRegions
        fill(ax, ...
            [12+(i-1)*24, 24+(i-1)*24, 24+(i-1)*24, 12+(i-1)*24], ...
            [yLimits(1), yLimits(1), yLimits(2), yLimits(2)], ...
            'k','FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
    end

    if ~isempty(tileShadingMatrix)
        for r = 1:size(tileShadingMatrix,1)
            tS  = tileShadingMatrix(r,1);
            tE  = tileShadingMatrix(r,2);
            gry = tileShadingMatrix(r,3);
            fill(ax, [tS tE tE tS], ...
                [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
                [gry gry gry],'FaceAlpha',0.4,'EdgeColor','none', ...
                'HandleVisibility','off');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- HELPER: addAnnotationBox --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addAnnotationBox(ax, seriesIdx, seriesNames, colors, ...
    seriesPeriodMean, seriesPeriodStd, ...
    seriesPhaseMean,  seriesPhaseStd, ...
    fontSize, annotRight)

    lines      = {};
    lineColors = {};

    for k = 1:numel(seriesIdx)
        i     = seriesIdx(k);
        hasP  = isfinite(seriesPeriodMean(i));
        hasPh = isfinite(seriesPhaseMean(i));

        parts = {};
        if hasP
            parts{end+1} = sprintf('T: %.2f+/-%.2f', ...
                seriesPeriodMean(i), seriesPeriodStd(i));
        end
        if hasPh
            parts{end+1} = sprintf('Ph: %.2f+/-%.2f', ...
                seriesPhaseMean(i), seriesPhaseStd(i));
        end
        if isempty(parts)
            lineText = seriesNames{i};
        else
            lineText = sprintf('%s %s h', seriesNames{i}, strjoin(parts,' '));
        end
        lines{end+1}      = lineText; %#ok
        lineColors{end+1} = colors(i,:); %#ok
    end

    if isempty(lines), return; end

    xLim = xlim(ax);
    yLim = ylim(ax);

    axPos    = ax.Position;
    figH_in  = ax.Parent.Position(4);
    axH_in   = axPos(4) * figH_in;
    axH_data = diff(yLim);
    lineStep = (fontSize + 15) / 72 * axH_data / axH_in;

    nLines = numel(lines);
    boxH   = nLines * lineStep + lineStep * 0.3;

    swatchW  = 0.045 * diff(xLim);
    textXOff = swatchW + 0.015 * diff(xLim);

    maxChars = max(cellfun(@numel, lines));
    charWidth = fontSize * 0.55 / 72;
    axW_in   = ax.Position(3) * ax.Parent.Position(3);
    axW_data = diff(xLim);
    boxW     = maxChars * charWidth / axW_in * axW_data ...
               + swatchW + 0.07 * diff(xLim);
    boxW     = min(boxW, 0.95 * diff(xLim));

    yTop = yLim(2) - 0.04 * diff(yLim);
    if annotRight
        xPos = xLim(2) - 0.02 * diff(xLim) - boxW;
    else
        xPos = xLim(1) + 0.02 * diff(xLim);
    end

    fill(ax, ...
        [xPos, xPos+boxW, xPos+boxW, xPos], ...
        [yTop, yTop, yTop-boxH, yTop-boxH], ...
        [1.00 1.00 1.00], ...
        'EdgeColor',[0.4 0.4 0.4],'LineWidth',0.8, ...
        'HandleVisibility','off');

    for k = 1:nLines
        yPos = yTop - (k - 0.5) * lineStep;
        line(ax, ...
            [xPos + 0.01*diff(xLim), xPos + 0.01*diff(xLim) + swatchW], ...
            [yPos, yPos], ...
            'Color',lineColors{k},'LineWidth',2.5,'HandleVisibility','off');
        text(ax, xPos + 0.01*diff(xLim) + textXOff, yPos, lines{k}, ...
            'VerticalAlignment','middle', ...
            'HorizontalAlignment','left', ...
            'FontSize',fontSize, ...
            'Color',[0 0 0], ...
            'Interpreter','none', ...
            'FontWeight','normal');
        fprintf('  [annotation] %s\n', lines{k});
    end
end
