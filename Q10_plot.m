function Q10_plot(plotDefs, strainStyles, cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q10_plot   Period boxplot + Q10 engine — called from Q10_config.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p = 1:numel(plotDefs)
    pd = plotDefs(p);
    fprintf('\n=== Figure %d / %d : %s ===\n', p, numel(plotDefs), pd.outputName);
    buildFigure(pd, strainStyles, cfg);
end

end

% =========================================================================
function buildFigure(pd, strainStyles, cfg)
% =========================================================================
nTiles = numel(pd.tiles);
nCols = pd.nCols;
nRows = ceil(nTiles / nCols);

tW = cfg.tileWidth;
tH = cfg.tileHeight;
figW = cfg.marginL + nCols*tW + (nCols-1)*cfg.hGap + cfg.marginR;
figH = cfg.marginB + nRows*tH + (nRows-1)*cfg.vGap + cfg.marginT;

fig = figure('Units','inches','Position',[1 1 figW figH],'Color','w');

% Store axes handles so we can compare y-limits across columns after drawing
axHandles = gobjects(nRows, nCols);

for t = 1:nTiles
row = ceil(t / nCols);
col = mod(t-1, nCols) + 1;

xL = (cfg.marginL + (col-1)*(tW+cfg.hGap)) / figW;
yB = (cfg.marginB + (nRows-row)*(tH+cfg.vGap)) / figH;
wN = tW / figW;
hN = tH / figH;

ax = axes('Units','normalized','Position',[xL yB wN hN],'Parent',fig); %#ok
axHandles(row, col) = ax;

tile = pd.tiles{t};
st = strainStyles.(tile.strainKey);
tileTitle = '';
if isfield(pd,'tileTitles') && t <= numel(pd.tileTitles)
tileTitle = pd.tileTitles{t};
end

drawTile(ax, tile, st, cfg, tileTitle);

if row == nRows
xlabel(ax,'Temperature (°C)','FontSize',cfg.axLabelSize);
else
set(ax,'XTickLabel',[]); % hide x tick labels on non-bottom tiles
end
% Always show y tick labels on leftmost column
if col == 1
ylabel(ax,'Period (h)','FontSize',cfg.axLabelSize);
end
end

% After all tiles are drawn: for each row, compare y-limits of adjacent
% columns and only suppress right-tile tick labels if they match left tile.
for row = 1:nRows
for col = 2:nCols
axL = axHandles(row, col-1);
axR = axHandles(row, col);
if ~isgraphics(axL) || ~isgraphics(axR)
continue;
end
ylL = ylim(axL);
ylR = ylim(axR);
if isequal(ylL, ylR)
set(axR, 'YTickLabel', []); % same limits — suppress right labels
else
 ylabel(axR, 'Period (h)', 'FontSize', cfg.axLabelSize); % different — show them
end
end
end

if ~isempty(pd.title)
annotation(fig,'textbox',[0, 1-cfg.marginT/figH, 1, cfg.marginT/figH], ...
'String',pd.title,'HorizontalAlignment','center', ...
'VerticalAlignment','middle','FontSize',15,'FontWeight','bold', ...
'EdgeColor','none','Interpreter','tex');
end

fig.PaperUnits = 'inches';
fig.PaperSize = [figW figH];
fig.PaperPosition = [0 0 figW figH];
exportgraphics(fig,[pd.outputName '.svg'],'ContentType','vector');
print(fig,[pd.outputName '.png'],'-dpng','-r300');
fprintf(' Saved: %s.pdf | %s.png\n', pd.outputName, pd.outputName);

end

% =========================================================================
function drawTile(ax, tile, st, cfg, tileTitle)
% =========================================================================
allTemps   = zeros(0,1);
allPeriods = zeros(0,1);

for k = 1:numel(tile.files)
    fpath = tile.files{k};
    if ~isfile(fpath)
        warning('Q10_plot:fileNotFound','Q10_plot: file not found:\n  %s', fpath);
        continue;
    end
    vals = readPeriods(fpath);
    if isempty(vals)
        warning('Q10_plot:noNumericData','Q10_plot: no numeric data in:\n  %s', fpath);
        continue;
    end
    n          = numel(vals);
    allTemps   = [allTemps;   repmat(double(tile.temps(k)), n, 1)]; %#ok
    allPeriods = [allPeriods; double(vals(:))];                     %#ok
end

if isempty(allTemps)
    error('Q10_plot:noData','No data loaded for strain "%s". Check file paths.', tile.strainKey);
end

% ------------------------------------------------------------------
% blankTemps: show x-axis tick but NO data box/scatter
% ------------------------------------------------------------------
if isfield(st,'blankTemps') && ~isempty(st.blankTemps)
    blank = st.blankTemps(:);
else
    blank = [];
end
keep            = ~ismember(allTemps, blank);
allTemps_plot   = allTemps(keep);
allPeriods_plot = allPeriods(keep);

% allTempsAll  — ALL temps incl. blanked  → x-axis ticks & plot positions
% uniqueTemps_plot — only temps with data → stats, Q10, regression
allTempsAll      = sort(unique([allTemps; blank]));
nT               = numel(allTempsAll);
positions        = (1:nT)';

uniqueTemps_plot = sort(unique(allTemps_plot));
nT_plot          = numel(uniqueTemps_plot);

% axis position index for each data temp
[~, axIdx_plot] = ismember(uniqueTemps_plot, allTempsAll);  % nT_plot x 1

% Stats over non-blank data only
stats_mean = zeros(nT_plot,1);
stats_std  = zeros(nT_plot,1);
for k = 1:nT_plot
    v             = allPeriods_plot(allTemps_plot == uniqueTemps_plot(k));
    stats_mean(k) = mean(v,'omitnan');
    stats_std(k)  = std(v,'omitnan');
end

% Regression (excludeTemps applied within non-blank data)
modelChoice  = st.model;
excludeTemps = st.excludeTemps;
if strcmp(modelChoice,'linearsubset') && ~isempty(excludeTemps)
    mask = ~ismember(allTemps_plot, excludeTemps);
else
    mask = true(size(allTemps_plot));
end
T_fit = allTemps_plot(mask);
P_fit = allPeriods_plot(mask);

switch modelChoice
    case {'linear','linearsubset'}, polyDeg = 1;
    case 'quadratic',               polyDeg = 2;
    case 'cubic',                   polyDeg = 3;
    otherwise
        polyDeg = 1;
        warning('Q10_plot:badModel','Unknown model "%s"; using linear.',modelChoice);
end

p_coef   = polyfit(T_fit, P_fit, polyDeg);
y_fitted = polyval(p_coef, T_fit);
SS_res   = sum((P_fit - y_fitted).^2);
SS_tot   = sum((P_fit - mean(P_fit)).^2);
R2       = 1 - SS_res/SS_tot;
n_fit    = numel(T_fit);
nCoef    = polyDeg + 1;
df       = n_fit - nCoef;
s2       = SS_res / max(df,1);

Xd = zeros(n_fit,nCoef);
for j = 1:nCoef
    Xd(:,j) = T_fit .^ (polyDeg-(j-1));
end
covBeta = s2 * pinv(Xd'*Xd);

se_slope   = sqrt(covBeta(polyDeg,polyDeg));
t_slope    = p_coef(polyDeg) / se_slope;
pval_slope = 2*(1 - tcdf(abs(t_slope), df));

pval_me = NaN;
try
    tblFit   = table(T_fit, P_fit,'VariableNames',{'Temperature','Period'});
    me_null  = fitlme(tblFit,'Period ~ 1 + (1|Temperature)');
    me_model = fitlme(tblFit,'Period ~ Temperature + (1|Temperature)');
    [~,pval_me] = compare(me_null, me_model);
catch
end

T_min   = min(T_fit);  T_max = max(T_fit);
tau_min = double(polyval(p_coef, T_min));
tau_max = double(polyval(p_coef, T_max));
Q10_reg = (tau_min/tau_max) .^ (10/(T_max-T_min));

exponent = 10/(T_max-T_min);
grad = zeros(1,nCoef);
for i = 1:nCoef
    pw       = polyDeg-(i-1);
    grad(i)  = exponent*(tau_min/tau_max)^(exponent-1) * ...
               (T_min^pw/tau_max - T_max^pw*tau_min/tau_max^2);
end
Q10_reg_std = sqrt(max(0, grad*covBeta*grad'));

pairT1    = uniqueTemps_plot(1:end-1);
pairT2    = uniqueTemps_plot(2:end);
Q10_pw    = zeros(numel(pairT1),1);
Q10_pw_sd = zeros(numel(pairT1),1);
for i = 1:numel(pairT1)
    tau1 = stats_mean(uniqueTemps_plot==pairT1(i));
    tau2 = stats_mean(uniqueTemps_plot==pairT2(i));
    sd1  = stats_std(uniqueTemps_plot==pairT1(i));
    sd2  = stats_std(uniqueTemps_plot==pairT2(i));
    Q10_pw(i)    = (tau1/tau2).^(10/(pairT2(i)-pairT1(i)));
    Q10_pw_sd(i) = Q10_pw(i)*sqrt((sd1/tau1)^2+(sd2/tau2)^2);
end

fprintf('  --- %s ---\n', tile.strainKey);
fprintf('  Model: %s  deg=%d  R2=%.4f\n', modelChoice, polyDeg, R2);
if polyDeg==1, fprintf('  Slope: %.4f h/degC  p=%.4g\n', p_coef(1), pval_slope); end
fprintf('  ME p-value: %.4g\n', pval_me);
fprintf('  Q10 (%d->%d C): %.3f +/- %.3f\n', T_min, T_max, Q10_reg, Q10_reg_std);
for i = 1:numel(pairT1)
    fprintf('    %d->%d C: %.3f +/- %.3f\n', pairT1(i),pairT2(i),Q10_pw(i),Q10_pw_sd(i));
end

% ------------------------------------------------------------------
% Colormap — nT entries (one per axis tick, including blanks)
% ------------------------------------------------------------------
hTmp = figure('Visible','off');
if isempty(st.colors)
    fullJet = flipud(colormap(hTmp,'jet'));
    subJet  = fullJet(24:232,:);
    cmapAll = interp1(linspace(0,1,size(subJet,1)), subJet, linspace(0,1,nT));
else
    cmapAll = st.colors;
    if size(cmapAll,1) ~= nT
        warning('Q10_plot:colorMismatch', ...
            'st.colors has %d rows but nT=%d; using jet.', size(cmapAll,1), nT);
        fullJet = flipud(colormap(hTmp,'jet'));
        subJet  = fullJet(24:232,:);
        cmapAll = interp1(linspace(0,1,size(subJet,1)), subJet, linspace(0,1,nT));
    end
end
close(hTmp);

% ------------------------------------------------------------------
% Draw — manual boxes so blank positions are truly empty gaps
% ------------------------------------------------------------------
axes(ax);
hold(ax,'on');
lw      = cfg.boxLineWidth;
bw      = 0.4;   % half-width of box in axis units

for k = 1:nT_plot
    xc   = axIdx_plot(k);          % correct axis position (respects gaps)
    vals = allPeriods_plot(allTemps_plot == uniqueTemps_plot(k));
    col  = cmapAll(xc,:);

    q1   = quantile(vals, 0.25);
    q3   = quantile(vals, 0.75);
    med  = median(vals);
    iqr_ = q3 - q1;
    wLo  = max(vals(vals >= q1 - 1.5*iqr_));
    wHi  = min(vals(vals <= q3 + 1.5*iqr_));
    outs = vals(vals < q1-1.5*iqr_ | vals > q3+1.5*iqr_);

    % whiskers first (drawn behind box)
    plot(ax,[xc xc],[wLo q1],'k-','LineWidth',lw);
    plot(ax,[xc xc],[q3 wHi],'k-','LineWidth',lw);
    % whisker caps
    plot(ax,[xc-bw/2 xc+bw/2],[wLo wLo],'k-','LineWidth',lw);
    plot(ax,[xc-bw/2 xc+bw/2],[wHi wHi],'k-','LineWidth',lw);
    % filled box on top of whiskers
    patch([xc-bw xc+bw xc+bw xc-bw xc-bw], ...
          [q1 q1 q3 q3 q1], col, ...
          'FaceAlpha',0.85,'EdgeColor','k','LineWidth',lw,'Parent',ax);
    % median line on top of box
    plot(ax,[xc-bw xc+bw],[med med],'k-','LineWidth',lw);
    % outliers
    if ~isempty(outs)
        scatter(ax, repmat(xc,size(outs)), outs, cfg.markerSize*0.8, ...
            'o','MarkerEdgeColor','k','LineWidth',1);
    end

    % jittered scatter
    rng(42+k);
    jitter = (rand(size(vals))-0.5)*2*cfg.jitterWidth;
    scatter(ax, xc+jitter, vals, cfg.markerSize, ...
        [0.45 0.45 0.45],'filled','MarkerFaceAlpha',0.55);
end

% Regression curve mapped onto full axis positions
x_temp_vec = linspace(T_min, T_max, 200);
x_pos_vec  = interp1(allTempsAll, positions, x_temp_vec,'linear','extrap');
y_curve    = polyval(p_coef, x_temp_vec);
plot(ax, x_pos_vec, y_curve, 'r--','LineWidth',2.5);

switch polyDeg
    case 1, eq_str = sprintf('y=%.3f+%.3f*T', p_coef(2), p_coef(1));
    case 2, eq_str = sprintf('y=%.2f+%.3f*T+%.4f*T^2', p_coef(3),p_coef(2),p_coef(1));
    case 3, eq_str = sprintf('y=%.2f+%.3f*T+%.4f*T^2+%.5f*T^3', p_coef(4),p_coef(3),p_coef(2),p_coef(1));
end


ann = sprintf('Q10(%d->%dC)=%.3f+/-%.3f\n%s\np=%.3f  R^2=%.3f', ...
    T_min, T_max, Q10_reg, Q10_reg_std, eq_str, pval_slope, R2);


text(ax, 0.03, 0.97, ann, 'Units','normalized','FontSize',cfg.annotFontSize, ...
    'VerticalAlignment','top','HorizontalAlignment','left', ...
    'BackgroundColor',[1 1 1 0.7],'EdgeColor',[0.7 0.7 0.7]);

% X-axis: ALL temps appear (blanked ones just have no box)
xticks(ax, 1:nT);
xticklabels(ax, arrayfun(@num2str, allTempsAll,'UniformOutput',false));
xlim(ax, [0.5, nT+0.5]);
ylim(ax, [st.periodMin, st.periodMax]);
ax.FontSize   = cfg.tickFontSize;
ax.LineWidth  = 1.8;
ax.TickDir    = 'out';
ax.TickLength = [0.02 0.02];
ax.Layer      = 'top';
box(ax,'on');
grid(ax,'off');
if ~isempty(tileTitle)
    title(ax, tileTitle,'FontSize',13,'FontWeight','bold','Interpreter','tex');
end
hold(ax,'off');

end

% =========================================================================
function periods = readPeriods(filepath)
% Read second column of a 2-column CSV with one header row.
% Uses fgetl — robust to hyphenated column headers like Segment_28-72h.
periods = [];
fid = fopen(filepath,'r');
if fid == -1
    error('Q10_plot:cannotOpen','Cannot open file: %s', filepath);
end
fgetl(fid);   % discard header
tline = fgetl(fid);
while ischar(tline)
    tline = strtrim(tline);
    if ~isempty(tline)
        comma = find(tline==',',1,'first');
        if ~isempty(comma)
            v = str2double(strtrim(tline(comma+1:end)));
            if ~isnan(v), periods(end+1,1) = v; end %#ok
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
end
