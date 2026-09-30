function [allPeriods, allPhases, allAdjustedPhases, allFittedCurves, allChirpRates, ...
          segmentResults, ...
          dampingCoeffs, meanDamping, stdDamping, ...
          phaseDriftSlopes, meanPhaseDrift, stdPhaseDrift] = ...
    analyzeAndPlotSegments(data, time, segments, samplesList, colors, ...
                           phaseOffsetHours, sampleName, ...
                           plotIndividualFittedCurves, plotAverageFittedCurve)
% analyzeAndPlotSegments  Fit chirped damped cosines to circadian data
%                         segments and report period, acrophase, damping,
%                         chirp rate, phase coherence, and drift.
%
%   [...] = analyzeAndPlotSegments(data, time, segments, samplesList,
%               colors, phaseOffsetHours, sampleName,
%               plotIndividualFittedCurves, plotAverageFittedCurve)
%
%   Model fitted per segment (fitCosineCurve):
%     y(t) = A*exp(-lambda*u)*cos(omega0*u + 0.5*kappa*u^2 + phi0) + B
%   where u = t - t0 (t0 = segment start).
%
%   Inputs
%   ------
%   data                       - Struct; each field is a sample time series.
%   time                       - Time vector (hours).
%   segments                   - Nx2 [segStart, segEnd] matrix (hours).
%   samplesList                - Cell array of field names in data.
%   colors                     - numSamples x 3 RGB colour matrix.
%   phaseOffsetHours           - Reference time (h) for acrophase reporting.
%                                  0           => relative to ZT0/CT0
%                                  segStart    => relative to segment start
%                                  segStart-24 => relative to prior day
%   sampleName                 - String prefix for all CSV output files.
%   plotIndividualFittedCurves - Logical: overlay each sample's fit curve.
%   plotAverageFittedCurve     - Logical: overlay mean +/- SD curve per
%                                segment with period/phase annotation.
%
%   Outputs
%   -------
%   allPeriods        - numSegments x numSamples fitted periods (h, at midpoint).
%   allPhases         - numSegments x numSamples initial phases phi0 (rad, [0,2pi)).
%   allAdjustedPhases - numSegments x numSamples acrophases (h, [0,T_mid)).
%   allFittedCurves   - numSegments x 1 cell; each cell is numSamples x 100
%                       fitted-curve values on a common 100-point grid.
%   allChirpRates     - numSegments x numSamples chirp rates kappa (rad/h^2).
%                         > 0 : period shortening
%                         < 0 : period lengthening
%                         = 0 : constant period
%   segmentResults    - Struct array per segment:
%                         segmentStart, segmentEnd
%                         meanPeriod, stdPeriod  (at midpoint)
%                         meanChirpRate, stdChirpRate
%                         meanPhase (circular), stdPhase (circular)
%                         meanLambda, stdLambda
%                         coherence, nValidSamples
%   dampingCoeffs     - numSamples x 1 amplitude-decay rate (h^-1).
%   meanDamping       - Mean of valid damping coefficients.
%   stdDamping        - SD of valid damping coefficients.
%   phaseDriftSlopes  - numSamples x 1 phase drift (h/h).
%   meanPhaseDrift    - Mean of valid drift slopes.
%   stdPhaseDrift     - SD of valid drift slopes.
%
%   CSV files written (prefix = sampleName)
%   ----------------------------------------
%   _individual_periods.csv      - Period at segment midpoint (T_mid) per sample.
%   _individual_periods_mean.csv  - Time-averaged period (T_mean) per sample.
%   _individual_deltaT.csv        - Period change DeltaT=T_end-T_start per sample.
%   _individual_phases.csv      - Acrophase (h) per sample per segment.
%   _individual_chirprates.csv  - Chirp rate kappa (rad/h^2) per sample per seg.
%   _individual_lambdas.csv     - Per-fit lambda (h^-1) per sample per segment.
%   _individual_peak_times.csv  - First two raw peaks per sample per segment.
%   _damping_coefficients.csv   - Cross-segment damping per sample.
%   _phase_coherence.csv        - Phase coherence (1/circ_std) per segment.
%   _phase_drift.csv            - Phase drift slope per sample.

% =========================================================================
% Initialise
% =========================================================================
numSegments = size(segments, 1);
numSamples  = length(samplesList);

% Pre-allocate result matrices  (NaN = fit not attempted or failed)
segmentPeriods        = nan(numSegments, numSamples);
segmentPhases         = nan(numSegments, numSamples);   % phi0, rad
segmentAdjustedPhases = nan(numSegments, numSamples);   % acrophase, h
segmentAmplitudes     = nan(numSegments, numSamples);
segmentDampings       = nan(numSegments, numSamples);   % lambda, h^-1
segmentChirpRates     = nan(numSegments, numSamples);   % kappa, rad/h^2
segmentPeriodMeans    = nan(numSegments, numSamples);   % T_mean (time-avg), h
segmentDeltaTs        = nan(numSegments, numSamples);   % DeltaT = T_end-T_start, h

% Common 100-point time grid per segment for averaged fitted curves
commonTimePerSegment = cell(numSegments, 1);
segmentFittedCurves  = cell(numSegments, 1);
for i = 1:numSegments
    commonTimePerSegment{i} = linspace(segments(i,1), segments(i,2), 100);
    segmentFittedCurves{i}  = nan(numSamples, 100);
end

% Segment midpoints -- used for damping/drift time axis
segmentMidTimes = mean(segments, 2);   % numSegments x 1

% Raw peak storage (first two peaks per segment per sample)
rawFirstTwoPeakTimesPerSegment = cell(numSegments, 1);
for i = 1:numSegments
    rawFirstTwoPeakTimesPerSegment{i} = nan(numSamples, 2);
end

% Phase trajectory for cross-segment drift analysis
samplePhaseTrajectories = nan(numSamples, numSegments);
validSegmentFlags       = false(numSamples, numSegments);

% Output pre-init
dampingCoeffs    = nan(numSamples, 1);
phaseDriftSlopes = nan(numSamples, 1);
meanDamping    = NaN;  stdDamping    = NaN;
meanPhaseDrift = NaN;  stdPhaseDrift = NaN;

% Phase coherence (computed after all samples are fitted)
phaseCoherence = nan(numSegments, 1);

% Storage for the interactive individual-curves floating window
segmentRawTimeStore  = cell(numSegments, numSamples);
segmentRawDataStore  = cell(numSegments, numSamples);
segmentFitLabelStore = cell(numSegments, numSamples);


% =========================================================================
% Auto-adjust segment end times: 66.7% into the half-cycle containing
% the segment end.
%   (A) Ends in PEAK-to-TROUGH: target = last_peak   + 0.667*(trough-peak)
%   (B) Ends in TROUGH-to-PEAK: target = last_trough + 0.667*(peak-trough)
% =========================================================================
segments_adj = segments;

for i = 1:numSegments
    segStart_i = segments(i, 1);
    segEnd_i   = segments(i, 2);
    mask_i     = (time > segStart_i) & (time <= segEnd_i);
    segTime_i  = time(mask_i);
    if numel(segTime_i) < 6
        continue;
    end

    % Step 1: median FFT period
    T_fft_all = nan(numSamples, 1);
    for jj = 1:numSamples
        sd   = data.(samplesList{jj}); sd = sd(mask_i);
        Fs_i = 1 / mean(diff(segTime_i));
        Nf   = 2^nextpow2(numel(sd));
        Yf   = fft(sd - mean(sd), Nf);
        Pf   = 2*abs(Yf(1:Nf/2+1)/Nf);
        ff   = Fs_i*(0:Nf/2)/Nf;
        vi   = ff >= 1/30 & ff <= 1/18;
        if any(vi)
            [~,im] = max(Pf(vi)); fv = ff(vi);
            T_fft_all(jj) = 1/fv(im);
        end
    end
    T_est = median(T_fft_all(~isnan(T_fft_all)));
    if isnan(T_est) || T_est <= 0
        continue;
    end

    % Step 2: average signal
    avgSig = zeros(numel(segTime_i), 1);
    for jj = 1:numSamples
        sd     = data.(samplesList{jj});
        avgSig = avgSig + sd(mask_i);
    end
    avgSig = avgSig / numSamples;

    % Step 3: find all peaks and troughs
    dt      = mean(diff(segTime_i));
    minDist = max(1, round(T_est * 0.4 / dt));
    [~, pkLocs] = findpeaks( avgSig, 'MinPeakDistance', minDist);
    [~, trLocs] = findpeaks(-avgSig, 'MinPeakDistance', minDist);

    % Step 4: last peak and last trough before segment end
    if isempty(pkLocs)
        t_last_peak = segStart_i - 1;
    else
        t_last_peak = segTime_i(pkLocs(end));
    end
    if isempty(trLocs)
        t_last_trough = segStart_i - 1;
    else
        t_last_trough = segTime_i(trLocs(end));
    end

    % Step 5: interval type and next landmark
    if t_last_peak > t_last_trough
        % Case A: PEAK-to-TROUGH interval
        afterPeak   = pkLocs(end)+1 : numel(avgSig);
        t_next_land = [];
        if numel(afterPeak) >= 3
            [~, trAfter] = findpeaks(-avgSig(afterPeak));
            if ~isempty(trAfter)
                t_next_land = segTime_i(afterPeak(trAfter(1)));
            end
        end
        if isempty(t_next_land)
            t_next_land = t_last_peak + T_est/2;
        end
        t_anchor     = t_last_peak;
        intervalType = 'peak->trough';
    elseif t_last_trough > segStart_i - 1
        % Case B: TROUGH-to-PEAK interval
        afterTrough = trLocs(end)+1 : numel(avgSig);
        t_next_land = [];
        if numel(afterTrough) >= 3
            [~, pkAfter] = findpeaks(avgSig(afterTrough));
            if ~isempty(pkAfter)
                t_next_land = segTime_i(afterTrough(pkAfter(1)));
            end
        end
        if isempty(t_next_land)
            t_next_land = t_last_trough + T_est/2;
        end
        t_anchor     = t_last_trough;
        intervalType = 'trough->peak';
    else
        % Fallback: no landmarks detected
        excess      = mod(segEnd_i - segStart_i, T_est);
        t_anchor    = segEnd_i - excess;
        t_next_land = t_anchor + T_est/2;
        intervalType = 'fallback';
    end

    % Step 6: 66.7% target
    t_target = t_anchor + 0.667 * (t_next_land - t_anchor);

    % Step 7: snap and validate
    [~, iSnap] = min(abs(time - t_target));
    t_adj = time(iSnap);
    if abs(t_adj - segEnd_i) <= T_est/2 && t_adj > segStart_i + T_est
        segments_adj(i, 2) = t_adj;
        fprintf('Segment %d end adjusted: %.1f -> %.1f h (%s)\n', ...
            i, segEnd_i, t_adj, intervalType);
    else
        fprintf('Segment %d end unchanged: %.1f h (%s, target=%.1f h)\n', ...
            i, segEnd_i, intervalType, t_target);
    end
end

% =========================================================================
% Rebuild common time grid using adjusted segment end times
% =========================================================================

for i = 1:numSegments
    commonTimePerSegment{i} = linspace(segments_adj(i,1), segments_adj(i,2), 100);
    segmentFittedCurves{i}  = nan(numSamples, 100);
end

% =========================================================================
% Main fitting loop
% =========================================================================
for j = 1:numSamples
    samplesName   = samplesList{j};
    currentSample = data.(samplesName);

    for i = 1:numSegments
        segmentStart = segments_adj(i, 1);
        segmentEnd   = segments_adj(i, 2);

        mask        = (time > segmentStart) & (time <= segmentEnd);
        timeSegment = time(mask);
        dataSegment = currentSample(mask);

        if numel(timeSegment) < 5
            continue;
        end

        % --- Raw peak detection (first two peaks in unfiltered data) ---
        [~, locs] = findpeaks(dataSegment);
        if numel(locs) >= 1
            rawFirstTwoPeakTimesPerSegment{i}(j,1) = ...
                mod(timeSegment(locs(1)) - phaseOffsetHours, 24);
        end
        if numel(locs) >= 2
            rawFirstTwoPeakTimesPerSegment{i}(j,2) = ...
                mod(timeSegment(locs(2)) - phaseOffsetHours, 24);
        end

        % --- Chirped damped cosine fit ---
        [amplitude, period, phase, adjustedPhase, offset, damping, chirpRate, ~] = ...
            fitCosineCurve(timeSegment, dataSegment, [18, 30], phaseOffsetHours);

        if isnan(amplitude)
            continue;
        end

        % Store fit parameters
        segmentAmplitudes(i,j)     = amplitude;
        segmentPeriods(i,j)        = period;       % T at midpoint, h
        segmentPhases(i,j)         = phase;         % phi0, rad [0,2pi)
        segmentAdjustedPhases(i,j) = adjustedPhase; % acrophase, h [0,T_mid)
        segmentDampings(i,j)       = damping;        % lambda, h^-1
        segmentChirpRates(i,j)     = chirpRate;      % kappa, rad/h^2

        % T_mean (time-averaged) and DeltaT for this sample/segment
        t0_s   = timeSegment(1);
        segD   = timeSegment(end) - t0_s;
        tmid_s = t0_s + segD/2;
        om_mid_s = 2*pi / period;
        om0_s    = max(om_mid_s - chirpRate*(tmid_s - t0_s), 2*pi/35);
        omEnd_s  = max(om0_s + chirpRate*segD, 2*pi/35);
        segmentDeltaTs(i,j) = 2*pi/omEnd_s - 2*pi/om0_s;
        if abs(chirpRate) > 1e-9
            segmentPeriodMeans(i,j) = (2*pi/(chirpRate*segD)) * log(omEnd_s/om0_s);
        else
            segmentPeriodMeans(i,j) = 2*pi / om0_s;
        end

        % Reconstruct fitted curve on segment time grid
        % omega0 recovered from: omega_mid = omega0 + kappa*(t_mid - t0)
        t0_seg    = timeSegment(1);
        t_mid_seg = mean(timeSegment);
        omega_mid = 2*pi / period;
        omega0    = omega_mid - chirpRate * (t_mid_seg - t0_seg);
        u         = timeSegment - t0_seg;
        fittedCurve = amplitude * exp(-damping * u) .* ...
                      cos(omega0*u + 0.5*chirpRate*u.^2 + phase) + offset;

        % R-squared
        residuals = dataSegment - fittedCurve;
        rSquared  = 1 - sum(residuals.^2) / ...
                        sum((dataSegment - mean(dataSegment)).^2);

        % Interpolate onto 100-point common grid for segment averaging
        segmentFittedCurves{i}(j,:) = interp1(timeSegment, fittedCurve, ...
            commonTimePerSegment{i}, 'linear', 'extrap');

        % Store for phase drift analysis
        samplePhaseTrajectories(j,i) = adjustedPhase;
        validSegmentFlags(j,i)       = true;

        % --- Store segment data + fit info (always needed for floating window) ---
        segmentRawTimeStore{i,j} = timeSegment;
        segmentRawDataStore{i,j} = dataSegment;
        labelStr = samplesName(max(1, end-4):end);
        dTdt     = -chirpRate * period^2 / (2*pi);
        segmentFitLabelStore{i,j} = sprintf( ...
            '%s  T=%.1fh  phi=%.1fh  dT/dt=%.3f  R^2=%.2f', ...
            labelStr, period, adjustedPhase, dTdt, rSquared);

        % --- Plot individual fit on the main axes when line 10 = 1 ---
        if plotIndividualFittedCurves
            plot(timeSegment, fittedCurve, '--', 'LineWidth', 1, ...
                'Color', colors(j,:), ...
                'DisplayName', sprintf( ...
                    '%s: T=%.1fh  phi=%.1fh  dT/dt=%.3f  R^2=%.2f', ...
                    labelStr, period, adjustedPhase, dTdt, rSquared));
        end

    end % segment loop
end % sample loop

% =========================================================================
% Phase coherence (computed after all fits)
% =========================================================================
for i = 1:numSegments
    validIdx    = ~isnan(segmentAdjustedPhases(i,:));
    validPhases = segmentAdjustedPhases(i, validIdx);

    if numel(validPhases) < 2
        continue;
    end

    meanPeriod_i = mean(segmentPeriods(i, ~isnan(segmentPeriods(i,:))));
    if isnan(meanPeriod_i); meanPeriod_i = 24; end

    angles           = 2*pi * validPhases / meanPeriod_i;
    phaseCoherence(i) = 1 / circ_std(angles);
end

% =========================================================================
% Average fitted curve with annotation  (plotAverageFittedCurve)
% =========================================================================
segmentColors = lines(numSegments);

sunYellow = [1, 0.843, 0];   % #FFD700

if plotAverageFittedCurve
    for i = 1:numSegments
        commonTime = commonTimePerSegment{i};
        curves     = segmentFittedCurves{i};
        color      = segmentColors(i,:);

        meanCurve = mean(curves, 1, 'omitnan');
        stdCurve  = std( curves, 0, 1, 'omitnan');

        % Mean curve - sun yellow
        plot(commonTime, meanCurve, '-', 'LineWidth', 2, 'Color', sunYellow, ...
            'DisplayName', sprintf('Chirped fit %d-%dh', ...
                segments(i,1), segments(i,2)));

        % SD shaded band - sun yellow
        fill([commonTime, fliplr(commonTime)], ...
             [meanCurve + stdCurve, fliplr(meanCurve - stdCurve)], ...
             sunYellow, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');

        % -----------------------------------------------------------------
        % Annotation statistics
        % -----------------------------------------------------------------
        validPeriods = segmentPeriods(i, ~isnan(segmentPeriods(i,:)));
        validKappas  = segmentChirpRates(i, ~isnan(segmentChirpRates(i,:)));
        validLambdas = segmentDampings(i, ~isnan(segmentDampings(i,:)));
        validAPhases = segmentAdjustedPhases(i, ~isnan(segmentAdjustedPhases(i,:)));

        meanPeriod = mean(validPeriods);   % used for circular acrophase normalisation
        stdPeriod = std(validPeriods);   % SD of T_mid across samples
        meanKappa  = mean(validKappas);   stdKappa  = std(validKappas);
        meanLambda = mean(validLambdas);  stdLambda = std(validLambdas);

        % Time-averaged period and period change per sample
        segDur   = segments(i,2) - segments(i,1);
        t0seg    = segments(i,1);
        tmidSeg  = segments(i,1) + segDur/2;
        vIdx     = find(~isnan(segmentPeriods(i,:)) & ~isnan(segmentChirpRates(i,:)));
        Tmean_v  = nan(size(vIdx));
        DeltaT_v = nan(size(vIdx));
        for kk = 1:numel(vIdx)
            jj     = vIdx(kk);
            kap    = segmentChirpRates(i,jj);
            om_mid = 2*pi / segmentPeriods(i,jj);
            om0    = max(om_mid - kap*(tmidSeg - t0seg), 2*pi/35);
            om_end = max(om0 + kap*segDur, 2*pi/35);
            T0     = 2*pi / om0;
            Tend   = 2*pi / om_end;
            DeltaT_v(kk) = Tend - T0;
            if abs(kap) > 1e-9
                % Exact time-average: (1/D)*integral_0^D 2pi/(om0+kap*u) du
                Tmean_v(kk) = (2*pi/(kap*segDur)) * log(om_end/om0);
            else
                Tmean_v(kk) = T0;   % kappa~0: constant period
            end
        end
        meanTmean  = mean(Tmean_v, 'omitnan');  stdTmean  = std(Tmean_v,  [], 'omitnan');
        meanDeltaT = mean(DeltaT_v,'omitnan');  stdDeltaT = std(DeltaT_v, [], 'omitnan');

        % Circular acrophase mean/SD
        angles       = 2*pi * validAPhases / meanPeriod;
        meanPhaseRad = atan2(mean(sin(angles)), mean(cos(angles)));
        meanPhase    = mod(meanPhaseRad * meanPeriod / (2*pi), meanPeriod);
        stdPhase     = circ_std(angles) * meanPeriod / (2*pi);

        % Raw peak statistics
        rawPeaks = rawFirstTwoPeakTimesPerSegment{i};
        peak1    = rawPeaks(~isnan(rawPeaks(:,1)), 1);
        peak2    = rawPeaks(~isnan(rawPeaks(:,2)), 2);
        if ~isempty(peak1); meanPeak1 = mean(peak1); stdPeak1 = std(peak1);
        else;               meanPeak1 = NaN;          stdPeak1 = NaN; end
        if ~isempty(peak2); meanPeak2 = mean(peak2); stdPeak2 = std(peak2);
        else;               meanPeak2 = NaN;          stdPeak2 = NaN; end

        % % Annotation position: just above mean-curve peak
        % [~, idxMax] = max(meanCurve);
        % xLabel = commonTime(idxMax);
        % yLabel = meanCurve(idxMax) + 0.08 * range(meanCurve);

        % Annotation position: just above the mean curve at the segment end
        xLabel = commonTime(end);                              % x = segment end time
        yLabel = max(meanCurve) + 0.08 * range(meanCurve);    % y = above the curve's global peak height



        phaseLabel = sprintf([ ...
            'T_{mid} = %.2f +/- %.2f h\n' ...   % <-- ADD THIS LINE
            'T_{mean} = %.2f +/- %.2f h\n' ...
            'dT       = %.2f +/- %.2f h\n' ...
            'kappa   = %.5f +/- %.5f rad/h^2\n' ...
            'phi     = %.2f +/- %.2f h\n' ...
            'lambda  = %.4f +/- %.4f h^-1\n' ...
            'Peak1   = %.2f +/- %.2f h\n' ...
            'Peak2   = %.2f +/- %.2f h'], ...
            meanPeriod, stdPeriod, ...           % <-- ADD THIS
            meanTmean,  stdTmean,   ...
            meanDeltaT, stdDeltaT,  ...
            meanKappa,  stdKappa,   ...
            meanPhase,  stdPhase,   ...
            meanLambda, stdLambda,  ...
            meanPeak1,  stdPeak1,   ...
            meanPeak2,  stdPeak2);

        text(xLabel, yLabel, phaseLabel, ...
            'Color', [0,0,0], 'FontSize', 10, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'BackgroundColor', 'w', 'EdgeColor', [0,0,0], 'Margin', 2);
    end
end


% =========================================================================
% Interactive floating window - individual fitted curves + raw data
% Always created regardless of plotIndividualFittedCurves value.
% Line 10 only controls overlay on the main figure / PDF.
% Axes tagged 'SegmentAxes_N' so findobj locates them after .fig reload.
% =========================================================================
listFrac = 0.22;
figH_in  = max(5, 3.0 * numSegments);

hFloat = figure( ...
    'Name',        'Individual Fitted Curves', ...
    'NumberTitle', 'off', ...
    'Units',       'inches', ...
    'Position',    [12, 1, 11, figH_in], ...
    'Color',       'w');

axLeft = listFrac + 0.04;
axW    = 1 - axLeft - 0.03;
axH    = (0.88 / numSegments) - 0.02;
for i = 1:numSegments
    axBot = 0.08 + (numSegments - i) * (0.88 / numSegments);
    axes('Parent',   hFloat, ...
         'Units',    'normalized', ...
         'Position', [axLeft, axBot, axW, axH], ...
         'Tag',      sprintf('SegmentAxes_%d', i));
end

ud.segmentRawTimeStore  = segmentRawTimeStore;
ud.segmentRawDataStore  = segmentRawDataStore;
ud.segmentFittedCurves  = segmentFittedCurves;
ud.commonTimePerSegment = commonTimePerSegment;
ud.segmentFitLabelStore = segmentFitLabelStore;
ud.samplesList          = samplesList;
ud.colors               = colors;
ud.segments             = segments;
ud.numSegments          = numSegments;
ud.fullTime             = time;
ud.fullData             = cell(numSamples, 1);
for jj = 1:numSamples
    ud.fullData{jj} = data.(samplesList{jj});
end
set(hFloat, 'UserData', ud);

uicontrol(hFloat, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [0.01, 0.97, listFrac-0.02, 0.025], ...
    'String', 'Select sample:', 'FontSize', 8, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

hList = uicontrol(hFloat, 'Style', 'listbox', ...
    'Units', 'normalized', 'Position', [0.01, 0.08, listFrac-0.02, 0.88], ...
    'String', samplesList, 'Value', 1, 'FontSize', 8, ...
    'BackgroundColor', [0.97 0.97 0.97], ...
    'Callback', @(src, ~) ifc_updatePlot(src));

ifc_updatePlot(hList);
savefig(hFloat, strcat(sampleName, '_individual_fits.fig'));

% =========================================================================
% Populate segmentResults struct
% =========================================================================
segmentResults = struct();
for i = 1:numSegments
    validPeriods = segmentPeriods(i, ~isnan(segmentPeriods(i,:)));
    validKappas  = segmentChirpRates(i, ~isnan(segmentChirpRates(i,:)));
    validLambdas = segmentDampings(i, ~isnan(segmentDampings(i,:)));
    validAPhases = segmentAdjustedPhases(i, ~isnan(segmentAdjustedPhases(i,:)));

    % Circular acrophase stats
    if ~isempty(validAPhases) && ~isempty(validPeriods)
        meanPer_i  = mean(validPeriods);
        if isnan(meanPer_i); meanPer_i = 24; end
        angs_i     = 2*pi * validAPhases / meanPer_i;
        cMeanRad   = atan2(mean(sin(angs_i)), mean(cos(angs_i)));
        cMeanPhase = mod(cMeanRad * meanPer_i / (2*pi), meanPer_i);
        cStdPhase  = circ_std(angs_i) * meanPer_i / (2*pi);
    else
        cMeanPhase = NaN;  cStdPhase = NaN;
    end

    % Time-averaged period and total period change across segment
    segDur_i  = segments(i,2) - segments(i,1);
    t0_i      = segments(i,1);
    tmid_i    = segments(i,1) + segDur_i/2;
    vIdx_i    = find(~isnan(segmentPeriods(i,:)) & ~isnan(segmentChirpRates(i,:)));
    Tmean_i   = nan(size(vIdx_i));
    DeltaT_i  = nan(size(vIdx_i));
    for kk = 1:numel(vIdx_i)
        jj_i    = vIdx_i(kk);
        kap_i   = segmentChirpRates(i, jj_i);
        om_mid_i = 2*pi / segmentPeriods(i, jj_i);
        om0_i   = max(om_mid_i - kap_i*(tmid_i - t0_i), 2*pi/35);
        omEnd_i = max(om0_i + kap_i*segDur_i, 2*pi/35);
        DeltaT_i(kk) = 2*pi/omEnd_i - 2*pi/om0_i;
        if abs(kap_i) > 1e-9
            Tmean_i(kk) = (2*pi/(kap_i*segDur_i)) * log(omEnd_i/om0_i);
        else
            Tmean_i(kk) = 2*pi / om0_i;
        end
    end

    segmentResults(i).segmentStart  = segments(i,1);
    segmentResults(i).segmentEnd    = segments(i,2);
    segmentResults(i).meanPeriodAvg = mean(Tmean_i,  'omitnan');  % time-averaged
    segmentResults(i).stdPeriodAvg  = std(Tmean_i,   [],  'omitnan');
    segmentResults(i).meanDeltaT    = mean(DeltaT_i, 'omitnan');  % T_end - T_start
    segmentResults(i).stdDeltaT     = std(DeltaT_i,  [],  'omitnan');
    segmentResults(i).meanChirpRate = mean(validKappas);
    segmentResults(i).stdChirpRate  = std(validKappas);
    segmentResults(i).meanPhase       = cMeanPhase;
    segmentResults(i).stdPhase        = cStdPhase;
    segmentResults(i).meanLambda      = mean(validLambdas);
    segmentResults(i).stdLambda       = std(validLambdas);
    segmentResults(i).coherence       = phaseCoherence(i);
    segmentResults(i).nValidSamples   = numel(validPeriods);
end

% =========================================================================
% Cross-segment damping analysis
% =========================================================================
% Model: log(amplitude_i) = A0 - gamma * t_mid_i  =>  gamma = decay rate
validDampingSamples = false(numSamples, 1);
dampingFile = fopen(strcat(sampleName, '_damping_coefficients.csv'), 'w');
fprintf(dampingFile, 'Sample,DampingCoefficient_per_h\n');

for j = 1:numSamples
    validSegs  = ~isnan(segmentAmplitudes(:,j)) & (segmentAmplitudes(:,j) > 0);
    amplitudes = segmentAmplitudes(validSegs, j);
    times      = segmentMidTimes(validSegs);

    if numel(amplitudes) < 2; continue; end

    lnA = log(amplitudes);
    t0  = times(1);
    try
        coeffs = robustfit(times - t0, lnA);
        dampingCoeffs(j)       = -coeffs(2);
        validDampingSamples(j) = true;
        fprintf(dampingFile, '%s,%.6f\n', samplesList{j}, dampingCoeffs(j));
    catch
        continue;
    end
end
fclose(dampingFile);

validCoeffs = dampingCoeffs(validDampingSamples);
if ~isempty(validCoeffs)
    meanDamping = mean(validCoeffs);
    stdDamping  = std(validCoeffs);
end

% =========================================================================
% Phase drift analysis  (acrophase vs segment midpoint time)
% =========================================================================
phaseDriftFile = fopen(strcat(sampleName, '_phase_drift.csv'), 'w');
fprintf(phaseDriftFile, 'Sample,DriftSlope_h_per_h,ResidualSD_h,RSquared\n');

for j = 1:numSamples
    validSegs = find(validSegmentFlags(j,:));
    if numel(validSegs) < 2; continue; end

    times  = segmentMidTimes(validSegs);
    phases = samplePhaseTrajectories(j, validSegs)';
    times  = times(:);
    phases = phases(:);

    try
        b         = robustfit(times, phases);
        predicted = b(1) + b(2)*times;
        residuals = phases - predicted;
        driftSD   = std(residuals);
        rsq       = 1 - sum(residuals.^2) / ...
                        sum((phases - mean(phases)).^2);
        phaseDriftSlopes(j) = b(2);
        fprintf(phaseDriftFile, '%s,%.6f,%.4f,%.4f\n', ...
            samplesList{j}, phaseDriftSlopes(j), driftSD, rsq);
    catch
        continue;
    end
end
fclose(phaseDriftFile);

validDrifts = phaseDriftSlopes(~isnan(phaseDriftSlopes));
if ~isempty(validDrifts)
    meanPhaseDrift = mean(validDrifts);
    stdPhaseDrift  = std(validDrifts);
end

% =========================================================================
% Phase coherence CSV
% =========================================================================
cohFile = fopen(strcat(sampleName, '_phase_coherence.csv'), 'w');
fprintf(cohFile, 'Time_h,Coherence_1perrad,ValidPhases\n');
for i = 1:numSegments
    if ~isnan(phaseCoherence(i))
        fprintf(cohFile, '%.2f,%.4f,%d\n', ...
            segmentMidTimes(i), phaseCoherence(i), ...
            sum(~isnan(segmentAdjustedPhases(i,:))));
    end
end
fclose(cohFile);

% =========================================================================
% CSV exports
% =========================================================================
% --- Individual periods (at midpoint) ---
fid = fopen(strcat(sampleName, '_individual_periods.csv'), 'w');
fprintf(fid, 'Sample');
for i = 1:numSegments
    fprintf(fid, ',Segment_%d-%dh', segments(i,1), segments(i,2));
end
fprintf(fid, '\n');
for j = 1:numSamples
    fprintf(fid, '%s', samplesList{j});
    for i = 1:numSegments
        fprintf(fid, ',%.6f', segmentPeriods(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- Individual time-averaged periods (T_mean, h) ---
fid = fopen(strcat(sampleName, '_individual_periods_mean.csv'), 'w');
fprintf(fid, 'Sample');
for i = 1:numSegments
    fprintf(fid, ',Segment_%d-%dh', segments(i,1), segments(i,2));
end
fprintf(fid, '\n');
for j = 1:numSamples
    fprintf(fid, '%s', samplesList{j});
    for i = 1:numSegments
        fprintf(fid, ',%.6f', segmentPeriodMeans(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- Individual period change DeltaT (h) ---
fid = fopen(strcat(sampleName, '_individual_deltaT.csv'), 'w');
fprintf(fid, 'Sample');
for i = 1:numSegments
    fprintf(fid, ',Segment_%d-%dh', segments(i,1), segments(i,2));
end
fprintf(fid, '\n');
for j = 1:numSamples
    fprintf(fid, '%s', samplesList{j});
    for i = 1:numSegments
        fprintf(fid, ',%.6f', segmentDeltaTs(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- Individual acrophases (hours) ---
fid = fopen(strcat(sampleName, '_individual_phases.csv'), 'w');
fprintf(fid, 'Sample');
for i = 1:numSegments
    fprintf(fid, ',Segment_%d-%dh', segments(i,1), segments(i,2));
end
fprintf(fid, '\n');
for j = 1:numSamples
    fprintf(fid, '%s', samplesList{j});
    for i = 1:numSegments
        fprintf(fid, ',%.6f', segmentAdjustedPhases(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- Individual chirp rates (rad/h^2) ---
fid = fopen(strcat(sampleName, '_individual_chirprates.csv'), 'w');
fprintf(fid, 'Sample');
for i = 1:numSegments
    fprintf(fid, ',Segment_%d-%dh', segments(i,1), segments(i,2));
end
fprintf(fid, '\n');
for j = 1:numSamples
    fprintf(fid, '%s', samplesList{j});
    for i = 1:numSegments
        fprintf(fid, ',%.8f', segmentChirpRates(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- Individual lambda per fit (h^-1) ---
fid = fopen(strcat(sampleName, '_individual_lambdas.csv'), 'w');
fprintf(fid, 'Sample');
for i = 1:numSegments
    fprintf(fid, ',Segment_%d-%dh', segments(i,1), segments(i,2));
end
fprintf(fid, '\n');
for j = 1:numSamples
    fprintf(fid, '%s', samplesList{j});
    for i = 1:numSegments
        fprintf(fid, ',%.6f', segmentDampings(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- Individual raw peak times ---
fid = fopen(strcat(sampleName, '_individual_peak_times.csv'), 'w');
fprintf(fid, 'Sample');
for segIdx = 1:numSegments
    fprintf(fid, ',Seg%d_Peak1_h,Seg%d_Peak2_h', segIdx, segIdx);
end
fprintf(fid, '\n');
for j = 1:numSamples
    fprintf(fid, '%s', samplesList{j});
    for i = 1:numSegments
        pt = rawFirstTwoPeakTimesPerSegment{i}(j,:);
        if ~isnan(pt(1)); fprintf(fid, ',%.4f', pt(1)); else; fprintf(fid, ','); end
        if ~isnan(pt(2)); fprintf(fid, ',%.4f', pt(2)); else; fprintf(fid, ','); end
    end
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('CSV files written with prefix: %s\n', sampleName);

% =========================================================================
% Console summary
% =========================================================================
fprintf('\n=== %s: Segment Analysis Summary ===\n', sampleName);
fprintf('Phase reference: t = %.1f h  ', phaseOffsetHours);
fprintf('(acrophases relative to this time)\n\n');
fprintf('%-20s  %-18s  %-16s  %-20s\n', ...
    'Segment', 'T_mean (h)', 'dT (h)', 'kappa (rad/h^2)');
fprintf('%s\n', repmat('-', 1, 80));
for i = 1:numSegments
    fprintf('  %4d-%4d h       %.2f+/-%.2f    %.2f+/-%.2f    %.5f+/-%.5f  (n=%d)\n', ...
        segments(i,1), segments(i,2), ...
        segmentResults(i).meanPeriodAvg, segmentResults(i).stdPeriodAvg, ...
        segmentResults(i).meanDeltaT,    segmentResults(i).stdDeltaT,    ...
        segmentResults(i).meanChirpRate, segmentResults(i).stdChirpRate, ...
        segmentResults(i).nValidSamples);
end

fprintf('\nAmplitude damping (cross-segment, h^-1):\n');
if any(validDampingSamples)
    for j = find(validDampingSamples)'
        fprintf('  %-20s  %.4f\n', samplesList{j}, dampingCoeffs(j));
    end
    fprintf('  Mean: %.4f +/- %.4f  (n=%d)\n', ...
        meanDamping, stdDamping, nnz(validDampingSamples));
else
    fprintf('  No valid damping estimates.\n');
end

fprintf('\nPhase coherence (1/rad):\n');
validCoh = phaseCoherence(~isnan(phaseCoherence));
for i = 1:numSegments
    if ~isnan(phaseCoherence(i))
        fprintf('  %d-%d h:  %.4f  (n=%d/%d)\n', ...
            segments(i,1), segments(i,2), phaseCoherence(i), ...
            sum(~isnan(segmentAdjustedPhases(i,:))), numSamples);
    end
end
if ~isempty(validCoh)
    fprintf('  Mean: %.4f +/- %.4f\n', mean(validCoh), std(validCoh));
end

fprintf('\nPhase drift (h/h):\n');
if ~isempty(validDrifts)
    fprintf('  Mean: %.4f +/- %.4f\n', meanPhaseDrift, stdPhaseDrift);
    for j = find(~isnan(phaseDriftSlopes))'
        fprintf('  %-20s  %.4f\n', samplesList{j}, phaseDriftSlopes(j));
    end
else
    fprintf('  No valid phase drift estimates.\n');
end

% =========================================================================
% Assign outputs
% =========================================================================
allPeriods        = segmentPeriods;
allPhases         = segmentPhases;
allAdjustedPhases = segmentAdjustedPhases;
allFittedCurves   = segmentFittedCurves;
allChirpRates     = segmentChirpRates;

end % main function

% =========================================================================
% Local helper: circular standard deviation
% =========================================================================
function s = circ_std(alpha)
% circ_std  Circular SD for a vector of angles (radians).
%   Based on Batschelet (1981): s = sqrt(-2*log(R)), R = |mean(exp(i*alpha))|
    R = abs(mean(exp(1i * alpha(:))));
    s = sqrt(-2 * log(max(R, eps)));
end

% =========================================================================
% Utility: annotate figure with group-level summary statistics
% =========================================================================
function addGlobalStatsLabel(meanDamp, stdDamp, meanDrift, stdDrift, ...
                              meanCoh, stdCoh, meanKappa, stdKappa)
% addGlobalStatsLabel  Add group-level annotation box to the current figure.
%
%   Usage:
%     addGlobalStatsLabel(meanDamping, stdDamping, ...
%                         meanPhaseDrift, stdPhaseDrift, ...
%                         mean(validCoh), std(validCoh), ...
%                         mean(kappas), std(kappas));
    txt = { ...
        sprintf('Mean lambda:   %.3f +/- %.3f h^-1',   meanDamp,  stdDamp)
        sprintf('Mean kappa:    %.5f +/- %.5f rad/h^2', meanKappa, stdKappa)
        sprintf('Phase drift:   %.3f +/- %.3f h/h',    meanDrift, stdDrift)
        sprintf('Coherence:     %.3f +/- %.3f 1/rad',  meanCoh,   stdCoh)
    };
    annotation(gcf, 'textbox', [0.60 0.82 0.36 0.15], ...
        'String',             txt,            ...
        'FitBoxToText',       'on',           ...
        'BackgroundColor',    [1 1 1 0.85],   ...
        'EdgeColor',          [0.4 0.4 0.4],  ...
        'LineWidth',          1,              ...
        'FontSize',           10,             ...
        'HorizontalAlignment','left');
end

% =========================================================================
% Local callback - redraw individual-fit axes for the selected sample.
% Axes found by Tag so the callback works after reopening the .fig file.
% =========================================================================
function ifc_updatePlot(hList)
    hFig = ancestor(hList, 'figure');
    ud   = get(hFig, 'UserData');
    j    = get(hList, 'Value');
    set(hFig, 'Name', sprintf('Individual Fitted Curves  -  %s', ud.samplesList{j}));

    fullT = ud.fullTime;
    fullD = ud.fullData{j};

    for i = 1:ud.numSegments
        ax = findobj(hFig, 'Type', 'axes', 'Tag', sprintf('SegmentAxes_%d', i));
        if isempty(ax); continue; end

        legend(ax, 'off');
        cla(ax);
        hold(ax, 'on');

        % Full time course - thin grey line
        plot(ax, fullT, fullD, '-', 'Color', [0.65 0.65 0.65], ...
            'LineWidth', 0.8, 'DisplayName', 'Full time course');

        % Fitted curve - thick coloured line
        fitT = ud.commonTimePerSegment{i};
        fitC = ud.segmentFittedCurves{i}(j, :);
        rawT = ud.segmentRawTimeStore{i, j};
        if ~isempty(rawT) && ~all(isnan(fitC))
            if ~isempty(ud.segmentFitLabelStore{i,j})
                fitLbl = ud.segmentFitLabelStore{i,j};
            else
                fitLbl = 'Fitted curve';
            end
            plot(ax, fitT, fitC, '-', 'Color', ud.colors(j,:), ...
                'LineWidth', 2.5, 'DisplayName', fitLbl);
        end

        title(ax, sprintf('Segment  %d-%d h', ud.segments(i,1), ud.segments(i,2)), ...
            'FontSize', 9);
        xlabel(ax, 'Time (h)', 'FontSize', 8);
        ylabel(ax, 'Luminescence (cps)', 'FontSize', 8);
        xlim(ax, [min(fullT), max(fullT)]);
ax.XTick = 0:24:max(fullT);
ax.XMinorTick = 'on';
ax.XAxis.MinorTickValues = 12:24:max(fullT)-12;
        lgd_i = legend(ax, 'show');
        lgd_i.Location = 'northeast'; lgd_i.FontSize = 7;
        lgd_i.Box = 'off'; lgd_i.Interpreter = 'none';
        grid(ax, 'on'); box(ax, 'on'); ax.FontSize = 8;
        hold(ax, 'off');
    end
end
