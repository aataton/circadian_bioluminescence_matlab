function [amplitude, period, phase, adjustedPhase, offset, damping, chirpRate, fitType] = ...
        fitCosineCurve(time, samples, periodRange, phaseOffsetHours)
% fitCosineCurve  Fit a chirped damped cosine to circadian time-series data.
%
%   Uses FFT-guided nonlinear least-squares (NLLS) to fit:
%
%     y(t) = A * exp(-lambda*u) * cos(omega0*u + 0.5*kappa*u^2 + phi0) + B
%
%   where  u = t - t0,  t0 = time(1)  (segment start).
%
%   Parameters
%   ----------
%   A      : amplitude at t0  (>= 0).
%   lambda : damping rate (h^-1).  lambda>0 => decay, lambda<0 => growth, lambda=0 => flat.
%   omega0 : angular frequency at t0 (rad/h).  period0 = 2*pi/omega0.
%   kappa  : chirp rate (rad/h^2).
%              kappa = 0  => constant period (standard cosine).
%              kappa > 0  => frequency increases => period shortens.
%              kappa < 0  => frequency decreases => period lengthens.
%            Instantaneous period at time t:
%              T(t) = 2*pi / (omega0 + kappa*(t - t0))
%   phi0   : initial phase at t0 (rad).
%   B      : vertical offset (mesor).
%
%   Inputs
%   ------
%   time             - Column/row vector of time points (hours).
%   samples          - Column/row vector of data values (same length as time).
%   periodRange      - [minPeriod, maxPeriod] in hours. Default [18, 30].
%   phaseOffsetHours - Reference time (h) for acrophase reporting:
%                        adjustedPhase = mod(t_peak - phaseOffsetHours, T_mid)
%                      where t_peak is the nearest peak to the segment
%                      midpoint and T_mid is the period at the midpoint.
%                      Use 0 for ZT0 / CT0.
%
%   Outputs
%   -------
%   amplitude     - Fitted amplitude at t0 (>= 0).
%   period        - Instantaneous period at segment midpoint (hours).
%                   period_at_start = 2*pi / omega0.
%                   Use chirpRate to compute T at any other time.
%   phase         - Raw fitted initial phase phi0 (rad, in [0, 2*pi)).
%                   Reconstruction:
%                     A*exp(-lambda*u).*cos(omega0*u+0.5*kappa*u^2+phase)+B
%                   where omega0 = 2*pi/periodAtStart (see below) and
%                   u = time - time(1).
%   adjustedPhase - Acrophase relative to phaseOffsetHours, in [0,T_mid) h.
%                   Computed numerically on a 10000-point dense grid.
%   offset        - Fitted mesor B.
%   damping       - Fitted lambda (h^-1). Negative => growing amplitude.
%   chirpRate     - Fitted kappa (rad/h^2).
%                   Period change rate at midpoint (h/h):
%                     dT/dt_mid = -kappa * T_mid^2 / (2*pi)
%   fitType       - String 'chirped_damped_cosine'.
%
%   Note: periodAtStart = 2*pi / (omega_mid + kappa*(t_mid - t0))
%         where omega_mid = 2*pi / period.  Store chirpRate alongside
%         period to fully characterise the fit.
%
%   Returns all-NaN outputs (with a warning) when the fit fails or when
%   there are fewer than 5 data points.

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
time    = time(:);
samples = samples(:);
if numel(time) ~= numel(samples)
    error('fitCosineCurve: time and samples must have the same length.');
end
if numel(time) < 5
    warning('fitCosineCurve: Too few data points (need >= 5). Returning NaNs.');
    [amplitude, period, phase, adjustedPhase, offset, damping, chirpRate] = deal(NaN);
    fitType = 'chirped_damped_cosine';
    return
end
periodRange = sort(periodRange(:)');

% -------------------------------------------------------------------------
% FFT-based initial frequency estimate (mean-subtracted for cleaner spectrum)
% -------------------------------------------------------------------------
Fs   = 1 / mean(diff(time));
NFFT = 2^nextpow2(numel(samples));

Y  = fft(samples - mean(samples), NFFT);
P2 = abs(Y / NFFT);
P1 = P2(1:NFFT/2+1);
P1(2:end-1) = 2 * P1(2:end-1);
f  = Fs * (0:NFFT/2) / NFFT;

freqRange = sort(1 ./ periodRange);          % [1/Tmax, 1/Tmin]
validIdx  = (f >= freqRange(1)) & (f <= freqRange(2));

if any(validIdx)
    [~, iRel]      = max(P1(validIdx));
    allValid       = find(validIdx);
    mainIdx        = allValid(iRel);
    dominantFreqHz = f(mainIdx);
    phase_init     = angle(Y(mainIdx));
    amp_init       = max(P1(validIdx));
else
    dominantFreqHz = 1 / mean(periodRange);
    [~, mainIdx]   = min(abs(f - dominantFreqHz));
    phase_init     = angle(Y(mainIdx));
    amp_init       = max(P1);
    warning('fitCosineCurve: No FFT peak found in [%.1f, %.1f] h. Using fallback.', ...
            periodRange(1), periodRange(2));
end

omega_init = 2*pi * dominantFreqHz;
t0_val     = time(1);
t_mid_val  = mean(time);

% -------------------------------------------------------------------------
% Chirped damped cosine NLLS fit
% -------------------------------------------------------------------------
% fittype expression with u = (x - t0):
%
%   a * exp(-e*(x-t0)) * cos(b*(x-t0) + 0.5*g*(x-t0)^2 + c) + d
%
% Coefficients (alphabetical order -- how MATLAB fittype assigns them):
%   a  amplitude at t0      bounds: (-Inf, Inf)   start: amp_init
%   b  omega0   rad/h       bounds: [omega_min, omega_max]
%   c  phi0     rad         bounds: (-Inf, Inf)   start: phase_init
%   d  offset               bounds: (-Inf, Inf)   start: mean(samples)
%   e  lambda   h^-1        bounds: [-0.02, 0.2]  start: 0.01  (negative = growing amplitude)
%   g  kappa    rad/h^2     bounds: [-kappa_b, kappa_b]  start: 0
%
% t0 is a fixed "problem" constant equal to time(1).

omega_min = 2*pi / periodRange(2);
omega_max = 2*pi / periodRange(1);

% Kappa bound: period must stay within periodRange over the segment.
segDur      = time(end) - t0_val;
kappa_bound = (omega_max - omega_min) / max(segDur, 1);
kappa_bound = min(kappa_bound, 0.03);   % hard cap: 0.03 rad/h^2 (relaxed for better chirp detection)

initialGuess = [amp_init,  omega_init, phase_init, mean(samples), 0.01,  0          ];
lowerBounds  = [-Inf,      omega_min,  -Inf,        -Inf,        -0.02, -kappa_bound];
upperBounds  = [ Inf,      omega_max,   Inf,          Inf,         0.2,   kappa_bound];

fitOptions = fitoptions( ...
    'Method',      'NonlinearLeastSquares', ...
    'StartPoint',  initialGuess, ...
    'Lower',       lowerBounds,  ...
    'Upper',       upperBounds,  ...
    'MaxIter',     5000,         ...
    'MaxFunEvals', 10000,        ...
    'TolFun',      1e-9,         ...
    'TolX',        1e-9);

chirpedModel = fittype( ...
    'a * exp(-e*(x-t0)) * cos(b*(x-t0) + 0.5*g*(x-t0)^2 + c) + d', ...
    'independent', 'x', ...
    'problem',     't0');

% -------------------------------------------------------------------------
% Multistart over kappa: try kappa = 0, +half-bound, -half-bound.
% Keep the result with the highest R^2 to avoid the local minimum at
% kappa = 0 that the single-start strategy sometimes converges to.
% -------------------------------------------------------------------------
kappaStarts  = [0,  kappa_bound/2, -kappa_bound/2];
bestR2       = -Inf;
fitResult    = [];

for ks = kappaStarts
    startK = initialGuess;
    startK(6) = ks;          % override kappa starting value
    fitOptionsK = fitoptions( ...
        'Method',      'NonlinearLeastSquares', ...
        'StartPoint',  startK, ...
        'Lower',       lowerBounds,  ...
        'Upper',       upperBounds,  ...
        'MaxIter',     5000,         ...
        'MaxFunEvals', 10000,        ...
        'TolFun',      1e-9,         ...
        'TolX',        1e-9);
    try
        fr   = fit(time, samples, chirpedModel, fitOptionsK, 'problem', {t0_val});
        yhat = feval(fr, time);
        ss_res = sum((samples - yhat).^2);
        ss_tot = sum((samples - mean(samples)).^2);
        r2k  = 1 - ss_res / ss_tot;
        if r2k > bestR2
            bestR2    = r2k;
            fitResult = fr;
        end
    catch
        % this starting point failed – try the next one
    end
end

if isempty(fitResult)
    warning('fitCosineCurve: all NLLS starting points failed. Returning NaNs.');
    [amplitude, period, phase, adjustedPhase, offset, damping, chirpRate] = deal(NaN);
    fitType = 'chirped_damped_cosine';
    return
end

% -------------------------------------------------------------------------
% Extract parameters
% -------------------------------------------------------------------------
a      = fitResult.a;
omega0 = fitResult.b;
c      = fitResult.c;
offset = fitResult.d;
lambda = fitResult.e;
kappa  = fitResult.g;

% Enforce positive amplitude (flip sign => add pi to phase)
if a < 0
    a = -a;
    c = c + pi;
end
amplitude = a;
damping   = lambda;   % negative lambda => growing amplitude (allowed)
chirpRate = kappa;

% Raw initial phase in [0, 2*pi)
phase = mod(c, 2*pi);

% Instantaneous period at segment midpoint (most representative single value)
omega_mid = omega0 + kappa * (t_mid_val - t0_val);
omega_mid = max(omega_mid, omega_min);   % safety clamp
period    = 2*pi / omega_mid;


% -------------------------------------------------------------------------
% Acrophase via numerical peak detection on a dense grid
% -------------------------------------------------------------------------
% No closed-form solution exists for the chirped case, so evaluate the
% model on a dense grid and choose the peak based on phaseOffsetHours.

% Dense time grid over the segment
tdense  = linspace(time(1), time(end), 10000);
udense  = tdense - t0_val;

% Reconstruct fitted curve on dense grid
curvedense = amplitude .* exp(-damping .* udense) .* ...
             cos(omega0 .* udense + 0.5 .* chirpRate .* udense.^2 + phase) + offset;

% Find all peaks of the fitted curve
[~, peakLocs] = findpeaks(curvedense);

if isempty(peakLocs)
    % Fallback: no peaks detected at all, use global maximum
    [~, maxIdx] = max(curvedense);
    tpeak = tdense(maxIdx);
else
    peakTimes = tdense(peakLocs);

    % Try primary rule: first peak at or after phaseOffsetHours
    afterIdx = find(peakTimes >= phaseOffsetHours, 1, 'first');

    if ~isempty(afterIdx)
        % Primary: first peak after offset
        tpeak = peakTimes(afterIdx);
    else
        % Fallback: no peak after offset in this segment,
        % revert to original behavior: peak closest to segment midpoint
        [~, closestIdx] = min(abs(peakTimes - t_mid_val));
        tpeak = peakTimes(closestIdx);
    end
end

% Wrap into one cycle defined by Tmid (=period)
adjustedPhase = mod(tpeak - phaseOffsetHours, period);



fitType = 'chirped_damped_cosine';
end
