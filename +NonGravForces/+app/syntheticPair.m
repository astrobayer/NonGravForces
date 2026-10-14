function [obs, ref, truth] = syntheticPair(options)
%SYNTHETICPAIR Truth generator: a 10 Hz model arc and its CRN-filtered 1 Hz observation, with every planted value.
%
%   [obs, ref, truth] = NonGravForces.app.syntheticPair()
%   [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 7, 'Scale', S, 'Quadratic', k2)
%
%   The ground truth the calibration tests share (N-02, CONCEPT ORC-7; the
%   v3 plan's WP-21). It builds, on the 10 Hz grid of ACT1A:
%
%       model(t)  orbital-band accelerations: along-track drag-like -2e-7 m/s^2
%                 modulated at 1/rev and by a slow solar-flux proxy (27-day
%                 period), cross-track 1/rev + 2/rev at 1e-8, radial at 5e-8
%       raw(t)    = S model + b(t) + K2 .* (model.^2 - mean) + noise + spikes,
%                 with b(t) = Bias + BiasDrift * (t / 1 day), noise coloured
%                 per axis (ASD NoiseAsd * sqrt(1 + (NoiseCorner/f)^2), the
%                 sensor.Electrostatic form) from timeSeriesPreparation's
%                 synthTimeSeries, spikes 1 s boxcars standing for thruster
%                 firings, and one gap of NaN
%
%   and passes raw through the REAL Level-1B chain: timeSeriesPreparation.
%   crnKernel("ACC") applied by crnApply and evaluated at the whole seconds
%   only, as ACT1B is formed. The observation therefore carries everything a
%   calibration meets in real data -- the filter, its 70.3 s edges, the
%   windows a gap removes, filtered impulses -- with every planted number
%   returned in `truth`. The reference is the model itself on the 10 Hz grid,
%   because the CRN stage needs the reference on its raw grid (CONCEPT 6.1).
%
%   ALL RANDOMNESS COMES FROM 'Seed': the spike epochs and axes from one
%   RandStream, the noise of axis m from synthTimeSeries with seed Seed*10+m,
%   which restores the global generator on exit. The same Seed gives the
%   same arcs bit for bit.
%
% INPUT ARGUMENTS
%   None.
%
% NAME-VALUE ARGUMENTS
%   Name            Size  Type    Units        Description
%   ----            ----  ----    -----        -----------
%   Days            1x1   double  day          Arc length, > 0; default 1
%   Start           1x1   int64   ns           First epoch, GPS ns past 2000-01-01 12:00:00;
%                                              default 2020-01-01 00:00:18 GPS
%   Seed            1x1   double  -            Seed of every random draw; default 1
%   Scale           3x3   double  -            True scale matrix S (a 3x1 is its diagonal);
%                                              default diag([0.96 0.98 0.95])
%   Bias            3x1   double  m/s^2        True bias at Start; default [2e-7; -1e-7; 5e-8]
%   BiasDrift       3x1   double  m/s^2/day    True linear bias drift; default [1e-9; 0; -5e-10]
%   Quadratic       3x1   double  s^2/m        True K2 per axis; default zeros
%   NoiseAsd        3x1   double  m/s^2/rtHz   White noise floor per axis; default [1e-10; 1e-9; 1e-10]
%   NoiseCorner     1x1   double  Hz           1/f corner of the noise; default 5e-3
%   Spikes          1x1   double  -            Number of planted 1 s spikes; default 3
%   SpikeAmplitude  1x1   double  m/s^2        Their amplitude; default 5e-8
%   Gap_s           1x1   double  s            Length of the one planted gap, 0 for none; default 30
%   OrbitPeriod_s   1x1   double  s            Orbital period; default 5640 (~490 km)
%
% OUTPUT ARGUMENTS
%   Name     Size  Type                        Units  Description
%   ----     ----  ----                        -----  -----------
%   obs      1x1   timeSeriesPreparation.Arc   m/s^2  CRN-filtered observation at 1 Hz, columns
%                                                     lin_accl_x/y/z; NaN where the filter window
%                                                     met an edge or the gap
%   ref      1x1   timeSeriesPreparation.Arc   m/s^2  The model at 10 Hz, same column names
%   truth    1x1   struct                      -      .scale .bias0 .biasDrift_per_day
%                                                     .quadratic .bias (3xK1 at obs.t) .noise
%                                                     (.asd .corner .seeds) .spikes (table: t
%                                                     int64, axis, amplitude) .gap (int64 [first
%                                                     last] on the 10 Hz grid, or empty)
%                                                     .kernel (crnKernel struct) .orbitPeriod_s
%                                                     .seed .provenance
%
% ACCURACY
%   Exact by construction: the planted values ARE the truth. The noise is a
%   single Gaussian realisation, so its sample variance scatters around the
%   target like any estimate (synthTimeSeries' own accuracy). The model is a
%   plausible orbital-band signal, not a force model; use NonGravForces.Model
%   for physics.
%
% ERRORS
%   NonGravForces:syntheticPair:depMissing    timeSeriesPreparation is not on the path
%   NonGravForces:syntheticPair:inputShape    Scale, Bias, BiasDrift, Quadratic or NoiseAsd of
%                                             the wrong size, or not finite
%   NonGravForces:syntheticPair:inputRange    the arc is too short for the CRN window plus the
%                                             spikes and the gap, or the spikes find no room
%                                             clear of each other
%   timeSeriesPreparation:*                   from crnKernel, crnApply, synthTimeSeries, Arc
%
% REFERENCES
%   Wu S.-C., Kruizinga G., Bertiger W. (2006) GRACE L1B ATBD, JPL D-27672,
%   App. B (the CRN chain); timeSeriesPreparation docs/CONCEPT.md ORC-7 and
%   6.2; the v3 plan's WP-21 (timeSeriesPreparation docs/inputs).
%
% EXAMPLE
%   [obs, ref, truth] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 3);
%   size(obs.obs), truth.scale
%
% See also NonGravForces.app.Calibration, NonGravForces.app.calibrate, timeSeriesPreparation.crnApply.
%
% KEYWORDS synthetic truth generator calibration accelerometer CRN test
% Created 2026-10-14 by Claude (Anthropic) for Matthias Weigelt.
arguments
    options.Days (1,1) double {mustBePositive, mustBeFinite} = 1
    options.Start (1,1) int64 = int64(7305 * 86400 - 43200 + 18) * int64(1e9)
    options.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    options.Scale double = diag([0.96 0.98 0.95])
    options.Bias double = [2e-7; -1e-7; 5e-8]
    options.BiasDrift double = [1e-9; 0; -5e-10]
    options.Quadratic double = zeros(3, 1)
    options.NoiseAsd double = [1e-10; 1e-9; 1e-10]
    options.NoiseCorner (1,1) double {mustBePositive, mustBeFinite} = 5e-3
    options.Spikes (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    options.SpikeAmplitude (1,1) double {mustBeFinite} = 5e-8
    options.Gap_s (1,1) double {mustBeNonnegative, mustBeFinite} = 30
    options.OrbitPeriod_s (1,1) double {mustBePositive, mustBeFinite} = 5640
end
unit = "syntheticPair";
if isempty(which('timeSeriesPreparation.crnApply'))
    NonGravForces.fail(unit, "depMissing", "the observation goes through timeSeriesPreparation.crnApply and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation");
end
S = options.Scale;
if isvector(S) && numel(S) == 3, S = diag(S(:)); end
vec3 = {options.Bias, options.BiasDrift, options.Quadratic, options.NoiseAsd};
if ~isequal(size(S), [3 3]) || any(cellfun(@numel, vec3) ~= 3) || ~all(isfinite([S(:); cell2mat(cellfun(@(v) v(:), vec3(:), "UniformOutput", false))]))
    NonGravForces.fail(unit, "inputShape", "Scale must be 3x3 (or 3x1) and Bias, BiasDrift, Quadratic, NoiseAsd 3x1, all finite -> check the sizes");
end
b0 = options.Bias(:); drift = options.BiasDrift(:); k2 = options.Quadratic(:); asd = options.NoiseAsd(:);

kernel = timeSeriesPreparation.crnKernel("ACC");
fs = kernel.fs;
n = round(options.Days * 86400 * fs);
edge = (kernel.Nf - 1) / 2;
gapN = round(options.Gap_s * fs);
if n < 2 * kernel.Nf + gapN + 4 * fs * max(options.Spikes, 1)
    NonGravForces.fail(unit, "inputRange", "%g day(s) is too short for the %d-sample CRN window plus the spikes and the gap -> lengthen 'Days'", options.Days, kernel.Nf);
end
stepNs = int64(round(1e9 / fs));
t = options.Start + int64(0:n - 1)' * stepNs;
ts = double(0:n - 1)' / fs;                                    % seconds since Start

% The model: an orbital-band signal, along-track drag-like with a solar-flux proxy.
u = 2 * pi * ts / options.OrbitPeriod_s;
flux = 1 + 0.15 * sin(2 * pi * ts / (27 * 86400) + 0.4);
model = [-2e-7 * flux .* (1 + 0.4 * cos(u)), ...
         1e-8 * (sin(u) + 0.3 * cos(2 * u)), ...
         5e-8 * (cos(u + 1) - 0.2 * sin(2 * u))];

% The observation on the raw grid.
biasRaw = b0' + drift' .* (ts / 86400);
quadRaw = k2' .* (model.^2 - mean(model.^2, 1));
fGrid = logspace(log10(fs / n), log10(fs / 2), 400)';
noise = zeros(n, 3);
seeds = options.Seed * 10 + (1:3);
for m = 1:3
    psd = asd(m)^2 * (1 + (options.NoiseCorner ./ fGrid).^2);
    noise(:, m) = timeSeriesPreparation.spectral.synthTimeSeries("psd", fGrid, psd, fs, n, Seed = seeds(m));
end
raw = model * S.' + biasRaw + quadRaw + noise;

% Spikes and the gap, placed by one stream, clear of the edges and of each other.
rs = RandStream("twister", "Seed", options.Seed);
gap = int64([]);
if gapN > 0
    g0 = edge + 1 + floor((n - 2 * edge - gapN) / 2);       % the gap sits mid-arc
    raw(g0:g0 + gapN - 1, :) = NaN;
    gap = t([g0, g0 + gapN - 1]).';
end
spikeT = zeros(options.Spikes, 1, 'int64'); spikeAxis = zeros(options.Spikes, 1);
taken = false(n, 1);
if gapN > 0, taken(max(1, g0 - kernel.Nf):min(n, g0 + gapN + kernel.Nf)) = true; end
for s = 1:options.Spikes
    placed = false;
    for attempt = 1:1000
        i0 = edge + 1 + randi(rs, n - 2 * edge - fs);
        win = max(1, i0 - kernel.Nf):min(n, i0 + fs + kernel.Nf);
        if ~any(taken(win)), placed = true; break; end
    end
    if ~placed
        NonGravForces.fail(unit, "inputRange", "no room for spike %d of %d clear of the others and of the gap -> lengthen 'Days' or plant fewer 'Spikes'", s, options.Spikes);
    end
    taken(win) = true;
    ax = randi(rs, 3);
    raw(i0:i0 + fs - 1, ax) = raw(i0:i0 + fs - 1, ax) + options.SpikeAmplitude;
    spikeT(s) = t(i0); spikeAxis(s) = ax;
end

names = {'lin_accl_x', 'lin_accl_y', 'lin_accl_z'};
day = string(datetime(2000, 1, 1, 12, 0, 0) + seconds(double(options.Start) / 1e9), 'yyyy-MM-dd');
rawArc = timeSeriesPreparation.Arc(t, raw, localAux(names, "ACT1A_" + day + "_C_04.txt", "ACT1A"));
out1Hz = t(1:round(fs):end);
obs = timeSeriesPreparation.crnApply(rawArc, kernel, OutputEpochs = out1Hz);
ref = timeSeriesPreparation.Arc(t, model, localAux(names, "MODEL_" + day + "_C_04.txt", "MODEL"));

t1 = double(obs.t - options.Start) / 1e9;
truth = struct('scale', S, 'bias0', b0, 'biasDrift_per_day', drift, 'quadratic', k2, ...
    'bias', b0 + drift .* (t1' / 86400), ...
    'noise', struct('asd', asd, 'corner', options.NoiseCorner, 'seeds', seeds), ...
    'spikes', table(spikeT, spikeAxis, repmat(options.SpikeAmplitude, options.Spikes, 1), ...
        'VariableNames', {'t', 'axis', 'amplitude'}), ...
    'gap', gap, 'kernel', kernel, 'orbitPeriod_s', options.OrbitPeriod_s, 'seed', options.Seed, ...
    'provenance', NonGravForces.provenance("app.syntheticPair", 'Options', options));
end

function aux = localAux(names, fileName, product)
%LOCALAUX The reader-shaped metadata timeSeriesPreparation.Arc builds from (its field names).
units = struct();
for k = 1:numel(names), units.(names{k}) = 'm/s2'; end
aux = struct("obsNames", {names}, "variableNames", {names}, "variableFieldNames", {names}, ...
    "variableUnits", units, "fileName", fileName, "product", product);
end
