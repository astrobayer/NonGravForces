function S = densityScale(files, sw, options)
%DENSITYSCALE Observed daily scale on a model density from Swarm's POD-only densities.
%
%   S = NonGravForces.env.densityScale(files, sw)
%   S = NonGravForces.env.densityScale(files, sw, 'Model', "msis00", 'Stride', 10)
%   f = S.at(epoch, h_m)                        % 1xK factor on the model density
%
%   Turns ESA's Swarm thermosphere densities (io.readSwarmDensity) into a factor on a
%   model density that a Model applies to its drag (Options drag.densityScale "swarm").
%   For each Swarm satellite and UTC day,
%
%       f_s(day) = sum(rho_Swarm) / sum(rho_model)
%
%   over the satellite's valid records of the day (every 'Stride'-th one), with
%   rho_model the named env.density backend at Swarm's own positions and epochs and
%   the same space weather the Model uses. A ratio of sums, not a mean of ratios:
%   the drag of a day is a sum of densities, and a ratio of sums is what scales it.
%
%   Swarm flies in two heights: A and C side by side (about 445 km in 2020), B above
%   (about 515 km). S.at(epoch, h) interpolates log f LINEARLY IN HEIGHT between the
%   lower pair (the geometric mean of A's and C's factors, at their mean height; one of
%   them alone when only one is given) and B, per day: log f = (1 - w) log f_low +
%   w log f_B, w = (h - h_low) / (h_B - h_low). A model error that varies like a scale
%   height is exactly log-linear in height; GRACE-FO lies between the two. Outside the
%   bracket the same line is extrapolated, and S.at refuses w outside [-1, 2] rather
%   than extrapolate further than the bracket is wide.
%
%   WHAT IT DOES NOT DO: Swarm's local-time planes are not GRACE-FO's, and a daily
%   factor removes only the part of the model's error that is common to all local
%   times of a day. On 2020-01 GRACE-C that left k_drag 0.90 where an independent model
%   (DTM2020) reached 0.96 (HANDOVER O-19).
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      Swarm DNSxPOD HAPI CSV files (io.readSwarmDensity), any
%                                 months, of B and of A and/or C. Files of one satellite are
%                                 read together.
%   sw       1x1   struct  -      io.spaceWeather result (celestrak family), the drivers of
%                                 the model density.
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   Model    1x1   string  -      env.density backend the factor is FOR: "msis00" (default),
%                                 "msis21" or "dtm2020". A Model refuses a scale formed for
%                                 another density model than its drag.density.
%   Stride   1x1   double  -      Every Stride-th valid 30 s record enters the sums (default
%                                 10, i.e. 5 min; the product's own resolution is ~20 min).
%   MinRecords 1x1 double  -      Fewer used records in a day make that day's factor NaN, which
%                                 S.at refuses (default 100, about a third of a day at Stride 10).
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   S     1x1   struct  -      .kind "swarm", .day (1xD UTC MJD of each day), .satellite (1xP "A"|"B"|"C"),
%                              .factor (PxD), .height_m (PxD mean geodetic height of the used
%                              records), .nUsed (PxD), .low (1xD) and .high (1xD) bracket factors,
%                              .hLow_m .hHigh_m (1xD), .model, .stride, .files, .provenance, and
%                              .at (function handle: epoch struct (any scale), 1xK heights [m]
%                              -> 1xK factors)
%
% ACCURACY
%   The factor carries Swarm's density error (stated as 30 % of the model variance by the
%   product specification) and the local-time mismatch named above; the sums are exact.
%   The log-linear height interpolation is exact for an error that scales like
%   exp(-h/H); across the 70 km bracket the measured factors of A/C and B differed by 10 %
%   on 2020-01. The day of an epoch is its UTC day.
%
% ERRORS
%   NonGravForces:densityScale:inputSatellites  no B, or neither A nor C, among the files
%   NonGravForces:densityScale:inputModel       Model is not a density backend with drivers
%   NonGravForces:densityScale:dataDay          S.at at an epoch whose day has no factor (no file
%                                               covers it, or fewer than MinRecords records)
%   NonGravForces:densityScale:inputRange       S.at at a height further outside the bracket than
%                                               its width, or h not 1xK
%   NonGravForces:readSwarmDensity:*, NonGravForces:density:*   propagated unwrapped
%
% REFERENCES
%   ESA Swarm Level 2 Product Specification, SW_DNSxPOD_2_; van den IJssel et al. (2020),
%   Adv. Space Res. 65, 1758-1771 (Swarm GPS-derived densities); HANDOVER O-19 (the month
%   that motivated it: rho_Swarm / rho_msis00 ~ 0.48 at GRACE-FO's height in 2020-01).
%
% EXAMPLE
%   root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
%   f = fullfile(root, 'data', 'samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv");
%   sw = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
%   S = NonGravForces.env.densityScale(f, sw);
%   S.at(referenceSystem.mjdToEpoch(58849.5, "UTC"), 490e3)
%
% See also NonGravForces.io.readSwarmDensity, NonGravForces.env.density, NonGravForces.Model, NonGravForces.Options.
%
% KEYWORDS density scale Swarm observed correction drag thermosphere
% Created 2026-10-27 by Claude (Anthropic) for Matthias Weigelt.
arguments
    files (1,:) string {mustBeNonempty}
    sw (1,1) struct
    options.Model (1,1) string = "msis00"
    options.Stride (1,1) double {mustBeInteger, mustBePositive} = 10
    options.MinRecords (1,1) double {mustBeInteger, mustBePositive} = 100
end
unit = "densityScale";
if ~any(options.Model == ["msis00" "msis21" "dtm2020"])
    NonGravForces.fail(unit, "inputModel", "Model %s is not a density backend driven by space weather -> msis00, msis21 or dtm2020", options.Model);
end
% read every file once, group by satellite
sat = strings(1, numel(files)); data = cell(1, numel(files));
for k = 1:numel(files)
    [ep, lla, rho, aux] = NonGravForces.io.readSwarmDensity(files(k));
    sat(k) = aux.satellite;
    data{k} = struct('day', ep.day(:)', 'sec', ep.sec(:)', 'lla', lla, 'rho', rho, 'valid', aux.valid);
end
P = unique(sat, 'stable');
if ~any(P == "B") || ~any(P == "A" | P == "C")
    NonGravForces.fail(unit, "inputSatellites", "the files hold Swarm %s; the height bracket needs B and A or C -> pass the DNSxPOD files of B and of A and/or C", strjoin(P, ", "));
end
P = intersect(["A" "B" "C"], P, 'stable');
merged = cell(1, numel(P));
for j = 1:numel(P)
    parts = [data{sat == P(j)}];
    m.day = [parts.day]; m.sec = [parts.sec]; m.lla = [parts.lla]; m.rho = [parts.rho]; m.valid = [parts.valid];
    [~, iu] = unique([m.day' m.sec'], 'rows');   % a day present in two files counts once
    m.day = m.day(iu); m.sec = m.sec(iu); m.lla = m.lla(:, iu); m.rho = m.rho(iu); m.valid = m.valid(iu);
    merged{j} = m;
end
allDays = cellfun(@(m) m.day, merged, 'UniformOutput', false);
days = unique([allDays{:}]);
D = numel(days);
factor = nan(numel(P), D); height = nan(numel(P), D); nUsed = zeros(numel(P), D);
for j = 1:numel(P)
    m = merged{j};
    for d = 1:D
        k = find(m.day == days(d) & m.valid);
        k = k(1:options.Stride:end);
        nUsed(j, d) = numel(k);
        if numel(k) < options.MinRecords, continue, end
        e = referenceSystem.daySecondToEpoch([m.day(k)' m.sec(k)'], "UTC");
        rhoM = NonGravForces.env.density(options.Model, e, m.lla(:, k), sw);
        factor(j, d) = sum(m.rho(k)) / sum(rhoM);
        height(j, d) = mean(m.lla(3, k));
    end
end
isLow = P == "A" | P == "C";
low = exp(mean(log(factor(isLow, :)), 1, 'omitnan'));   % geometric mean of A and C, or the one with a factor
hLow = mean(height(isLow, :), 1, 'omitnan');
high = factor(P == "B", :); hHigh = height(P == "B", :);
S.kind = "swarm";
S.day = days;
S.satellite = P;
S.factor = factor;
S.height_m = height;
S.nUsed = nUsed;
S.low = low; S.high = high; S.hLow_m = hLow; S.hHigh_m = hHigh;
S.model = options.Model;
S.stride = options.Stride;
S.files = files;
S.at = @(epoch, h) localAt(unit, days, low, high, hLow, hHigh, epoch, h);
S.provenance = NonGravForces.provenance("env.densityScale", 'Files', files, 'Options', options);
end

function f = localAt(unit, days, low, high, hLow, hHigh, epoch, h)
%LOCALAT The factor at each epoch's UTC day and height: log-linear in height between the bracket.
utc = referenceSystem.convertTimeScale(epoch, "UTC");
day = reshape(utc.day, 1, []);
if ~isequal(size(h), size(day))
    NonGravForces.fail(unit, "inputRange", "h must be 1x%d, one height per epoch -> pass lla(3, :)", numel(day));
end
[ok, d] = ismember(day, days);
d(~ok) = 1;                                          % a placeholder index; those epochs are refused below
bad = ~ok | isnan(low(d)) | isnan(high(d));
if any(bad)
    bad = unique(day(bad));
    NonGravForces.fail(unit, "dataDay", "no Swarm density factor for UTC day(s) MJD %s -> fetch the Swarm DNSxPOD files of those days (setupNonGravForces('FetchPool'))", mat2str(bad));
end
w = (h - hLow(d)) ./ (hHigh(d) - hLow(d));
if any(w < -1 | w > 2)
    NonGravForces.fail(unit, "inputRange", "height %.0f km is further outside the Swarm bracket [%.0f %.0f] km than the bracket is wide -> the scale is not defined there", ...
        h(find(w < -1 | w > 2, 1)) / 1e3, hLow(d(1)) / 1e3, hHigh(d(1)) / 1e3);
end
f = exp((1 - w) .* log(low(d)) + w .* log(high(d)));
end
