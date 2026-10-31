function S = densityDriverScale(sw, cls, options)
%DENSITYDRIVERSCALE Daily scale on NRLMSISE-00 from solar and geomagnetic drivers alone, for any date.
%
%   S = NonGravForces.env.densityDriverScale(sw, cls)
%   S = NonGravForces.env.densityDriverScale(sw, cls, 'Offset', -0.041)   % GRACE-FO's own level
%   f = S.at(epoch, h_m)                        % 1xK factor on the model density
%
%   The extrapolable counterpart of env.densityScale (HANDOVER O-19, step 4). env.densityScale
%   measures the factor f = rho_observed / rho_msis00 per day from Swarm's densities and so
%   exists only where Swarm flew; this unit PREDICTS it from drivers that exist for every date
%   since 1951, with the regression tools/densityDriverFit.m selected by yearly blocked
%   cross-validation over 2000-2026 and six satellites (Swarm A/B/C POD-only, TU Delft's
%   CHAMP, GRACE-A and GRACE-FO 1 accelerometer densities):
%
%       log f = a0 + cF Fbar + cdF dF + cap ap + cs1 sin(p) + cc1 cos(p) + cs2 sin(2p)
%               + cc2 cos(2p) + cQ ((Fbar - Fc)/100)^2 + offset
%
%   Fbar   the 81-day mean of CLS's adjusted 30 cm flux F30 over the 81 UTC days ENDING THE DAY
%          BEFORE (no future data), sfu
%   dF     F30 of the day before minus Fbar, sfu
%   ap     the daily Ap of the day (CelesTrak SW-All)
%   p      2 pi (day of year - 1) / 365.25 of the UTC day
%   Fc     the centre of the square, the data set's mean Fbar
%
%   One factor per UTC day; it does not depend on height or local time (adding them did not
%   lower the cross-validated error). log f peaks near Fbar ~ 99 sfu of F30 and falls on
%   both sides: NRLMSISE-00 overestimates most at solar minimum (f ~ 0.5) and least at medium
%   activity (f ~ 0.9). The coefficients are the fit of run 40676249479 ("F30 + season +
%   square", trimmed days); 'Offset' adds a satellite's own level in log f, -0.041 for
%   GRACE-FO 1 (its accelerometer densities read 4 % below Swarm's in the same years).
%
%   A Model applies S as its 'DensityScale' under Options drag.densityScale "driver"; under
%   "combined" (the standard preset since 0.17.0) S is the 'Fallback' of env.densityScale and
%   serves the days without a Swarm factor.
%
% INPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   sw    1x1   struct  -      io.spaceWeather result, celestrak family (SW-All): the daily Ap
%   cls   1x1   struct  -      io.spaceWeather result, f30 family, of CLS's ADJUSTED file
%                              (radio_flux_adjusted_observation.txt): the daily F30
%
% NAME-VALUE ARGUMENTS
%   Name    Size  Type    Units  Description
%   ----    ----  ----    -----  -----------
%   Model   1x1   string  -      density backend the factor is for; only "msis00" (default), the
%                                model the coefficients were fitted against
%   Offset  1x1   double  log f  added to log f: a satellite's own level (default 0, Swarm's;
%                                -0.041 for GRACE-FO 1, +0.036 CHAMP, +0.045 GRACE-A, run 40676249479)
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   S     1x1   struct  -      .kind "driver", .model "msis00", .proxy "F30", .coef (struct of the
%                              transcribed coefficients, .a0 .Fbar .dF .ap .sinAnnual .cosAnnual
%                              .sinSemiannual .cosSemiannual .FbarSq .centre_sfu .run), .offset,
%                              .drivers (function handle: epoch -> 8xK drivers in the order above),
%                              .at (function handle: epoch struct (any scale), 1xK heights [m] ->
%                              1xK factors; the heights are checked for size and otherwise unused),
%                              .provenance
%
% ACCURACY
%   The regression's held-out error is 0.17-0.18 in log f per day (blocked CV by year 0.172,
%   forward across 2020 0.181, backward 0.180), against 0.23-0.25 for a constant factor: one
%   day's predicted factor is good to about 18 %. Mission levels differ by 3-5 % ('Offset').
%   The transcription is exact: log f equals the fit's own prediction at its reference day to
%   rounding (tests). Outside the drivers the fit saw (2000-08..2026-08, F30 Fbar from ~45 sfu in the
%   2008/2019 minima to the cycle-23 maximum) the square extrapolates and is not validated.
%
% ERRORS
%   NonGravForces:densityDriverScale:inputModel   Model is not msis00
%   NonGravForces:densityDriverScale:inputFamily  sw is not the celestrak family, or cls not f30, or cls is the
%                                                 flux not adjusted to 1 AU (.adjusted false)
%   NonGravForces:densityDriverScale:dataDay      S.at at a UTC day whose 81 previous days of F30 or
%                                                 whose Ap are not all in the files
%   NonGravForces:densityDriverScale:inputRange   h not 1xK, one height per epoch
%
% REFERENCES
%   HANDOVER O-19 (steps 1-3; the fit, its validation and the coefficient record of
%   run 40676249479); Dudok de Wit, Bruinsma & Shibasaki (2014), J. Space Weather Space Clim.
%   4, A06 (F30 as the better thermospheric proxy); Picone et al. (2002), J. Geophys. Res. 107(A12),
%   1468 (NRLMSISE-00); March et al. (2019), Adv. Space Res. 63, 213-238 (TU Delft densities).
%
% EXAMPLE
%   root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
%   sw  = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
%   cls = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
%   S = NonGravForces.env.densityDriverScale(sw, cls);
%   S.at(referenceSystem.mjdToEpoch(58849.5, "UTC"), 490e3)
%
% See also NonGravForces.env.densityScale, NonGravForces.env.density, NonGravForces.Model, NonGravForces.io.spaceWeather.
%
% KEYWORDS density scale driver regression F30 extrapolation correction drag thermosphere
% Created 2026-10-28 by Claude (Anthropic) for Matthias Weigelt.
arguments
    sw (1,1) struct
    cls (1,1) struct
    options.Model (1,1) string = "msis00"
    options.Offset (1,1) double {mustBeFinite} = 0
end
unit = "densityDriverScale";
if options.Model ~= "msis00"
    NonGravForces.fail(unit, "inputModel", "the coefficients were fitted against msis00, not %s -> use drag.density msis00, or refit with tools/densityDriverFit.m 'Model'", options.Model);
end
if ~isfield(sw, 'index') || ~isfield(sw.index, 'apDaily') || ~isfield(cls, 'index') || ~isfield(cls.index, 'f30')
    NonGravForces.fail(unit, "inputFamily", "sw must be io.spaceWeather of SW-All (celestrak) and cls of CLS's adjusted radio flux (f30) -> read the two files with io.spaceWeather");
end
if ~isfield(cls, 'adjusted') || ~isequal(cls.adjusted, true)
    % the coefficients were fitted to the flux ADJUSTED to 1 AU; the absolute flux differs by +-3.4 % over the year
    NonGravForces.fail(unit, "inputFamily", "cls is CLS's flux NOT adjusted to 1 AU, the coefficients were fitted to the adjusted one -> read radio_flux_adjusted_observation.txt");
end
% THE COEFFICIENT RECORD, transcribed from tools/densityDriverFit.m's full fit (trimmed days, F30), every digit.
c = struct('run', "40676249479", 'centre_sfu', 75.712339474163883, 'a0', -0.70150825262940808, ...
    'Fbar', 0.0065937301051911179, 'dF', 0.0016375596570947545, 'ap', 0.0012805481671233857, ...
    'sinAnnual', 0.03526596776917159, 'cosAnnual', 0.01288830194051793, ...
    'sinSemiannual', 0.029905119164486137, 'cosSemiannual', -0.022318436087852906, ...
    'FbarSq', -1.4039269195633037);
beta = [c.Fbar; c.dF; c.ap; c.sinAnnual; c.cosAnnual; c.sinSemiannual; c.cosSemiannual; c.FbarSq];
mjdF = double(cls.mjd(:)); f30 = double(cls.index.f30(:));
mjdA = double(sw.mjd(:)); ap = double(sw.index.apDaily(:));
S.kind = "driver";
S.model = options.Model;
S.proxy = "F30";
S.coef = c;
S.offset = options.Offset;
S.drivers = @(epoch) localDrivers(unit, mjdF, f30, mjdA, ap, c.centre_sfu, epoch);
S.at = @(epoch, h) localAt(unit, mjdF, f30, mjdA, ap, c, beta, options.Offset, epoch, h);
S.provenance = NonGravForces.provenance("env.densityDriverScale", 'Options', options);
end

function f = localAt(unit, mjdF, f30, mjdA, ap, c, beta, offset, epoch, h)
%LOCALAT The factor of each epoch's UTC day: exp(a0 + sum(beta .* drivers) + offset).
z = localDrivers(unit, mjdF, f30, mjdA, ap, c.centre_sfu, epoch);
if ~isequal(size(h), [1 size(z, 2)])
    NonGravForces.fail(unit, "inputRange", "h must be 1x%d, one height per epoch -> pass lla(3, :)", size(z, 2));
end
f = exp(c.a0 + sum(beta .* z, 1) + offset);        % column-wise: a BLAS product sums in a batch-size-dependent order
end

function z = localDrivers(unit, mjdF, f30, mjdA, ap, centre, epoch)
%LOCALDRIVERS 8xK drivers of the epochs' UTC days, in the order of the coefficient record.
utc = referenceSystem.convertTimeScale(epoch, "UTC");
utcDay = double(reshape(utc.day, 1, []));
[u, ~, back] = unique(utcDay);
Z = nan(8, numel(u));
for k = 1:numel(u)
    [ok, i] = ismember(u(k) - (81:-1:1)', mjdF);     % the 81 days ending the day before
    [okA, j] = ismember(u(k), mjdA);
    if ~all(ok) || ~okA || any(~isfinite(f30(i))) || ~isfinite(ap(j))
        NonGravForces.fail(unit, "dataDay", "UTC day MJD %d needs F30 for MJD %d..%d and the day's Ap -> pass space-weather files that cover them", ...
            u(k), u(k) - 81, u(k) - 1);
    end
    Fbar = mean(f30(i));
    t = datetime(u(k) + 2400000.5, 'ConvertFrom', 'juliandate');
    p = 2 * pi * (day(t, 'dayofyear') - 1) / 365.25;
    Z(:, k) = [Fbar; f30(i(end)) - Fbar; ap(j); sin(p); cos(p); sin(2 * p); cos(2 * p); ((Fbar - centre) / 100) ^ 2];
end
z = Z(:, back);
end
