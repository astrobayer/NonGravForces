function [A, info] = attributeRevolutions(T, options)
%ATTRIBUTEREVOLUTIONS Attribute per-revolution scale variations to drivers: density, geometry or sensor.
%
%   A = NonGravForces.app.attributeRevolutions(T, 'Orbit', orb, 'SpaceWeather', sw)
%   [A, info] = NonGravForces.app.attributeRevolutions(T, 'Drivers', D, 'Families', struct('myDriver', "geometry"))
%
%   The attribution of CONCEPT 6.2 and 7 (N-03c; the v3 plan's WP-17
%   `driverRegression`): the per-revolution scales of app.perRevolution are
%   regressed on what could move them, and the regression says which family
%   of causes does. A scale that follows the solar flux is a density (drag
%   model) error showing through; one that follows beta' or the eclipse
%   fraction is a radiation-geometry error; one that follows time alone is
%   the sensor. The regression, its VIF gate and its block-bootstrap sigmas
%   are timeSeriesPreparation.regressDrivers' (CONCEPT 3.2: no estimator is
%   written twice); this unit supplies the drivers and reads the verdict.
%
%   THE DRIVERS, one value per revolution (row of T):
%     time              days since the first revolution's tMid       sensor
%     betaPrime         mean solar beta angle over the window, rad   geometry
%                       ('Orbit'; referenceSystem.orbitSunGeometry)
%     eclipseFraction   mean of 1 - nu over the window ('Orbit';      geometry
%                       env.eclipseFactor, conical)
%     f107, f107a, ap   at tMid ('SpaceWeather'; io.spaceWeather)     density
%     amp_<effect>      T's per-effect RMS on the parameter's axis    drag: density,
%                                                                     others: geometry
%     <caller's>        each variable of 'Drivers'                   'Families', else "other"
%   A driver that is constant or missing on every usable revolution carries no
%   information and is dropped before the regression; info.dropped names it.
%
%   THE VERDICT, per parameter: the families of the drivers whose bootstrap
%   |coef / sigma| exceeds 'Significance' (3, not 2: measured on 60 synthetic
%   revolutions, the block bootstrap's sigma of a driver with NO effect was
%   0.54 of the OLS sigma -- a block length of 14 leaves about four independent
%   blocks -- so a null driver read |z| = 2.7 and two verdicts out of two came
%   out "mixed" at 2. A few revolutions do not support a 2-sigma claim). One family names the verdict
%   ("density", "geometry", "sensor", or a caller's); several give "mixed";
%   none gives "unresolved" with reason "no significant driver". When
%   regressDrivers refuses for collinearity -- a VIF above 'MaxVif', the
%   CONCEPT's rule -- the verdict is "unresolved" with reason "collinear" and
%   its message: two drivers that move together cannot be told apart, and
%   saying so is the answer, not a failure.
%
% INPUT ARGUMENTS
%   Name  Size  Type   Units  Description
%   ----  ----  ----   -----  -----------
%   T     RxN   table  -      app.perRevolution's table (t0, t1, tMid, scale, optional amp_*)
%
% NAME-VALUE ARGUMENTS
%   Name          Size  Type                      Units  Description
%   ----          ----  ----                      -----  -----------
%   Orbit         1x1   struct                    m      .epoch .r .v (GCRS, as io.readGNV1B
%                                                        gives them); adds betaPrime and
%                                                        eclipseFraction; default none
%   SpaceWeather  1x1   NonGravForces.io.spaceWeather -  Adds f107, f107a, ap; default none
%   Drivers       RxN   table                     -      Further drivers, one numeric column each;
%                                                        default none
%   Families      1x1   struct                    -      Field = driver name, value = family; overrides
%                                                        and extends the defaults above
%   Significance  1x1   double                    -      |coef/sigma| above which a driver counts;
%                                                        default 3
%   MaxVif        1x1   double                    -      regressDrivers' VIF gate; default 10
%   NBoot         1x1   double                    -      regressDrivers' replicates; default 1000
%   Seed          1x1   double                    -      regressDrivers' seed; default 0
%   BlockLength   1x1   double                    -      regressDrivers' block length in revolutions;
%                                                        default [] = automatic (needs >= 20 rows)
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   A     Jx6   table   -      One row per estimated per-revolution scale: parameter ("S_xx", ...),
%                              verdict, reason, significant (string, the significant drivers),
%                              r2, nRows
%   info  1x1   struct  -      .drivers (R-row table of every driver offered) .families (struct)
%                              .dropped (string) .regression (1xJ cell of regressDrivers results,
%                              [] where it refused)
%
% ACCURACY
%   The regression's: OLS by QR, sigmas by moving-block bootstrap over
%   revolutions (regressDrivers). A verdict is a statement at 'Significance'
%   bootstrap sigmas about THIS arc; with R revolutions and K drivers it rests
%   on R - K - 1 degrees of freedom. THE AUTOMATIC BLOCK LENGTH NEEDS AT LEAST 20
%   REVOLUTIONS (about 1.3 days; timeSeriesPreparation's Politis-White rule):
%   one day (16) is refused with timeSeriesPreparation:blockLength:tooShort --
%   found by N-05, after this help had claimed a day supports a few drivers.
%   Below 20, pass 'BlockLength' explicitly and own that choice. betaPrime and
%   eclipseFraction are window means of exact geometry; the eclipse model is
%   conical.
%   DAILY INDICES NEED MANY DAYS. f107, f107a and ap change once a day, so over
%   a few days they are step functions with a few distinct values; a bootstrap
%   replicate built from blocks inside one day sees them constant, and its solve
%   inside regressDrivers is rank deficient (MATLAB:rankDeficientMatrix, measured
%   on 40 revolutions = 2.6 days). The point estimate is unaffected, the sigmas
%   are not to be trusted; attribute space weather over weeks (N-09 campaigns).
%
% ERRORS
%   NonGravForces:attributeRevolutions:inputTable    T lacks t0, t1, tMid or scale
%   NonGravForces:attributeRevolutions:inputOrbit    'Orbit' without .epoch, .r, .v of matching size
%   NonGravForces:attributeRevolutions:inputDrivers  'Drivers' without one numeric row per revolution
%   NonGravForces:attributeRevolutions:dataNoDriver  no driver varies over the usable revolutions
%   timeSeriesPreparation:regressDrivers:tooFewRows  fewer revolutions than drivers + 2
%   timeSeriesPreparation:blockLength:tooShort       fewer than 20 revolutions and no 'BlockLength'
%   referenceSystem:*, NonGravForces:spaceWeather:*  from the geometry and the index lookup
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md 6.2 and 7; the v3 plan
%   (docs/inputs/HANDOVER_AccelCal_v3.md) WP-17; Kutner M. et al. (2005)
%   Applied Linear Statistical Models, sect. 10.5 (VIF).
%
% EXAMPLE
%   rs = RandStream("twister", "Seed", 1); R = 40;
%   f = 70 + 10 * randn(rs, R, 1);
%   T = table(zeros(R, 1, 'int64'), zeros(R, 1, 'int64'), int64(1:R)' * int64(5640e9), ...
%       [0.9 + 0.005 * (f - 70) + 1e-3 * randn(rs, R, 1), nan(R, 2)], ...
%       'VariableNames', {'t0', 't1', 'tMid', 'scale'});
%   A = NonGravForces.app.attributeRevolutions(T, 'Drivers', table(f), 'Families', struct('f', "density"), 'NBoot', 200)
%
% See also NonGravForces.app.perRevolution, timeSeriesPreparation.regressDrivers, referenceSystem.orbitSunGeometry.
%
% KEYWORDS calibration attribution per-revolution drivers regression VIF density geometry sensor
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
arguments
    T table
    options.Orbit struct = struct([])
    options.SpaceWeather = []
    options.Drivers table = table()
    options.Families struct = struct()
    options.Significance (1,1) double {mustBePositive} = 3
    options.MaxVif (1,1) double {mustBeGreaterThanOrEqual(options.MaxVif, 1)} = 10
    options.NBoot (1,1) double {mustBeInteger, mustBePositive} = 1000
    options.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.BlockLength double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
end
unit = "attributeRevolutions";
if ~all(ismember(["t0" "t1" "tMid" "scale"], string(T.Properties.VariableNames))) || size(T.scale, 2) ~= 3
    NonGravForces.fail(unit, "inputTable", "T must be app.perRevolution's table, with t0, t1, tMid and a three-column scale -> got %s", strjoin(string(T.Properties.VariableNames), ", "));
end
R = height(T);

% --- the drivers -------------------------------------------------------------------------
fam = struct('time', "sensor", 'betaPrime', "geometry", 'eclipseFraction', "geometry", ...
    'f107', "density", 'f107a', "density", 'ap', "density", 'amp_drag', "density", ...
    'amp_radiation', "geometry", 'amp_thermal', "geometry");
D = table(double(T.tMid - T.tMid(1)) / 86400e9, 'VariableNames', {'time'});
if ~isempty(options.Orbit)
    D = [D, localOrbitDrivers(T, options.Orbit, unit)];
end
if ~isempty(options.SpaceWeather)
    d = options.SpaceWeather.query(referenceSystem.gpsNanosecondsToEpoch(T.tMid));
    D.f107 = d.f107(:); D.f107a = d.f107a(:); D.ap = d.ap(:);
end
if width(options.Drivers) > 0
    if height(options.Drivers) ~= R || ~all(varfun(@(x) isnumeric(x) && iscolumn(x), options.Drivers, 'OutputFormat', 'uniform'))
        NonGravForces.fail(unit, "inputDrivers", "'Drivers' must hold one numeric column per driver and %d rows, one per revolution -> got %d rows", R, height(options.Drivers));
    end
    D = [D, options.Drivers];
end
for f = string(fieldnames(options.Families))'
    fam.(f) = string(options.Families.(f));
end
ampNames = string(T.Properties.VariableNames(startsWith(T.Properties.VariableNames, "amp_")));

% --- one regression per estimated per-revolution scale -----------------------------------
ax = ["x" "y" "z"];
axList = find(any(isfinite(T.scale), 1));
J = numel(axList);
param = strings(J, 1); verdict = strings(J, 1); reason = strings(J, 1); sig = strings(J, 1);
r2 = nan(J, 1); nRows = nan(J, 1); regs = cell(1, J); dropped = strings(1, 0);
for jj = 1:J
    i = axList(jj);
    param(jj) = "S_" + ax(i) + ax(i);
    Zt = D;
    for a = ampNames, Zt.(a) = T.(a)(:, i); end
    P = T.scale(:, i);
    rowsP = isfinite(P);
    names = string(Zt.Properties.VariableNames);
    Z = Zt{:, :};
    keep = false(1, numel(names));
    for k = 1:numel(names)
        z = Z(rowsP & isfinite(Z(:, k)), k);
        keep(k) = numel(z) >= 2 && any(z ~= z(1));
    end
    dropped = unique([dropped, names(~keep)]);
    if ~any(keep)
        NonGravForces.fail(unit, "dataNoDriver", "no driver varies over the %d usable revolutions of %s -> supply 'Orbit', 'SpaceWeather' or 'Drivers'", nnz(rowsP), param(jj));
    end
    names = names(keep); Z = Z(:, keep);
    try
        res = timeSeriesPreparation.regressDrivers(P, Z, Names = param(jj), DriverNames = names, MaxVif = options.MaxVif, ...
            NBoot = options.NBoot, Seed = options.Seed, Stage = "attributeRevolutions." + param(jj), BlockLength = options.BlockLength);
    catch err
        if err.identifier ~= "timeSeriesPreparation:regressDrivers:collinear", rethrow(err); end
        verdict(jj) = "unresolved"; reason(jj) = "collinear: " + string(err.message);
        continue
    end
    regs{jj} = res;
    r2(jj) = res.r2; nRows(jj) = res.nRows;
    z = abs(res.coef(2:end) ./ res.sigma(2:end));
    hit = names(z(:)' > options.Significance);
    sig(jj) = strjoin(hit, ", ");
    fams = strings(1, numel(hit));
    for k = 1:numel(hit)
        if isfield(fam, hit(k)), fams(k) = fam.(hit(k)); else, fams(k) = "other"; end
    end
    fams = unique(fams);
    if isempty(fams), verdict(jj) = "unresolved"; reason(jj) = "no significant driver";
    elseif isscalar(fams), verdict(jj) = fams; reason(jj) = "";
    else, verdict(jj) = "mixed"; reason(jj) = strjoin(fams, " + ");
    end
end
A = table(param, verdict, reason, sig, r2, nRows, ...
    'VariableNames', {'parameter', 'verdict', 'reason', 'significant', 'r2', 'nRows'});
info = struct('drivers', D, 'families', fam, 'dropped', dropped, 'regression', {regs});
end

% =========================================================================================
function Z = localOrbitDrivers(T, orb, unit)
%LOCALORBITDRIVERS Window means of beta' and of the eclipse fraction along the orbit.
if ~all(isfield(orb, {'epoch', 'r', 'v'})) || size(orb.r, 1) ~= 3 || ~isequal(size(orb.r), size(orb.v)) ...
        || numel(orb.epoch.day) ~= size(orb.r, 2)
    NonGravForces.fail(unit, "inputOrbit", "'Orbit' needs .epoch, .r and .v (3xN GCRS) of one size -> e.g. io.readGNV1B");
end
ep = orb.epoch;
if string(ep.scale) ~= "GPS", ep = referenceSystem.convertTimeScale(ep, "GPS"); end
t = referenceSystem.epochTo(ep, "gpsNanoseconds");
beta = referenceSystem.orbitSunGeometry(orb.r, orb.v, ep);
sun = referenceSystem.geocentricBodyState("sun", ep);
ecl = 1 - NonGravForces.env.eclipseFactor(orb.r, sun);
R = height(T);
betaPrime = nan(R, 1); eclipseFraction = nan(R, 1);
for w = 1:R
    in = t >= T.t0(w) & t < T.t1(w);
    if any(in)
        betaPrime(w) = mean(beta(in));
        eclipseFraction(w) = mean(ecl(in));
    end
end
Z = table(betaPrime, eclipseFraction);
end
