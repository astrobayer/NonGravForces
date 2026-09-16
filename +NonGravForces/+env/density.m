function [rho, species, T, info] = density(model, epoch, lla, sw, options)
%DENSITY Thermosphere density and composition: exponential, NRLMSISE-00 or NRLMSIS 2.1 (both compiled).
%
%   [rho, species, T, info] = NonGravForces.env.density("msis21", epoch, lla, sw)
%   [rho, species, T, info] = NonGravForces.env.density("msis00", epoch, lla, sw)
%   [rho, species, T, info] = NonGravForces.env.density("msis00", epoch, lla, [], 'Drivers', d)
%   rho = NonGravForces.env.density("exponential", epoch, lla, [], 'Exponential', p)
%
%   K-vectorised evaluation of the neutral atmosphere at geodetic positions.
%   The position is GEODETIC ON WGS-84 (latitude, longitude, height), which
%   the caller obtains from referenceSystem.cartesianToGeodetic with
%   'Ellipsoid', "WGS84" passed explicitly (its default is GRS80): a spherical
%   Earth is 21.4 km wrong at the poles and a factor 1.40 in density (HANDOVER
%   C.1, defect D2).
%
%   BACKENDS
%     msis00       NRLMSISE-00 (Picone et al. 2002) through the MEX of the
%                  Brodowski C fetched by setup (NonGravForces.buildMex). Drivers from
%                  a celestrak-family spaceWeather struct queried at the
%                  epochs: previous-day observed F10.7, 81-day centred mean,
%                  the 7-element ap history (switch 9 = -1); or from
%                  'Drivers'. gtd7d by default: the total mass density
%                  includes anomalous oxygen, as the model authors recommend
%                  for drag above ~500 km; 'AnomalousOxygen', false gives gtd7.
%     msis21       NRLMSIS 2.1 (NRL 2022; formulation Emmert et al. 2021) through the MEX of NRL's
%                  unmodified Fortran, compiled in double precision
%                  (NonGravForces.buildMex; setupNonGravForces fetches the package).
%                  The same drivers as msis00, with the same meaning (MSISCALC's
%                  header: previous-day F10.7, 81-day centred mean, ap history read
%                  when 7xK). The model forms local time from UT and longitude
%                  itself, so 'LocalSolarTime_h' is refused. The day of year is the
%                  integer day, as the legacy interface NRL tests the model with.
%                  The mass density includes anomalous oxygen unless
%                  'AnomalousOxygen' is false, which is passed to the model as its
%                  own mass-density flag, not subtracted here. Where NRLMSIS marks a
%                  species as not modelled (its value 9.999e-38, e.g. O, H and N in
%                  the fully mixed lower atmosphere of NRL's test cases at 0.2 km),
%                  that species is returned as 0.
%                  NRLMSIS 2.1 is NRLMSIS 2.0 plus nitric oxide: its temperature and
%                  every species this toolbox reads are those of 2.0.
%     exponential  rho = rho0 exp(-(h - h0)/H) with the caller's rho0, h0, H
%                  ('Exponential'); sizing only, no composition, no drivers.
%   DTM2020 and JB2008 are not provided (owner ruling 2026-09-15, HANDOVER Part I).
%
%   NaN POLICY: a NaN anywhere in lla or in the drivers raises numericNaN,
%   because a NaN density silently zeroes the drag (HANDOVER E.2 rule 7).
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units      Description
%   ----     ----  ----    -----      -----------
%   model    1x1   string  -          "exponential" | "msis00" | "msis21"
%   epoch    1x1   struct  -          referenceSystem epoch struct, K entries, any scale
%   lla      3xK   double  rad,rad,m  Geodetic latitude, longitude, height (WGS-84)
%   sw       1x1   struct  -          NonGravForces.io.spaceWeather result of the
%                                     celestrak family, or [] with 'Drivers' or
%                                     for "exponential".
%
% NAME-VALUE ARGUMENTS
%   Name              Size  Type     Units    Description
%   ----              ----  ----     -----    -----------
%   Drivers           1x1   struct   -        .f107 (1xK, previous-day observed
%                                             F10.7), .f107a (1xK, 81-day
%                                             centred), .ap (1xK daily Ap, or 7xK
%                                             NRLMSISE history). Overrides sw.
%   AnomalousOxygen   1x1   logical  -        true (default): the mass density includes
%                                             anomalous oxygen (msis00: gtd7d; msis21: the
%                                             model's default); false: it does not (gtd7).
%   Exponential       1x1   struct   -        .rho0_kg_m3, .h0_m, .H_m; required
%                                             for "exponential", ignored otherwise.
%   LocalSolarTime_h  1xK   double   h        msis00 only. Override of the local solar time handed
%                                             to MSIS. Default [] = sec/3600 + lon/15,
%                                             the consistent choice the reference
%                                             recommends. The published test cases
%                                             set an inconsistent 16 h and are
%                                             reproduced only with this override.
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units    Description
%   ----     ----  ----    -----    -----------
%   rho      1xK   double  kg/m^3   Total mass density.
%   species  8xK   double  1/m^3    Number densities He O N2 O2 Ar H N
%                                   O_anomalous (const.speciesName order); NaN
%                                   for "exponential".
%   T        2xK   double  K        [exospheric temperature; temperature at h];
%                                   NaN for "exponential".
%   info     1x1   struct  -        .backend, .drivers (struct actually used:
%                                   f107 f107a ap, and .source "spaceWeather" |
%                                   "caller" | "none"), .anomalousOxygen,
%                                   .msisInput (struct doy sec alt_km lat_deg
%                                   lon_deg, 1xK each, and lst_h for msis00;
%                                   both msis backends), .provenance.
%
% ACCURACY
%   msis00 is bit-identical to the reference C for the same inputs (same
%   object code); the SI conversion (x1e6 for number densities, x1e3 for the
%   mass density) is one multiplication. The published 17-case table is
%   reproduced to its printed 7 digits every run (tests/tNgfEnv.m). Local
%   solar time is sec/3600 + lon/15 as the reference documents; the equation
%   of time (< 17 min) is neglected, as the reference recommends. The model
%   itself is uncertain by 10-30 % (HANDOVER C.4).
%   msis21 returns the model's SI outputs unconverted. NRL's 200 reference cases
%   (msis2.1_test_ref_dp.txt, printed to 4 significant digits) are reproduced
%   within half a unit of the last printed digit every run. The model error of
%   the climatology itself is not measured here; the gas-surface interaction and
%   the drivers' daily and 3-hourly resolution bound what either MSIS can give.
%
% ERRORS
%   NonGravForces:density:inputShape        lla not 3xK or K differs from epoch
%   NonGravForces:density:inputDrivers      no sw and no Drivers; sw not of the
%                                           celestrak family; Drivers or
%                                           LocalSolarTime_h malformed
%   NonGravForces:density:inputExponential  "exponential" without a complete
%                                           Exponential struct
%   NonGravForces:density:numericNaN        NaN in lla or in the drivers
%   NonGravForces:density:depMissing        the backend's MEX is not built
%   NonGravForces:density:inputOption       LocalSolarTime_h given to msis21, which forms local time itself
%   NonGravForces:msis21Gateway:*           path, parameter file and screening errors of the gateway, unwrapped
%   NonGravForces:spaceWeather:*            from sw.query, unwrapped
%   referenceSystem:*                       propagated unwrapped
%
% REFERENCES
%   Picone, Hedin, Drob, Aikin (2002) JGR 107(A12) 1468, NRLMSISE-00.
%   Brodowski, NRLMSISE-00 C release 20041227, DOCUMENTATION §5 (inputs) and §7 (table).
%   Emmert J.T. et al. (2021) NRLMSIS 2.0, Earth Space Sci. 8, e2020EA001321;
%   NRL, nrlmsis2.1 package (2022-04-04): readme.txt and the header of msis_calc.F90.
%
% EXAMPLE
%   ep = referenceSystem.mjdToEpoch(58849.5, "UTC");
%   lla = [deg2rad(60); deg2rad(-70); 400e3];
%   d = struct('f107', 150, 'f107a', 150, 'ap', 4);
%   rho = NonGravForces.env.density("msis00", ep, lla, [], 'Drivers', d, 'AnomalousOxygen', false)
%
% See also NonGravForces.io.spaceWeather, NonGravForces.buildMex, referenceSystem.cartesianToGeodetic.
%
% KEYWORDS density thermosphere NRLMSISE-00 NRLMSIS exponential drag
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic); msis21 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    model (1,1) string {mustBeMember(model, ["exponential","msis00","msis21"])}
    epoch (1,1) struct
    lla double
    sw = []
    options.Drivers = []
    options.AnomalousOxygen (1,1) logical = true
    options.Exponential = []
    options.LocalSolarTime_h = []
end
unit = "density";
if ~all(isfield(epoch, ["day","sec","scale"]))
    NonGravForces.fail(unit, "inputShape", "epoch must be a referenceSystem epoch struct (.day .sec .scale) -> build one with daySecondToEpoch or mjdToEpoch");
end
K = numel(epoch.day);
if ~(ismatrix(lla) && size(lla, 1) == 3 && size(lla, 2) == K && isreal(lla))
    NonGravForces.fail(unit, "inputShape", "lla must be 3x%d [lat; lon; h] for %d epochs, got %s -> transpose or match the epoch count", K, K, mat2str(size(lla)));
end
if any(~isfinite(lla(:)))
    NonGravForces.fail(unit, "numericNaN", "lla contains NaN/Inf at column %d -> screen the positions first; a NaN density would silently zero the drag", find(any(~isfinite(lla), 1), 1));
end
info.backend = model;
switch model
    case "exponential"
        p = options.Exponential;
        if ~(isstruct(p) && all(isfield(p, ["rho0_kg_m3","h0_m","H_m"]))) || ~all(isfinite([p.rho0_kg_m3 p.h0_m p.H_m])) || p.rho0_kg_m3 <= 0 || p.H_m <= 0
            NonGravForces.fail(unit, "inputExponential", "model exponential needs 'Exponential', struct('rho0_kg_m3', r, 'h0_m', h0, 'H_m', H) with positive rho0 and H -> pass the scale height model explicitly (sizing only)");
        end
        rho = p.rho0_kg_m3 * exp(-(lla(3, :) - p.h0_m) / p.H_m);
        species = nan(8, K); T = nan(2, K);
        info.drivers = struct('source', "none");
        info.anomalousOxygen = false;
        info.exponential = p;
    case "msis00"
        drv = localDrivers(unit, epoch, sw, options.Drivers, K);
        gw = fullfile(fileparts(mfilename('fullpath')), 'private', "msisGateway." + mexext);
        if ~isfile(gw)
            NonGravForces.fail(unit, "depMissing", "the NRLMSISE-00 MEX %s is not built -> run NonGravForces.buildMex once on this machine", gw);
        end
        utc = referenceSystem.convertTimeScale(epoch, "UTC");
        dt = referenceSystem.epochTo(utc, "datetime");
        in.doy     = double(day(dt, 'dayofyear'))';
        in.sec     = double(utc.sec(:))';
        in.alt_km  = lla(3, :) / 1000;
        in.lat_deg = rad2deg(lla(1, :));
        in.lon_deg = rad2deg(lla(2, :));
        if isempty(options.LocalSolarTime_h)
            in.lst_h = in.sec / 3600 + in.lon_deg / 15;      % the consistent choice the reference recommends
        else
            lst = double(options.LocalSolarTime_h(:)');
            if isscalar(lst), lst = repmat(lst, 1, K); end
            if numel(lst) ~= K || any(~isfinite(lst))
                NonGravForces.fail(unit, "inputDrivers", "LocalSolarTime_h must be 1x%d finite hours -> match the epoch count", K);
            end
            in.lst_h = lst;
        end
        [d, t] = msisGateway(in.doy, in.sec, in.alt_km, in.lat_deg, in.lon_deg, in.lst_h, drv.f107a, drv.f107, drv.ap, double(options.AnomalousOxygen));
        rho = d(6, :) * 1e3;                                   % g/cm^3 -> kg/m^3
        species = d([1 2 3 4 5 7 8 9], :) * 1e6;               % cm^-3 -> m^-3, He O N2 O2 Ar H N O_anom
        T = t;
        info.drivers = drv;
        info.anomalousOxygen = options.AnomalousOxygen;
        info.msisInput = in;
    case "msis21"
        if ~isempty(options.LocalSolarTime_h)
            NonGravForces.fail(unit, "inputOption", "msis21 forms the local solar time from UT and longitude itself -> remove 'LocalSolarTime_h' (it serves msis00's published test cases only)");
        end
        drv = localDrivers(unit, epoch, sw, options.Drivers, K);
        here = fullfile(fileparts(mfilename('fullpath')), 'private');
        gw = fullfile(here, "msis21Gateway." + mexext);
        if ~isfile(gw)
            NonGravForces.fail(unit, "depMissing", "the NRLMSIS 2.1 MEX %s is not built -> run setupNonGravForces (fetches NRL's package) and NonGravForces.buildMex once on this machine", gw);
        end
        utc = referenceSystem.convertTimeScale(epoch, "UTC");
        dt = referenceSystem.epochTo(utc, "datetime");
        in = struct('doy', double(day(dt, 'dayofyear'))', 'sec', double(utc.sec(:))', 'alt_km', lla(3, :) / 1000, ...
                    'lat_deg', rad2deg(lla(1, :)), 'lon_deg', rad2deg(lla(2, :)));
        [d, t] = msis21Gateway([here filesep], double(options.AnomalousOxygen), in.doy, in.sec, in.alt_km, in.lat_deg, in.lon_deg, drv.f107a, drv.f107, drv.ap);
        dmissing = 9.999e-38;                                  % msis_constants.F90: the model's marker for "not modelled here"
        d(d == dmissing) = 0;
        rho = d(1, :);                                         % kg/m^3, SI from the model
        species = d([5 4 2 3 7 6 8 9], :);                     % MSIS 2.1 order rho N2 O2 O He H Ar N Oanom NO -> He O N2 O2 Ar H N Oanom
        T = t;
        info.drivers = drv;
        info.anomalousOxygen = options.AnomalousOxygen;
        info.msisInput = in;
end
info.provenance = NonGravForces.provenance("env.density", 'Options', struct('model', model, 'anomalousOxygen', options.AnomalousOxygen, 'lstOverride', ~isempty(options.LocalSolarTime_h)));
end

function drv = localDrivers(unit, epoch, sw, drivers, K)
if ~isempty(drivers)
    if ~(isstruct(drivers) && all(isfield(drivers, ["f107","f107a","ap"])))
        NonGravForces.fail(unit, "inputDrivers", "Drivers must have fields f107, f107a, ap -> see the help");
    end
    f107 = double(drivers.f107(:)'); f107a = double(drivers.f107a(:)'); ap = double(drivers.ap);
    if isscalar(f107), f107 = repmat(f107, 1, K); end
    if isscalar(f107a), f107a = repmat(f107a, 1, K); end
    if isscalar(ap), ap = repmat(ap, 1, K); end
    if isvector(ap) && numel(ap) == 7 && K == 1, ap = ap(:); end
    if isrow(ap) && numel(ap) == K && size(ap, 1) == 1, ap = reshape(ap, 1, K); end
    if ~(numel(f107) == K && numel(f107a) == K && size(ap, 2) == K && any(size(ap, 1) == [1 7]))
        NonGravForces.fail(unit, "inputDrivers", "Drivers must be 1x%d (f107, f107a) and 1x%d or 7x%d (ap) -> match the epoch count", K, K, K);
    end
    src = "caller";
elseif isstruct(sw) && isfield(sw, 'format') && isfield(sw, 'query')
    if sw.format ~= "celestrak"
        NonGravForces.fail(unit, "inputDrivers", "spaceWeather family %s does not drive NRLMSIS -> load the CelesTrak SW-All file (family celestrak)", sw.format);
    end
    q = sw.query(epoch);
    f107 = q.f107; f107a = q.f107a; ap = q.apHistory;
    src = "spaceWeather";
else
    NonGravForces.fail(unit, "inputDrivers", "the msis backends need drivers -> pass a celestrak spaceWeather struct as sw, or 'Drivers', struct('f107', .., 'f107a', .., 'ap', ..)");
end
if any(~isfinite([f107(:); f107a(:); ap(:)]))
    NonGravForces.fail(unit, "numericNaN", "a driver is NaN/Inf -> spaceWeather never returns NaN; check the caller's Drivers");
end
drv = struct('f107', f107, 'f107a', f107a, 'ap', ap, 'source', src);
end
