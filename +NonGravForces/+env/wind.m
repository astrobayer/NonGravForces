function [vrel, info] = wind(model, epoch, r_gcrs, v_gcrs, options)
%WIND Velocity of the satellite relative to the air, GCRS: co-rotation, optional thermospheric wind.
%
%   [vrel, info] = NonGravForces.env.wind("corotation", epoch, r, v, 'Eop', eop)
%   [vrel, info] = NonGravForces.env.wind("none", epoch, r, v)
%   [vrel, info] = NonGravForces.env.wind("corotation", epoch, r, v, 'Eop', eop, 'Wind_m_s', w)
%   [vrel, info] = NonGravForces.env.wind("hwm14", epoch, r, v, 'Eop', eop, 'SpaceWeather', sw)
%   [vrel, info] = NonGravForces.env.wind("hwm14", epoch, r, v, 'Eop', eop, 'Disturbed', false)
%
%   The drag kernel needs the velocity of the satellite THROUGH THE GAS, not
%   through the inertial frame. At 7.6 km/s the co-rotation of the atmosphere is
%   up to 470 m/s at the equator, so omitting it is a 6 % error in the
%   cross-track drag and a 0.4 % error in its magnitude (Doornbos 2011).
%
%   MODELS
%     "none"        vrel = v. The caller states that the gas is inertial, which
%                   is never true at LEO; provided for isolating the effect.
%     "corotation"  The air co-rotates rigidly with the Earth. This unit does NOT
%                   form omega x r: it takes the Earth's rotation from the SAME
%                   authority the frames come from,
%                       vrel = v - Rdot * (R' * r),
%                   with R and Rdot from referenceSystem.rotationItrsToGcrs.
%                   Rdot carries the polar motion, the length-of-day variation and
%                   the precession-nutation rate, none of which a nominal
%                   omega x r contains (rule R1, one authority per fact; the
%                   predecessor's hand-written term was defect D4).
%     "hwm14"       The co-rotation PLUS the Horizontal Wind Model 14 (Drob et al.
%                   2015): quiet-time winds (HWM14) and, with 'Disturbed' true, the
%                   storm-induced disturbance winds (DWM07, Emmert et al. 2008) driven
%                   by the current 3-hour ap. NRL's unmodified Fortran through the MEX
%                   NonGravForces.buildMex builds; setupNonGravForces fetches it. The
%                   model is evaluated at the geodetic position (WGS-84, from the ITRS
%                   position through R) and UT; its meridional (+north) and zonal
%                   (+east) components are rotated from the local geodetic frame to
%                   the ITRS and then to the GCRS with R, and subtracted as well:
%                       vrel = v - Rdot (R' r) - R (w_north n + w_east e).
%                   Quiet-time winds are 40-80 m/s RMS, 1 % of the relative speed and
%                   the leading term in cross-track drag.
%
%   A caller who has a wind field from elsewhere passes it through 'Wind_m_s' and
%   it is subtracted after the co-rotation, whatever the model.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   model    1x1   string  -      "none" | "corotation" | "hwm14"
%   epoch    1x1   struct  -      referenceSystem epoch struct with K entries
%   r_gcrs   3xK   double  m      Satellite position, GCRS
%   v_gcrs   3xK   double  m/s    Satellite velocity, GCRS
%
% NAME-VALUE ARGUMENTS
%   Name       Size  Type    Units  Description
%   ----       ----  ----    -----  -----------
%   Eop           1x1      struct   -    referenceSystem.readEop result, or a file path;
%                                        required by "corotation" and "hwm14" (no default:
%                                        the rotation cannot be guessed, rule R5)
%   Wind_m_s      3xK      double   m/s  Air velocity in GCRS to subtract in addition,
%                                        default [] (none)
%   SpaceWeather  1x1      struct   -    io.spaceWeather (celestrak family): the current
%                                        3-hour ap at each epoch. "hwm14" only.
%   Ap            1x1|1xK  double   2 nT The current 3-hour ap, given instead of
%                                        SpaceWeather. "hwm14" only.
%   Disturbed     1x1      logical  -    "hwm14": add DWM07 (default true, which needs an
%                                        ap source); false: quiet-time winds only
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   vrel     3xK   double  m/s    Satellite velocity relative to the air, GCRS
%   info     1x1   struct  -      .model, .corotation_m_s (3xK, the term removed),
%                                 .hwm14 (hwm14: .northEast_m_s 2xK, .gcrs_m_s 3xK,
%                                 .ap 1xK or [], .disturbed, .geodetic 3xK [lat; lon; h]
%                                 rad rad m), .speed_m_s (1xK, |vrel|), .windGiven
%                                 (logical), .K, .provenance
%
% ACCURACY
%   Exact given R and Rdot: one matrix product and one subtraction, rounding
%   1e-16 relative. The co-rotation magnitude was measured at 255.8 m/s at
%   -59.5 deg latitude on the 2020-01-01 GRACE-C arc, against a nominal
%   omega*r*cos(lat) of 255.9 m/s; the difference is the polar motion and the
%   length-of-day term that a nominal rate omits. The model error of "corotation"
%   is the rigid co-rotation assumption itself: real thermospheric winds reach
%   40-80 m/s RMS in quiet conditions and several hundred m/s in storms.
%   "hwm14" reproduces NRL's published check output (Check/gfortran.txt, printed to
%   1 mm/s) through this unit; the model is single precision (NRL quotes +-1 mm/s
%   between compilers) and is a climatology: an individual day's wind departs from
%   it by an amount this toolbox has not measured.
%
% ERRORS
%   NonGravForces:wind:inputShape      r or v not 3xK with K matching the epoch
%   NonGravForces:wind:inputEop        "corotation" or "hwm14" without 'Eop'
%   NonGravForces:wind:inputAp         "hwm14" disturbed without an ap source, both sources, a
%                                      negative or malformed Ap, or an ap input to another model
%   NonGravForces:wind:depMissing      "hwm14" MEX not built
%   NonGravForces:wind:inputRange      "hwm14" at a position below the ellipsoid (the model does not return there)
%   NonGravForces:hwm14Gateway:*       the gateway's refusals, unwrapped
%   NonGravForces:wind:numericNaN      NaN/Inf in r, v or the supplied wind
%   referenceSystem:*                  propagated unwrapped
%
% REFERENCES
%   Doornbos E. (2011) Thermospheric density and wind determination from satellite
%   dynamics, Sect. 2.3; Drob D.P. et al. (2015) An update to the Horizontal Wind
%   Model (HWM): The quiet time thermosphere, Earth and Space Science 2:301;
%   IERS Conventions (2010) Ch. 5 through referenceSystem.
%
% EXAMPLE
%   eop = referenceSystem.readEop(file);
%   [ep, r, v] = NonGravForces.io.readGNV1B(gnv, 'Eop', eop);
%   [vrel, info] = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
%   vecnorm(v - vrel)          % 0 .. 470 m/s over an orbit
%
% See also NonGravForces.kernel.dragPanelForce, NonGravForces.env.density, referenceSystem.rotationItrsToGcrs.
%
% KEYWORDS wind co-rotation relative velocity thermosphere drag HWM14 DWM07
% Created 2026-09-08 by Claude Opus 5 (Anthropic); hwm14 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    model (1,1) string {mustBeMember(model, ["none","corotation","hwm14"])}
    epoch (1,1) struct
    r_gcrs double
    v_gcrs double
    options.Eop = []
    options.Wind_m_s double = []
    options.SpaceWeather = []
    options.Ap double = []
    options.Disturbed (1,1) logical = true
end
unit = "wind";
if model ~= "hwm14" && (~isempty(options.SpaceWeather) || ~isempty(options.Ap))
    NonGravForces.fail(unit, "inputAp", "'SpaceWeather' and 'Ap' drive hwm14 only -> remove them for model %s", model);
end
K = numel(epoch.sec);
if ~(ismatrix(r_gcrs) && isequal(size(r_gcrs), [3 K]) && isequal(size(v_gcrs), [3 K]))
    NonGravForces.fail(unit, "inputShape", "r_gcrs and v_gcrs must be 3xK with K = %d epochs -> got %s and %s", K, mat2str(size(r_gcrs)), mat2str(size(v_gcrs)));
end
if any(~isfinite(r_gcrs(:))) || any(~isfinite(v_gcrs(:)))
    NonGravForces.fail(unit, "numericNaN", "r_gcrs or v_gcrs contain NaN/Inf -> screen the state");
end
w = options.Wind_m_s;
if ~isempty(w)
    if ~isequal(size(w), [3 K])
        NonGravForces.fail(unit, "inputShape", "Wind_m_s must be 3xK with K = %d -> got %s", K, mat2str(size(w)));
    end
    if any(~isfinite(w(:)))
        NonGravForces.fail(unit, "numericNaN", "Wind_m_s contains NaN/Inf -> screen the wind field");
    end
end
corot = zeros(3, K);
info.hwm14 = [];
hw = zeros(3, K);
if model ~= "none"
    if isempty(options.Eop)
        NonGravForces.fail(unit, "inputEop", "%s needs 'Eop' -> pass referenceSystem.readEop(file) or the file path; the Earth rotation is not assumed", model);
    end
    [R, Rdot] = referenceSystem.rotationItrsToGcrs(epoch, 'Eop', options.Eop);
    rItrs = reshape(pagemtimes(pagetranspose(R), reshape(r_gcrs, 3, 1, K)), 3, K);
    corot = reshape(pagemtimes(Rdot, reshape(rItrs, 3, 1, K)), 3, K);   % = omega x r, carried by the dependency
    if model == "hwm14"
        [hw, info.hwm14] = localHwm14(unit, epoch, rItrs, R, options, K);
    end
end
vrel = v_gcrs - corot - hw;
if ~isempty(w)
    vrel = vrel - w;
end
info.model = model;
info.corotation_m_s = corot;
info.speed_m_s = vecnorm(vrel);
info.windGiven = ~isempty(w);
info.K = K;
info.provenance = NonGravForces.provenance("env.wind", 'Options', struct('model', model, 'windGiven', ~isempty(w), 'disturbed', options.Disturbed));
end

function [wG, hi] = localHwm14(unit, epoch, rItrs, R, options, K)
%LOCALHWM14 HWM14 (+ DWM07) at the geodetic positions, rotated to the GCRS.
ap = [];
if options.Disturbed
    if ~isempty(options.SpaceWeather) && ~isempty(options.Ap)
        NonGravForces.fail(unit, "inputAp", "give 'SpaceWeather' or 'Ap', not both -> one source for the 3-hour ap");
    end
    if ~isempty(options.Ap)
        ap = double(options.Ap(:)');
        if isscalar(ap), ap = repmat(ap, 1, K); end
        if numel(ap) ~= K || any(~isfinite(ap)) || any(ap < 0)
            NonGravForces.fail(unit, "inputAp", "'Ap' must be 1x%d (or scalar), finite and >= 0 -> the 3-hour ap in units of 2 nT", K);
        end
    elseif isstruct(options.SpaceWeather) && isfield(options.SpaceWeather, 'format') && options.SpaceWeather.format == "celestrak"
        q = options.SpaceWeather.query(epoch);
        ap = q.apHistory(2, :);                                     % the current 3-hour ap
    else
        NonGravForces.fail(unit, "inputAp", "hwm14 with 'Disturbed' true needs the 3-hour ap -> pass 'SpaceWeather', io.spaceWeather(SW-All.csv), or 'Ap'; or 'Disturbed', false for quiet-time winds");
    end
elseif ~isempty(options.SpaceWeather) || ~isempty(options.Ap)
    NonGravForces.fail(unit, "inputAp", "'Disturbed' is false, so no ap is read -> remove 'SpaceWeather' and 'Ap'");
end
here = fullfile(fileparts(mfilename('fullpath')), 'private');
if ~isfile(fullfile(here, "hwm14Gateway." + mexext))
    NonGravForces.fail(unit, "depMissing", "the HWM14 MEX is not built -> run setupNonGravForces (fetches the NRL files) and NonGravForces.buildMex once on this machine");
end
lla = referenceSystem.cartesianToGeodetic(rItrs', 'Ellipsoid', "WGS84")';   % [lat; lon; h], rad rad m
if any(lla(3, :) < 0)
    % HWM14 never returns for a negative altitude (measured 2026-09-15); the gateway refuses too.
    NonGravForces.fail(unit, "inputRange", "the position at column %d is %.3f m below the ellipsoid -> hwm14 is defined above the surface only", find(lla(3, :) < 0, 1), -min(lla(3, :)));
end
utc = referenceSystem.convertTimeScale(epoch, "UTC");
dt = referenceSystem.epochTo(utc, "datetime");
iyd = mod(double(year(dt(:)')), 100) * 1000 + double(day(dt(:)', 'dayofyear'));
ap2 = -ones(1, K);                                                % negative: quiet-time winds only
if options.Disturbed, ap2 = ap; end
wNE = hwm14Gateway(here, iyd, double(utc.sec(:)'), lla(3, :) / 1000, rad2deg(lla(1, :)), rad2deg(lla(2, :)), ap2);
sl = sin(lla(1, :)); cl = cos(lla(1, :)); so = sin(lla(2, :)); co = cos(lla(2, :));
north = [-sl .* co; -sl .* so; cl];                                % local geodetic north, ITRS
east = [-so; co; zeros(1, K)];
wI = wNE(1, :) .* north + wNE(2, :) .* east;
wG = reshape(pagemtimes(R, reshape(wI, 3, 1, K)), 3, K);
hi = struct('northEast_m_s', wNE, 'gcrs_m_s', wG, 'ap', ap, 'disturbed', options.Disturbed, 'geodetic', lla);
end
