function [nu, info] = eclipseFactor(r, r_sun, options)
%ECLIPSEFACTOR Shadow function nu in [0,1]: cylindrical, conical or SOLAARS-CF Earth shadow, optional lunar shadow.
%
%   nu = NonGravForces.env.eclipseFactor(r, r_sun)                       % conical
%   [nu, info] = NonGravForces.env.eclipseFactor(r, r_sun, 'Model', "cylindrical")
%   [nu, info] = NonGravForces.env.eclipseFactor(r, r_sun, 'Model', "solaarsCF", 'EclipticNormal', e)
%   [nu, info] = NonGravForces.env.eclipseFactor(r, r_sun, 'Moon', r_moon)  % Earth x Moon factor
%
%   Fraction of the solar radiation reaching the satellite, K-vectorised. For
%   the geometric models the occulting bodies are spheres: the Earth of radius
%   referenceSystem.ellipsoidParameters("WGS84").a (the equatorial radius, the
%   conservative, widest shadow) and, with 'Moon', the Moon of radius
%   NonGravForces.const.Rmoon. The Earth and Moon factors multiply (the lunar
%   factor is always geometric-conical).
%
%   MODELS
%     cylindrical  nu = 0 when the satellite is behind the body (r - r_B has a
%                  negative component along the Sun direction) and inside the
%                  cylinder of radius R_B; else 1. C^-1 at the edge.
%     conical      Montenbruck & Gill (2000) 3.4.2: with the apparent radii
%                  a = asin(R_S/|r_S - r|), b = asin(R_B/|r_B - r|) and the
%                  angular separation c: full light for c >= a + b, umbra 0 for
%                  c <= b - a, annular 1 - b^2/a^2 for c <= a - b, else
%                  x = (c^2 + a^2 - b^2)/(2c), y = sqrt(a^2 - x^2),
%                  A = a^2 acos(x/a) + b^2 acos((c-x)/b) - c y, nu = 1 - A/(pi a^2),
%                  evaluated with the arcs as atan2(y, x) and atan2(y, c-x) and the
%                  separation c as atan2(|d x s|, d.s) (see ACCURACY).
%                  Continuous, C^0 at the umbra edge (HANDOVER C.3).
%     solaarsCF    Robertson (2015) Sect. 3.7: the SOLAARS curve fit, a sum of
%                  three hyperbolic-tangent sigmoids in the distance rE of the
%                  satellite from the Sun-Earth line, with coefficients that
%                  depend on the distance rR behind the terminator plane
%                  (Table 3.2), fitted to a full ray-tracing model of the
%                  penumbra WITH Earth oblateness (rE's component along the
%                  ecliptic normal, the projected spin axis, stretched by a/b
%                  onto the fit's sphere of radius a; until 2026-10-31 the
%                  in-ecliptic part was, an inverted flattening that widened
%                  the shadow ~8 s per edge), atmospheric refraction, absorption and
%                  scattering. Smooth (C^inf). Needs 'EclipticNormal'. The day
%                  side (rR <= 0) is full light; rR is clamped to the fit
%                  domain 2-45 Mm and the clamp count is reported.
%
% INPUT ARGUMENTS
%   Name     Size      Type    Units  Description
%   ----     ----      ----    -----  -----------
%   r        3xK       double  m      Satellite position, Earth-centred (GCRS).
%   r_sun    3xK|3x1   double  m      Earth -> Sun vector (same for all K if 3x1).
%
% NAME-VALUE ARGUMENTS
%   Name            Size      Type    Units  Description
%   ----            ----      ----    -----  -----------
%   Model           1x1       string  -      "conical" (default) | "cylindrical" | "solaarsCF"
%   Moon            3xK|3x1   double  m      Earth -> Moon vector; adds the lunar factor.
%                                            Default [] (no lunar shadow).
%   EclipticNormal  3xK|3x1   double  -      Unit normal of the ecliptic in the frame of
%                                            r (e.g. r_sun x v_sun normalised, as
%                                            radiationSources measures it). Required
%                                            for solaarsCF, ignored otherwise. Default [].
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   nu       1xK   double  -      Fraction of the solar radiation, in [0, 1].
%   info     1x1   struct  -      .model, .earth (1xK Earth factor), .moon (1xK
%                                 or []), .radiusEarth_m, .radiusMoon_m,
%                                 .radiusSun_m, .solaars (struct rR_Mm, rE_Mm
%                                 1xK, oblatenessScale, clampedCount; solaarsCF
%                                 only), .provenance.
%
% ACCURACY
%   The conical formula treats the two discs in the tangent plane (small
%   angles for the Sun; the Earth's disc is 1.2 rad at LEO), which is the
%   Montenbruck & Gill approximation; against a numerical integration of the
%   same disc overlap (tools/pyxcheck/eclipse_check.py) it agrees to 1e-12
%   (measured). CONDITIONING: written with acos, the Earth-disc arc
%   acos((c-x)/b) has its argument within 2e-3 of 1 and the lens area is a
%   difference of terms 1e4 times larger than itself, which lost 3e-11 in nu;
%   the atan2 forms restore full precision (tested as rotation invariance to
%   1e-12, the HANDOVER accuracy figure). What the geometric models do not
%   model: Earth oblateness and atmospheric refraction and extinction (the
%   penumbra is in reality wider and darker: Robertson 2015). The oblateness
%   is NOT small in time: the geometric models take a sphere of the EQUATORIAL
%   radius a, while on a polar orbit the limb at shadow entry and exit lies at
%   high latitude, where the radius is close to b (21 km less). At 490 km and
%   beta 0 every nu level is crossed about 7.5 s early at entry and late at exit,
%   about the duration of the penumbra itself (7.8 s); measured with an independent
%   re-implementation of "conical" in the review of 2026-11-12 (HANDOVER O-57).
%   solaarsCF models the oblateness. The Sun is a uniform disc.
%   SOLAARS-CF: the thesis quotes a curve-fit residual below 3 % of nu against
%   the full SOLAARS model (Fig. 3.30); the coefficients here were transcribed
%   from Table 3.2 by three agreeing text extractions and reproduce the
%   sigmoid centres of Fig. 3.28 (tools/pyxcheck/solaars_check.py is the
%   independent restatement); no published numeric nu value exists. Outside
%   the fit domain (rR < 2 Mm, i.e. very low orbits near the terminator, or
%   rR > 45 Mm beyond GEO) the coefficients are held at the domain edge.
%
% ERRORS
%   NonGravForces:eclipseFactor:inputShape           r not 3xK; r_sun, Moon or EclipticNormal not 3xK / 3x1
%   NonGravForces:eclipseFactor:numericNaN           NaN/Inf in a position
%   NonGravForces:eclipseFactor:inputEclipticNormal  solaarsCF without EclipticNormal, or not unit
%   referenceSystem:*                                propagated unwrapped
%
% REFERENCES
%   Montenbruck O., Gill E. (2000) Satellite Orbits, Springer, sect. 3.4.2, eqs 3.79-3.87.
%   Robertson R.V. (2015) Highly Physical Solar Radiation Pressure Modeling During Penumbra
%   Transitions, PhD thesis, Virginia Tech, sect. 3.7, eqs 3.80-3.86, Table 3.2.
%
% EXAMPLE
%   r = [-(6378137 + 490e3) * cosd(68.25); 0; (6378137 + 490e3) * sind(68.25)];
%   nu = NonGravForces.env.eclipseFactor(r, [149597870700; 0; 0])   % penumbra: 0 < nu < 1
%
% See also NonGravForces.env.radiationSources, referenceSystem.geocentricBodyState.
%
% KEYWORDS eclipse shadow penumbra umbra conical cylindrical SOLAARS
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    r double
    r_sun double
    options.Model (1,1) string {mustBeMember(options.Model, ["cylindrical","conical","solaarsCF"])} = "conical"
    options.Moon double = []
    options.EclipticNormal double = []
end
unit = "eclipseFactor";
if ~(ismatrix(r) && size(r, 1) == 3)
    NonGravForces.fail(unit, "inputShape", "r must be 3xK, got %s -> one column per epoch", mat2str(size(r)));
end
if any(~isfinite(r(:)))
    NonGravForces.fail(unit, "numericNaN", "r contains NaN/Inf at column %d -> screen the positions", find(any(~isfinite(r), 1), 1));
end
K = size(r, 2);
r_sun = localExpand(unit, r_sun, K, "r_sun");
c = NonGravForces.const;
ell = referenceSystem.ellipsoidParameters("WGS84");
info.model = options.Model;
info.radiusEarth_m = ell.a;
info.radiusSun_m = c.Rsun;
info.radiusMoon_m = c.Rmoon;
if options.Model == "solaarsCF"
    if isempty(options.EclipticNormal)
        NonGravForces.fail(unit, "inputEclipticNormal", "solaarsCF needs 'EclipticNormal' (3x1 or 3xK unit vector, e.g. r_sun x v_sun normalised) -> pass it, or use conical");
    end
    e = localExpand(unit, options.EclipticNormal, K, "EclipticNormal");
    if any(abs(vecnorm(e) - 1) > 1e-6)
        NonGravForces.fail(unit, "inputEclipticNormal", "EclipticNormal must be unit vectors (norm off by %.2e) -> normalise", max(abs(vecnorm(e) - 1)));
    end
    [info.earth, info.solaars] = localSolaars(r, r_sun, e, ell.a, ell.b);
else
    info.earth = localBody(r, r_sun, zeros(3, K), ell.a, c.Rsun, options.Model);
end
nu = info.earth;
if isempty(options.Moon)
    info.moon = [];
else
    r_moon = localExpand(unit, options.Moon, K, "Moon");
    moonModel = options.Model;
    if moonModel == "solaarsCF", moonModel = "conical"; end
    info.moon = localBody(r, r_sun, r_moon, c.Rmoon, c.Rsun, moonModel);
    nu = nu .* info.moon;
end
info.provenance = NonGravForces.provenance("env.eclipseFactor", 'Options', struct('Model', options.Model, 'moon', ~isempty(options.Moon)));
end

function v = localExpand(unit, v, K, name)
if ~(ismatrix(v) && size(v, 1) == 3 && any(size(v, 2) == [1 K]))
    NonGravForces.fail(unit, "inputShape", "%s must be 3xK (K = %d) or 3x1, got %s -> one column per epoch or one for all", name, K, mat2str(size(v)));
end
if size(v, 2) == 1, v = repmat(v, 1, K); end
if any(~isfinite(v(:)))
    NonGravForces.fail(unit, "numericNaN", "%s contains NaN/Inf -> screen the inputs", name);
end
end

function nu = localBody(r, rs, rb, Rb, Rsun, model)
% Visible fraction of the Sun with body B (centre rb, radius Rb) as the occulter.
d = rs - r;  nd = vecnorm(d);          % satellite -> Sun
s = rb - r;  ns = vecnorm(s);          % satellite -> body centre
K = size(r, 2);
nu = ones(1, K);
canOcc = ns < nd;                      % a body beyond the Sun casts no shadow here
switch model
    case "cylindrical"
        u = d ./ nd;
        along = sum(s .* u, 1);
        perp = vecnorm(s - along .* u);
        nu(canOcc & along > 0 & perp < Rb) = 0;
    case "conical"
        a = asin(min(1, Rsun ./ nd));
        b = asin(min(1, Rb ./ ns));
        % separation angle by atan2: well conditioned at the small angles of the penumbra, where
        % acos of a cosine near 1 loses ~1e-12 (measured 8e-12 in nu, 2026-09-07)
        cc = atan2(vecnorm(cross(d, s, 1)), sum(d .* s, 1));
        full = cc >= a + b;
        umbra = cc <= b - a;
        annular = ~umbra & (cc <= a - b);
        partial = canOcc & ~full & ~umbra & ~annular;
        nu(canOcc & umbra) = 0;
        nu(canOcc & annular) = 1 - (b(canOcc & annular) ./ a(canOcc & annular)).^2;
        if any(partial)
            ap = a(partial); bp = b(partial); cp = cc(partial);
            x = (cp.^2 + ap.^2 - bp.^2) ./ (2 * cp);
            y = sqrt(max(0, ap.^2 - x.^2));
            % lens area with the two arcs written as atan2 of the shared half-chord y: acos((c-x)/b) with
            % (c-x)/b within 2e-3 of 1 (Earth's disc 250x the Sun's) amplified rounding to 3e-11 in nu
            % (measured 2026-09-07); atan2(y, c-x) is well conditioned there and equal in exact arithmetic.
            A = ap.^2 .* atan2(y, x) + bp.^2 .* atan2(y, cp - x) - cp .* y;
            nu(partial) = 1 - A ./ (pi * ap.^2);
        end
        nu = min(1, max(0, nu));
end
end

function [nu, sol] = localSolaars(r, rs, e, aE, bE)
% Robertson (2015) SOLAARS-CF, Sect. 3.7: eqs 3.80-3.86 and Table 3.2. Distances in 1e6 m.
Rhat = rs ./ vecnorm(rs);                          % Earth -> Sun unit vector, 3xK
rR = -sum(r .* Rhat, 1);                           % eq. 3.80: distance behind the terminator plane [m]
rEvec = r + rR .* Rhat;                            % eq. 3.81: r - (r.R)R
rEperp = sum(e .* rEvec, 1);                       % eq. 3.82: component along the ecliptic normal
rEpar2 = max(0, sum(rEvec.^2, 1) - rEperp.^2);     % eq. 3.83 (squared): the part in the ecliptic plane
sO = aE / bE;                                      % eq. 3.84
% eq. 3.85: the oblate silhouette (semi-axes a in the ecliptic, ~b along the normal, the projected spin axis) is
% mapped onto the fit's sphere of radius a by stretching the NORMAL component by a/b. Until 2026-10-31 the
% in-ecliptic part was stretched instead -- an Earth flattened at the equator -- which widened the shadow by
% ~8 s per edge on a polar orbit at small beta and shifted the edges asymmetrically at larger beta (HANDOVER,
% "the x lead read at the shadow edges"; the 50 % point at the polar limb now coincides with conical's).
rE = sqrt(sO^2 * rEperp.^2 + rEpar2) / 1e6;        % eq. 3.85, in Mm
rRm = rR / 1e6;
lo = 2.0; hi = 45.0;                               % the fit domain of the thesis (2000-45000 km)
clamped = rRm > 0 & (rRm < lo | rRm > hi);
rRc = min(hi, max(lo, rRm));
T = localSolaarsTable();
a = zeros(8, numel(rRm));
for k = 1:8
    b = T(k, 2:5);
    switch T(k, 1)
        case 1, a(k, :) = b(1) * exp(b(2) * rRc) + b(3) * exp(b(4) * rRc);
        case 2, a(k, :) = b(1) * rRc + b(2);
        case 3, a(k, :) = b(1) * rRc.^b(2) + b(3);
    end
end
num = 1 + a(1,:) + a(2,:) + a(1,:) .* tanh(a(3,:) .* (rE - a(4,:))) + a(2,:) .* tanh(a(5,:) .* (rE - a(6,:))) + tanh(a(7,:) .* (rE - a(8,:)));
nu = num ./ (2 + 2 * a(1,:) + 2 * a(2,:));         % eq. 3.86
nu(rR <= 0) = 1;                                   % day side: the fit describes the night side only
nu = min(1, max(0, nu));
sol = struct('rR_Mm', rRm, 'rE_Mm', rE, 'oblatenessScale', sO, 'clampedCount', nnz(clamped));
end

function T = localSolaarsTable()
% Robertson (2015) Table 3.2: form (1: b1 e^{b2 rR} + b3 e^{b4 rR}; 2: b1 rR + b2; 3: b1 rR^b2 + b3), b1..b4.
% Independently restated in tools/pyxcheck/solaars_cf_coefficients.csv (tested to agree).
T = [1,  0.1715,    -0.1423,     0.01061,  -0.01443;   % a1
     2,  0.008162,   0.3401,     0,         0;         % a2
     1,  260.9,     -0.4661,     27.81,    -0.009437;  % a3
     3, -0.006119,   1.176,      6.385,     0;         % a4
     1,  87.56,     -0.09188,    19.30,    -0.01089;   % a5
     2,  0.002047,   6.409,      0,         0;         % a6
     1,  61.98,     -0.1629,     27.87,    -0.02217;   % a7
     1,  6.413,     -0.0002593, -0.01479,  -0.1318];   % a8
end
