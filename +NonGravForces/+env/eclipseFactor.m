function [nu, info] = eclipseFactor(r, r_sun, options)
%ECLIPSEFACTOR Shadow function nu in [0,1]: cylindrical or conical Earth shadow, optional lunar shadow.
%
%   nu = NonGravForces.env.eclipseFactor(r, r_sun)                       % conical
%   [nu, info] = NonGravForces.env.eclipseFactor(r, r_sun, 'Model', "cylindrical")
%   [nu, info] = NonGravForces.env.eclipseFactor(r, r_sun, 'Moon', r_moon)  % Earth x Moon factor
%
%   Fraction of the solar disc visible from the satellite, K-vectorised. The
%   occulting bodies are spheres: the Earth of radius
%   referenceSystem.ellipsoidParameters("WGS84").a (the equatorial radius, the
%   conservative, widest shadow) and, with 'Moon', the Moon of radius
%   NonGravForces.const.Rmoon. The two factors multiply.
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
%     solaarsCF    refuses with notImplemented until the SOLAARS-CF coefficient
%                  table is in hand (HANDOVER O-13).
%
% INPUT ARGUMENTS
%   Name     Size      Type    Units  Description
%   ----     ----      ----    -----  -----------
%   r        3xK       double  m      Satellite position, Earth-centred (GCRS).
%   r_sun    3xK|3x1   double  m      Earth -> Sun vector (same for all K if 3x1).
%
% NAME-VALUE ARGUMENTS
%   Name     Size      Type    Units  Description
%   ----     ----      ----    -----  -----------
%   Model    1x1       string  -      "conical" (default) | "cylindrical" | "solaarsCF"
%   Moon     3xK|3x1   double  m      Earth -> Moon vector; adds the lunar factor.
%                                     Default [] (no lunar shadow).
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   nu       1xK   double  -      Visible fraction of the Sun, in [0, 1].
%   info     1x1   struct  -      .model, .earth (1xK Earth factor), .moon (1xK
%                                 or []), .radiusEarth_m, .radiusMoon_m,
%                                 .radiusSun_m, .provenance.
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
%   1e-12, the HANDOVER accuracy figure). What it
%   does not model: Earth oblateness (up to 21 km of limb height, ~2 % of the
%   penumbra width) and atmospheric refraction and extinction (the penumbra
%   is in reality wider and darker: Robertson 2015) - both are SOLAARS-CF's
%   business. The Sun is a uniform disc (no limb darkening).
%
% ERRORS
%   NonGravForces:eclipseFactor:inputShape        r not 3xK; r_sun or Moon not 3xK / 3x1
%   NonGravForces:eclipseFactor:numericNaN        NaN/Inf in a position
%   NonGravForces:eclipseFactor:notImplemented    Model "solaarsCF" (O-13)
%   referenceSystem:*                             propagated unwrapped
%
% REFERENCES
%   Montenbruck O., Gill E. (2000) Satellite Orbits, Springer, sect. 3.4.2, eqs 3.79-3.87.
%   Robertson R.V. (2015) PhD thesis, Virginia Tech (SOLAARS).
%
% EXAMPLE
%   r = [-(6378137 + 490e3) * cosd(68.25); 0; (6378137 + 490e3) * sind(68.25)];
%   nu = NonGravForces.env.eclipseFactor(r, [149597870700; 0; 0])   % penumbra: 0 < nu < 1
%
% See also NonGravForces.env.radiationSources, referenceSystem.geocentricBodyState.
%
% KEYWORDS eclipse shadow penumbra umbra conical cylindrical
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    r double
    r_sun double
    options.Model (1,1) string {mustBeMember(options.Model, ["cylindrical","conical","solaarsCF"])} = "conical"
    options.Moon double = []
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
if options.Model == "solaarsCF"
    NonGravForces.fail(unit, "notImplemented", "solaarsCF needs the SOLAARS-CF coefficient table (HANDOVER O-13) -> use conical or cylindrical");
end
c = NonGravForces.const;
ell = referenceSystem.ellipsoidParameters("WGS84");
info.model = options.Model;
info.radiusEarth_m = ell.a;
info.radiusSun_m = c.Rsun;
info.radiusMoon_m = c.Rmoon;
info.earth = localBody(r, r_sun, zeros(3, K), ell.a, c.Rsun, options.Model);
nu = info.earth;
if isempty(options.Moon)
    info.moon = [];
else
    r_moon = localExpand(unit, options.Moon, K, "Moon");
    info.moon = localBody(r, r_sun, r_moon, c.Rmoon, c.Rsun, options.Model);
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
