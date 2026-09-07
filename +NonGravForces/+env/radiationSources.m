function [src, info] = radiationSources(epoch, r, q_body2gcrs, options, nv)
%RADIATIONSOURCES Source list {dir, P, band} in the body frame: Sun (with shadow) plus Earth albedo and longwave tiles.
%
%   [src, info] = NonGravForces.env.radiationSources(epoch, r, q, options)
%   [src, info] = NonGravForces.env.radiationSources(epoch, r, q, options, 'Ephemeris', eph)
%
%   The one place where solar and Earth radiation become SOURCES (principle
%   P3): the kernel later computes the force on each panel from every source
%   without knowing which effect it is. For K epochs the list has a FIXED
%   number of sources M = 1 + 2 * tiles: the Sun, one albedo (visible band)
%   source per tile and one longwave (infrared band) source per tile;
%   invisible or unlit tiles carry pressure 0, so everything is a plain
%   3xMxK / MxK array.
%
%   SUN. Position from referenceSystem.geocentricBodyState ("sofa" analytic
%   for Options.sun.source "rs.analytic"; "mice" with a furnished ephemeris
%   for "rs.de440", which must be passed as 'Ephemeris'). Pressure
%   P = tsi/c (au/|r_sun - r|)^2 nu with nu from eclipseFactor (Options
%   srp.eclipse; the lunar factor when srp.lunarEclipse is true).
%
%   EARTH (Options.erp.backend). knocke: earthRadiance's cap tiles (erp.grid =
%   [rings sectors] around the sub-satellite point); for
%   tile j with normal n_j, area A_j, at distance rho from the satellite with
%   outgoing cosine mu = n_j . rho_hat (visible when mu > 0) and solar cosine
%   mu_s = n_j . s_hat (lit when mu_s > 0):
%       P_albedo = a_j S mu_s mu A_j / (pi rho^2 c),   P_ir = e_j (S/4) mu A_j / (pi rho^2 c)
%   with S = tsi (au/|r_sun|)^2. pointSource: the far-field Lambert sphere,
%   P_albedo = (2/3) a0 S (R/|r|)^2 psi(theta) / c with the Lambert phase law
%   psi = (sin theta + (pi - theta) cos theta)/pi, theta the Sun-Earth-satellite
%   angle, and P_ir = e0 (S/4) (R/|r|)^2 / c; both directed at the Earth's
%   centre. none: no Earth sources.
%
%   DIRECTIONS point FROM THE SATELLITE TO THE SOURCE, in the body frame:
%   dir_body = R(q)' dir_gcrs with R from referenceSystem.rotationFromQuaternion
%   of the BODY->GCRS quaternion (io.readSCA1B's output).
%
% INPUT ARGUMENTS
%   Name          Size  Type                    Units  Description
%   ----          ----  ----                    -----  -----------
%   epoch         1x1   struct                  -      referenceSystem epoch struct, K entries
%   r             3xK   double                  m      Satellite position, GCRS
%   q_body2gcrs   4xK   double                  -      Unit quaternions BODY->GCRS, scalar first
%   options       1x1   NonGravForces.Options   -      srp.*, erp.*, sun.source are read
%
% NAME-VALUE ARGUMENTS
%   Name          Size  Type    Units  Description
%   ----          ----  ----    -----  -----------
%   Ephemeris     1x1   struct  -      referenceSystem.readEphemeris result; required
%                                     for sun.source "rs.de440". Default [].
%   Sun           3xK   double  m      Caller-supplied Earth -> Sun vector (GCRS) that
%                                     REPLACES the ephemeris (info.sunSource "caller"):
%                                     for tests and external ephemerides. Default [].
%   Moon          3xK   double  m      Caller-supplied Earth -> Moon vector, used with
%                                     'Sun' when srp.lunarEclipse is true. Default [].
%
% OUTPUT ARGUMENTS
%   Name     Size    Type    Units  Description
%   ----     ----    ----    -----  -----------
%   src      1x1     struct  -      .dir (3xMxK, body, unit, satellite -> source),
%                                   .P (MxK, Pa), .band (1xM string "vis"|"ir"),
%                                   .nu (1xK Sun shadow factor), .M, .K,
%                                   .frame "body", .meta (see info).
%   info     1x1     struct  -      .sunSource (Options value or "caller"), .sunDistance_m (1xK), .sunUnit_gcrs
%                                   (3xK), .eclipseModel, .lunarEclipse,
%                                   .erpBackend, .grid, .tiles (earthRadiance
%                                   struct or []), .tsi_W_m2, .earthRadius_m,
%                                   .speedOfLight_m_s, .provenance.
%
% ACCURACY
%   Solar irradiance uses the satellite-Sun distance for the Sun source and the
%   Earth-Sun distance for the tiles (difference 4e-5 relative, below the 1e-3
%   uncertainty of tsi itself). The tile sum converges to the continuous cap
%   integral as the grid is refined (measured 2026-09-07 against
%   tools/pyxcheck/knocke_check.py: see tests/tNgfEnv.m for the numbers per
%   grid: about 1e-3 at [12 24], the standard tier, and 3e-4 at [24 48], the precise tier). The
%   point source is a far-field model: at 490 km it differs from the tile sum
%   by the finite-distance geometry (measured and recorded per run).
%   geocentricBodyState's analytic Sun is accurate to 8e3 m RMS (1e-7 of the
%   distance): irrelevant here.
%
% ERRORS
%   NonGravForces:radiationSources:inputShape      r not 3xK, q not 4xK, or K differs from epoch
%   NonGravForces:radiationSources:numericNaN      NaN/Inf in r or q
%   NonGravForces:radiationSources:depKernel       sun.source rs.de440 without 'Ephemeris'
%   NonGravForces:eclipseFactor:*, NonGravForces:earthRadiance:*   from the callees
%   referenceSystem:*                              propagated unwrapped
%
% REFERENCES
%   Knocke, Ries, Tapley (1988) AIAA-88-4292; Montenbruck & Gill (2000) 3.4;
%   HANDOVER D.4 (source list), C.3 (approximations).
%
% EXAMPLE
%   ep = referenceSystem.mjdToEpoch(58849.5, "GPS"); r = [6.88e6; 0; 0]; q = [1; 0; 0; 0];
%   src = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard"));
%   size(src.dir), src.band(1), src.P(1)          % [3 577 1], "vis", ~4.5e-6 Pa
%
% See also NonGravForces.env.eclipseFactor, NonGravForces.env.earthRadiance, NonGravForces.Options.
%
% KEYWORDS radiation sources SRP albedo longwave source list body frame
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    epoch (1,1) struct
    r double
    q_body2gcrs double
    options (1,1) NonGravForces.Options
    nv.Ephemeris = []
    nv.Sun = []
    nv.Moon = []
end
unit = "radiationSources";
if ~all(isfield(epoch, ["day","sec","scale"]))
    NonGravForces.fail(unit, "inputShape", "epoch must be a referenceSystem epoch struct -> build one with mjdToEpoch or daySecondToEpoch");
end
K = numel(epoch.day);
if ~(ismatrix(r) && isequal(size(r), [3 K]))
    NonGravForces.fail(unit, "inputShape", "r must be 3x%d for %d epochs, got %s -> one column per epoch", K, K, mat2str(size(r)));
end
if ~(ismatrix(q_body2gcrs) && isequal(size(q_body2gcrs), [4 K]))
    NonGravForces.fail(unit, "inputShape", "q_body2gcrs must be 4x%d, got %s -> io.readSCA1B or geom.AttitudeProvider output", K, mat2str(size(q_body2gcrs)));
end
if any(~isfinite(r(:))) || any(~isfinite(q_body2gcrs(:)))
    NonGravForces.fail(unit, "numericNaN", "r or q contains NaN/Inf -> screen the state first");
end
c = NonGravForces.const;
cLight = referenceSystem.iersConstants().speedOfLight;
% --- Sun (and Moon) positions
sunSource = options.sun.source;
if ~isempty(nv.Sun)
    if ~(ismatrix(nv.Sun) && isequal(size(nv.Sun), [3 K])) || any(~isfinite(nv.Sun(:)))
        NonGravForces.fail(unit, "inputShape", "Sun must be a finite 3x%d array -> one Earth->Sun vector per epoch", K);
    end
    rs = nv.Sun;
    if options.srp.lunarEclipse
        if ~(ismatrix(nv.Moon) && isequal(size(nv.Moon), [3 K]))
            NonGravForces.fail(unit, "inputShape", "srp.lunarEclipse with a caller Sun needs 'Moon' 3x%d as well -> pass both or neither", K);
        end
        rm = nv.Moon;
    else
        rm = [];
    end
    sunSource = "caller";
else
switch options.sun.source
    case "rs.analytic"
        rs = referenceSystem.geocentricBodyState("sun", epoch);
        if options.srp.lunarEclipse, rm = referenceSystem.geocentricBodyState("moon", epoch); else, rm = []; end
    case "rs.de440"
        if ~(isstruct(nv.Ephemeris) && ~isempty(nv.Ephemeris))
            NonGravForces.fail(unit, "depKernel", "sun.source rs.de440 needs furnished kernels -> referenceSystem.downloadKernels once, then pass 'Ephemeris', referenceSystem.readEphemeris(files)");
        end
        rs = referenceSystem.geocentricBodyState("sun", epoch, 'Source', "mice", 'Ephemeris', nv.Ephemeris);
        if options.srp.lunarEclipse, rm = referenceSystem.geocentricBodyState("moon", epoch, 'Source', "mice", 'Ephemeris', nv.Ephemeris); else, rm = []; end
end
end
% --- shadow and solar pressure
[nu, ecl] = NonGravForces.env.eclipseFactor(r, rs, 'Model', options.srp.eclipse, 'Moon', rm);
d = rs - r; nd = vecnorm(d);
sunUnitSat = d ./ nd;                                              % 3xK, satellite -> Sun
Psun = options.srp.tsi / cLight * (c.au ./ nd).^2 .* nu;            % 1xK
S_earth = options.srp.tsi * (c.au ./ vecnorm(rs)).^2;               % 1xK irradiance at the Earth
sHat = rs ./ vecnorm(rs);                                           % 3xK
% --- Earth sources
R = ecl.radiusEarth_m;
switch options.erp.backend
    case "none"
        tiles = []; M = 0;
        dirE = zeros(3, 0, K); PE = zeros(0, K); bandE = string.empty(1, 0);
    case "pointSource"
        tiles = NonGravForces.env.earthRadiance("pointSource", epoch, [], r);
        nr = vecnorm(r); rhat = r ./ nr;
        cosT = max(-1, min(1, sum(rhat .* sHat, 1)));
        theta = acos(cosT);
        psi = (sin(theta) + (pi - theta) .* cos(theta)) / pi;
        Ealb = (2/3) * tiles.albedo .* S_earth .* (R ./ nr).^2 .* psi;
        Eir  = tiles.emissivity .* (S_earth / 4) .* (R ./ nr).^2;
        M = 2;
        dirE = repmat(reshape(-rhat, 3, 1, K), 1, 2, 1);
        PE = [Ealb; Eir] / cLight;
        bandE = ["vis", "ir"];
    case "knocke"
        tiles = NonGravForces.env.earthRadiance("knocke", epoch, options.erp.grid, r);
        Mt = tiles.M;
        p = R * tiles.n;                                            % 3xMtxK tile centres
        rho = reshape(r, 3, 1, K) - p;                              % 3xMtxK tile -> satellite
        dist2 = sum(rho.^2, 1);                                     % 1xMtxK
        dist = sqrt(dist2);
        rhoHat = rho ./ dist;
        mu = max(0, sum(tiles.n .* rhoHat, 1));                     % 1xMtxK outgoing cosine
        muS = max(0, reshape(sum(tiles.n .* reshape(sHat, 3, 1, K), 1), Mt, K));   % MtxK solar cosine
        g = reshape(mu, Mt, K) .* tiles.area_m2 ./ reshape(pi * dist2, Mt, K);    % MtxK
        Ealb = tiles.albedo .* S_earth .* muS .* g;
        Eir  = tiles.emissivity .* (S_earth / 4) .* g;
        M = 2 * Mt;
        dirE = cat(2, -rhoHat, -rhoHat);                            % satellite -> tile, twice (vis, ir)
        PE = [Ealb; Eir] / cLight;
        bandE = [repmat("vis", 1, Mt), repmat("ir", 1, Mt)];
    case "ceres"
        tiles = NonGravForces.env.earthRadiance("ceres", epoch, options.erp.grid, r);   % refuses (WP-27)
end
% --- to the body frame, assemble
Rq = referenceSystem.rotationFromQuaternion(q_body2gcrs);            % 3x3xK BODY->GCRS
dirG = cat(2, reshape(sunUnitSat, 3, 1, K), dirE);                   % 3x(1+M)xK in the GCRS
src.dir = pagemtimes(Rq, 'transpose', dirG, 'none');
src.P = [Psun; PE];
src.band = ["vis", bandE];
src.nu = nu;
src.M = 1 + M;
src.K = K;
src.frame = "body";
info.sunSource = sunSource;
info.sunDistance_m = nd;
info.sunUnit_gcrs = sHat;
info.eclipseModel = options.srp.eclipse;
info.lunarEclipse = options.srp.lunarEclipse;
info.erpBackend = options.erp.backend;
info.grid = options.erp.grid;
info.tiles = tiles;
info.tsi_W_m2 = options.srp.tsi;
info.earthRadius_m = R;
info.speedOfLight_m_s = cLight;
info.provenance = NonGravForces.provenance("env.radiationSources", 'Options', options.toStruct());
src.meta = info;
end
