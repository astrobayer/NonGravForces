function tiles = earthRadiance(backend, epoch, grid, r)
%EARTHRADIANCE Earth albedo and longwave emissivity on a cap of tiles around the sub-satellite point (Knocke), or as a point source.
%
%   tiles = NonGravForces.env.earthRadiance("knocke", epoch, [12 24], r)
%   tiles = NonGravForces.env.earthRadiance("pointSource", epoch, [], r)
%
%   Produces the radiative surface that NonGravForces.env.radiationSources
%   turns into pressure sources (principle P3: sources, not effects). The
%   visible cap of a spherical Earth (radius ellipsoidParameters("WGS84").a,
%   central half-angle acos(R/|r|)) is tiled in RINGS OF EQUAL CENTRAL ANGLE
%   around the sub-satellite point and SECTORS of equal azimuth, as Knocke
%   (1988) does: the near-field weight 1/rho^2 varies by orders of magnitude
%   across the cap at LEO, and a fixed global latitude-longitude grid of the
%   same tile count misses the flux by half (measured 2026-09-07: -55 % for
%   12x24 global tiles at 490 km, against -1e-3 for the cap grid below).
%
%   knocke  Knocke, Ries & Tapley (1988) Table 1, zonal with an annual term:
%             a(phi,t) = a0 + a1 cos(w(t - t0)) P1(sin phi) + a2 P2(sin phi)
%             e(phi,t) = e0 + e1 cos(w(t - t0)) P1(sin phi) + e2 P2(sin phi)
%           a0 = 0.34, a1 = 0.10, a2 = 0.29, e0 = 0.68, e1 = -0.07, e2 = -0.18,
%           w = 2 pi / 365.25 d, t0 = 22 DECEMBER (1981-12-22 00:00 UTC here;
%           the C.1 pitfall was a toolbox that counted from 1 January). phi is
%           the geocentric latitude of the tile centre in the frame of r (the
%           GCRS): the model is zonal, so no Earth rotation and no EOP enter.
%   pointSource  the whole Earth as one Lambertian sphere with a0 and e0; no
%           tiles (radiationSources applies the far-field phase law).
%   ceres   refuses with notImplemented (LP-3b, WP-27).
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   backend  1x1   string  -      "knocke" | "pointSource" | "ceres"
%   epoch    1x1   struct  -      referenceSystem epoch struct, K entries, any scale
%   grid     1x2   double  -      [nRings nSectors], positive integers (Options
%                                 erp.grid; ignored for pointSource, pass []).
%   r        3xK   double  m      Satellite position (Earth-centred, GCRS); the
%                                 cap is centred on r/|r|.
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   tiles    1x1   struct  -      .backend, .grid, .radius_m, .M (tiles, 0 for
%                                 pointSource), .n (3xMxK unit normals of the
%                                 tile centres, GCRS), .area_m2 (MxK), .lat_rad
%                                 (MxK geocentric latitude of the centres),
%                                 .capAngle_rad (1xK), .albedo (MxK; 1xK for
%                                 pointSource), .emissivity (same),
%                                 .daysSinceT0 (1xK), .coefficients, .provenance.
%
% ACCURACY
%   Knocke's zonal fit to Nimbus-7 data is itself good to a few percent of the
%   flux (~0.2 nm/s^2 in acceleration, HANDOVER C.3). Tile values are the
%   formula at the tile-centre latitude (exact to eps; tested against an
%   independent restatement of Table 1). Tile areas are exact spherical
%   ring-sector areas (they sum to the cap area to eps). The tiling error of
%   the integrated flux is measured in radiationSources' convergence test
%   against tools/pyxcheck/knocke_check.py (0.25 deg global integration) and
%   recorded per run for [6 12], [12 24], [24 48].
%
% ERRORS
%   NonGravForces:earthRadiance:inputGrid        grid not two positive integers (knocke)
%   NonGravForces:earthRadiance:inputShape       r not 3xK or epoch not an epoch struct
%   NonGravForces:earthRadiance:numericNaN       NaN/Inf in r
%   NonGravForces:earthRadiance:inputRange       |r| <= R (inside the Earth)
%   NonGravForces:earthRadiance:notImplemented   backend "ceres" (WP-27)
%   referenceSystem:*                            propagated unwrapped
%
% REFERENCES
%   Knocke P.C., Ries J.C., Tapley B.D. (1988) Earth radiation pressure effects on
%   satellites, AIAA/AAS Astrodynamics Conf., AIAA-88-4292, Table 1, eqs 3-4, sect. III.
%
% EXAMPLE
%   ep = referenceSystem.mjdToEpoch(58849.5, "UTC"); r = [6.88e6; 0; 0];
%   t = NonGravForces.env.earthRadiance("knocke", ep, [12 24], r);
%   sum(t.area_m2) / (2*pi*t.radius_m^2*(1 - cos(t.capAngle_rad)))     % 1
%
% See also NonGravForces.env.radiationSources, NonGravForces.Options.
%
% KEYWORDS albedo emissivity Knocke Earth radiation tiles cap
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    backend (1,1) string {mustBeMember(backend, ["pointSource","knocke","ceres"])}
    epoch (1,1) struct
    grid double
    r double
end
unit = "earthRadiance";
if ~all(isfield(epoch, ["day","sec","scale"]))
    NonGravForces.fail(unit, "inputShape", "epoch must be a referenceSystem epoch struct -> build one with mjdToEpoch or daySecondToEpoch");
end
K = numel(epoch.day);
if ~(ismatrix(r) && isequal(size(r), [3 K]))
    NonGravForces.fail(unit, "inputShape", "r must be 3x%d for %d epochs, got %s -> one column per epoch", K, K, mat2str(size(r)));
end
if any(~isfinite(r(:)))
    NonGravForces.fail(unit, "numericNaN", "r contains NaN/Inf -> screen the positions");
end
if backend == "ceres"
    tiles = [];
    NonGravForces.fail(unit, "notImplemented", "backend ceres is specified in LP-3b (WP-27, needs the SYN1deg variable names, O-3) -> use knocke");
end
co = localCoefficients();
utc = referenceSystem.convertTimeScale(epoch, "UTC");
t0 = referenceSystem.datetimeToEpoch(datetime(1981, 12, 22), "UTC");
days = (double(utc.day(:)') - double(t0.day)) + double(utc.sec(:)') / 86400;
seasonal = cos(2 * pi / co.period_days * days);            % 1xK
ell = referenceSystem.ellipsoidParameters("WGS84");
R = ell.a;
nr = vecnorm(r);
if any(nr <= R)
    NonGravForces.fail(unit, "inputRange", "|r| = %.0f m is not above the Earth's surface (R = %.0f m) at column %d -> check the position units (metres, GCRS)", nr(find(nr <= R, 1)), R, find(nr <= R, 1));
end
tiles.backend = backend;
tiles.radius_m = R;
tiles.daysSinceT0 = days;
tiles.coefficients = co;
tiles.capAngle_rad = acos(R ./ nr);
switch backend
    case "pointSource"
        tiles.grid = [0 0]; tiles.M = 0;
        tiles.n = zeros(3, 0, K); tiles.area_m2 = zeros(0, K); tiles.lat_rad = zeros(0, K);
        tiles.albedo = repmat(co.a0, 1, K);
        tiles.emissivity = repmat(co.e0, 1, K);
    case "knocke"
        if ~(isequal(size(grid), [1 2]) && all(isfinite(grid)) && all(grid == round(grid)) && all(grid >= 1))
            NonGravForces.fail(unit, "inputGrid", "grid must be [nRings nSectors] positive integers, got %s -> e.g. [12 24]", mat2str(grid));
        end
        nR = grid(1); nS = grid(2); M = nR * nS;
        alphaMax = tiles.capAngle_rad;                                    % 1xK
        edge = (0:nR)' / nR;                                              % (nR+1)x1 fractions
        aLo = edge(1:end-1) * alphaMax;  aHi = edge(2:end) * alphaMax;    % nRxK
        aC = 0.5 * (aLo + aHi);
        ringArea = R^2 * (2 * pi / nS) * (cos(aLo) - cos(aHi));          % nRxK
        phi = (2 * pi / nS) * ((1:nS)' - 0.5);                            % nSx1 azimuths
        % orthonormal basis of the tangent plane at the sub-satellite point, per epoch
        rhat = r ./ nr;                                                   % 3xK
        zAxis = repmat([0; 0; 1], 1, K);
        e1 = cross(zAxis, rhat, 1);
        nearPole = vecnorm(e1) < 1e-9;
        e1(:, nearPole) = repmat([1; 0; 0], 1, nnz(nearPole));
        e1 = e1 ./ vecnorm(e1);
        e2 = cross(rhat, e1, 1);
        % tile centres: ring-major, sector-minor ordering (tile (i,j) -> index (i-1)*nS + j)
        [AC, PH] = deal(repmat(reshape(aC, nR, 1, K), 1, nS, 1), repmat(phi', nR, 1, K));   % nRxnSxK
        cosA = reshape(cos(AC), 1, M, K); sinA = reshape(sin(AC), 1, M, K);
        cosP = reshape(cos(PH), 1, M, K); sinP = reshape(sin(PH), 1, M, K);
        n = cosA .* reshape(rhat, 3, 1, K) + sinA .* (cosP .* reshape(e1, 3, 1, K) + sinP .* reshape(e2, 3, 1, K));
        tiles.grid = grid; tiles.M = M;
        tiles.n = n;
        tiles.area_m2 = reshape(repmat(reshape(ringArea, nR, 1, K), 1, nS, 1), M, K);
        x = reshape(n(3, :, :), M, K);                                    % sin(latitude) of each tile centre
        tiles.lat_rad = asin(max(-1, min(1, x)));
        P2 = 0.5 * (3 * x.^2 - 1);
        tiles.albedo     = co.a0 + co.a1 * (x .* seasonal) + co.a2 * P2;     % MxK
        tiles.emissivity = co.e0 + co.e1 * (x .* seasonal) + co.e2 * P2;
end
tiles.provenance = NonGravForces.provenance("env.earthRadiance", 'Options', struct('backend', backend, 'grid', grid));
end

function co = localCoefficients()
% Knocke, Ries & Tapley (1988) Table 1. The test compares these against
% tools/pyxcheck/knocke_coefficients.csv, an independent restatement.
co = struct('a0', 0.34, 'a1', 0.10, 'a2', 0.29, 'e0', 0.68, 'e1', -0.07, 'e2', -0.18, ...
            'period_days', 365.25, 't0', "1981-12-22");
end
