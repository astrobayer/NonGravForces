function tiles = earthRadiance(backend, epoch, grid, r, nv)
%EARTHRADIANCE Earth albedo and longwave emission on a cap of tiles around the sub-satellite point (Knocke or CERES), or as a point source.
%
%   tiles = NonGravForces.env.earthRadiance("knocke", epoch, [12 24], r)
%   tiles = NonGravForces.env.earthRadiance("ceres", epoch, [24 48], r, 'Ceres', cer, 'RotationItrsToGcrs', R)
%   tiles = NonGravForces.env.earthRadiance("pointSource", epoch, [], r)
%   tiles = NonGravForces.env.earthRadiance("ceresClimatology", epoch, [24 48], r, 'Climatology', clm, 'RotationItrsToGcrs', R)
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
%   ceres   the measured monthly CERES EBAF fluxes of io.importCERES on the same
%           cap tiles, on the sphere of the product's 20 km flux reference level
%           (radius a + 20 km). Each tile centre is rotated into the ITRS with
%           'RotationItrsToGcrs' (ITRS->GCRS, 3x3xK, as rotationItrsToGcrs returns
%           it), its geocentric direction converted to the grid's GEODETIC
%           latitude, tan(phi_d) = tan(phi_c)/(1 - f)^2, and every field is
%           interpolated bilinearly between cell centres (periodic in longitude,
%           held constant poleward of the outermost centres) and linearly in time
%           between the two monthly nodes that bracket the epoch. Then
%             albedo    = SWup / SWdown   (0 where the interpolated SWdown is 0)
%             exitance  = LWup            (W m^-2, emitted as given)
%           The albedo is a monthly TOA albedo that radiationSources applies to the
%           INSTANTANEOUS insolation, which restores the day side; the longwave is
%           a monthly mean with no diurnal cycle. EBAF albedos above 1 near the
%           terminator are the product's own (its DQS says so) and are kept.
%   ceresClimatology  the GROOPS monthly CERES climatology, two grids read by
%           io.importGroopsGrid and passed together as 'Climatology':
%             .albedo    CERES_SYN1deg_albedo_climatology.dat      dimensionless
%             .longwave  CERES_SYN1deg_longwaveFlux_climatology.dat  W m^-2
%           1 degree, 2412 monthly grids bounding 2413 nodes from 1900 to 2101. Like
%           ceres it supplies a MEASURED EXITANCE rather than an emissivity, so
%           .emissivity is empty and radiationSources uses the same branch. Four
%           things differ from ceres and each is asserted rather than assumed:
%             - TIME IS PIECEWISE CONSTANT, not interpolated. The container bounds T
%               grids by T+1 nodes and each grid IS that month's mean; interpolating
%               between two would invent a value the product does not carry.
%             - THE GRID'S LATITUDES ARE GEOCENTRIC (io.importGroopsGrid's help said
%               geodetic until 2026-09-17). The grid is laid out uniformly in GEODETIC
%               latitude and the file stores the GEOCENTRIC latitude of those points,
%               so the centres are 89.4966, 88.4899, ... and their spacing runs 0.9933
%               to 1.0067. The tile latitude is therefore NOT converted here, where the
%               ceres backend does convert; the two conventions differ by up to 0.19
%               deg, a fifth of a cell.
%             - LATITUDE IS SEARCHED, longitude indexed. The latitude axis is not
%               uniform (above) and the longitude axis carries single-precision jitter
%               (1.37e-05 deg between neighbouring gaps, first centre -179.499994),
%               so longitude is indexed with the AVERAGE spacing and its uniformity is
%               asserted to 1e-3 rather than to float noise.
%             - THE REFERENCE LEVEL IS READ OUT OF THE FILE, not declared: the per-row
%               radii are the ellipsoid radius at the stored geocentric latitude plus a
%               constant, 30000.1 m for this product, to a spread of 0.168 m over 180
%               rows. Measuring it is what caught the latitude convention.
%           It is a CLIMATOLOGY: the file name says SYN1deg but there is no diurnal
%           cycle in a monthly mean, and the 3-hourly product is a different thing
%           (HANDOVER O-3, still open).
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   backend  1x1   string  -      "knocke" | "pointSource" | "ceres" | "ceresClimatology"
%   epoch    1x1   struct  -      referenceSystem epoch struct, K entries, any scale
%   grid     1x2   double  -      [nRings nSectors], positive integers (Options
%                                 erp.grid; ignored for pointSource, pass []).
%   r        3xK   double  m      Satellite position (Earth-centred, GCRS); the
%                                 cap is centred on r/|r|.
%
% NAME-VALUE ARGUMENTS
%   Name                Size   Type    Units  Description
%   ----                ----   ----    -----  -----------
%   Ceres               1x1    struct  -      io.importCERES result; required by "ceres", refused by the others
%   Climatology         1x1    struct  -      .albedo and .longwave, each an io.importGroopsGrid result;
%                                             required by "ceresClimatology", refused by the others
%   RotationItrsToGcrs  3x3xK  double  -      ITRS->GCRS rotation at the epochs; required by both gridded
%                                             backends, refused by knocke and pointSource
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   tiles    1x1   struct  -      .backend, .grid, .radius_m, .M (tiles, 0 for
%                                 pointSource), .n (3xMxK unit normals of the
%                                 tile centres, GCRS), .area_m2 (MxK), .lat_rad
%                                 (MxK geocentric latitude of the centres, GCRS;
%                                 for ceres the ITRS geodetic latitude the grid was read at,
%                                 for ceresClimatology the ITRS geocentric one),
%                                 .capAngle_rad (1xK), .albedo (MxK; 1xK for
%                                 pointSource), .emissivity (same; [] for ceres and
%                                 ceresClimatology), .exitanceLw_W_m2 (MxK, ceres and
%                                 ceresClimatology; [] otherwise),
%                                 .ceres (ceres: .product .version .months (the two
%                                 bracketing nodes per epoch, 2xK datetime) .timeWeight (1xK)),
%                                 .climatology (ceresClimatology: .referenceHeight_m measured from
%                                 the file, .gridIndex (1xK, into the grids PASSED), .month (1xK
%                                 datetime of each grid's opening node), .nodes (1x2, first and
%                                 last node of the grids passed), .gridCount),
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
%   ceres: the tile values are the interpolated product (exact to rounding for a
%   field that is linear across a cell, tested); the tile sum is checked against
%   tools/pyxcheck/ceres_check.py, an independent standard-library integration of
%   the same continuous field over the visible cap. The model error is the
%   product's (a few W m^-2 per cell and month) plus the missing diurnal cycle of
%   the longwave and of the albedo itself. TWILIGHT: EBAF's SWup includes a
%   twilight flux, 0.3-0.6 W m^-2 in polar-night cells where SWdown is 0 (all
%   7920 such cells of December 2019, measured by the oracle); those cells get
%   albedo 0 and the twilight flux is lost, and at the polar-night edge the
%   ratio exceeds 1 (11.8 in one December cell) where the instantaneous
%   insolation it multiplies is at or near zero.
%   ceresClimatology: the same interpolation error as ceres, plus the climatology's
%   own -- it is a monthly mean of many years, so it carries neither the diurnal cycle
%   nor the particular year. Measured against knocke over the reference day's cap
%   (2026-09-17, 24x48 tiles, December 2019): the tile albedo spans 0.131-0.692 against
%   knocke's smooth 0.189-0.522 and the exitance 175.7-275.3 W m^-2 against knocke's
%   implied 195.7-263.5, with cap means 0.468 against 0.416 and 211.1 against 217.6
%   W m^-2. So it is the same physics with the spatial structure restored, not a
%   different magnitude.
%
% ERRORS
%   NonGravForces:earthRadiance:inputClimatology the Climatology grids are missing, mismatched, not a
%                                                time series, not uniform in longitude, not decreasing in
%                                                latitude, not one constant height above the ellipsoid,
%                                                given to a backend that does not read them, or
%                                                ceresClimatology without a usable 'RotationItrsToGcrs'
%   NonGravForces:earthRadiance:inputGrid        grid not two positive integers (every backend but
%                                                pointSource)
%   NonGravForces:earthRadiance:inputShape       r not 3xK or epoch not an epoch struct
%   NonGravForces:earthRadiance:numericNaN       NaN/Inf in r
%   NonGravForces:earthRadiance:inputRange       |r| <= R (inside the Earth, or below the reference level
%                                                of ceres or ceresClimatology)
%   NonGravForces:earthRadiance:inputCeres       ceres without a usable 'Ceres' or 'RotationItrsToGcrs',
%                                                'Ceres' given to another backend, or
%                                                'RotationItrsToGcrs' given to knocke or pointSource
%   NonGravForces:earthRadiance:dataRange        an epoch outside the months the 'Ceres' struct holds,
%                                                or outside the nodes of the 'Climatology' grids
%   referenceSystem:*                            propagated unwrapped
%
% REFERENCES
%   Knocke P.C., Ries J.C., Tapley B.D. (1988) Earth radiation pressure effects on
%   satellites, AIAA/AAS Astrodynamics Conf., AIAA-88-4292, Table 1, eqs 3-4, sect. III.
%   NASA LaRC, CERES_EBAF_Ed4.2 Data Quality Summary (2025) (20 km reference level,
%   geodetic weights, albedos above 1). Rodriguez-Solano C.J. et al. (2012) Adv. Space
%   Res. 49:1113 (CERES-based albedo and emission grids for Earth radiation pressure).
%
% EXAMPLE
%   ep = referenceSystem.mjdToEpoch(58849.5, "UTC"); r = [6.88e6; 0; 0];
%   t = NonGravForces.env.earthRadiance("knocke", ep, [12 24], r);
%   sum(t.area_m2) / (2*pi*t.radius_m^2*(1 - cos(t.capAngle_rad)))     % 1
%
% See also NonGravForces.env.radiationSources, NonGravForces.Options.
%
% KEYWORDS albedo emissivity Knocke CERES Earth radiation tiles cap
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic); ceres 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    backend (1,1) string {mustBeMember(backend, ["pointSource","knocke","ceres","ceresClimatology"])}
    epoch (1,1) struct
    grid double
    r double
    nv.Ceres = []
    nv.RotationItrsToGcrs = []
    nv.Climatology = []
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
gridded = backend == "ceres" || backend == "ceresClimatology";
if ~gridded && ~isempty(nv.RotationItrsToGcrs)
    NonGravForces.fail(unit, "inputCeres", "'RotationItrsToGcrs' belongs to the gridded backends -> remove it, or choose ceres or ceresClimatology");
end
if backend ~= "ceres" && ~isempty(nv.Ceres)
    NonGravForces.fail(unit, "inputCeres", "'Ceres' belongs to backend ceres -> remove it, or choose ceres");
end
if backend ~= "ceresClimatology" && ~isempty(nv.Climatology)
    NonGravForces.fail(unit, "inputClimatology", "'Climatology' belongs to backend ceresClimatology -> remove it, or choose ceresClimatology");
end
co = localCoefficients();
utc = referenceSystem.convertTimeScale(epoch, "UTC");
t0 = referenceSystem.datetimeToEpoch(datetime(1981, 12, 22), "UTC");
days = (double(utc.day(:)') - double(t0.day)) + double(utc.sec(:)') / 86400;
seasonal = cos(2 * pi / co.period_days * days);            % 1xK
ell = referenceSystem.ellipsoidParameters("WGS84");
R = ell.a;
if backend == "ceres"
    cer = localCheckCeres(unit, nv.Ceres, nv.RotationItrsToGcrs, K);
    R = R + cer.referenceHeight_m;                          % the product's flux reference level
elseif backend == "ceresClimatology"
    clm = localCheckClimatology(unit, nv.Climatology, nv.RotationItrsToGcrs, K, ell);
    R = R + clm.referenceHeight_m;                          % MEASURED from the file, not declared
end
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
        tiles.grid = [0 0]; tiles.M = 0; tiles.exitanceLw_W_m2 = [];
        tiles.n = zeros(3, 0, K); tiles.area_m2 = zeros(0, K); tiles.lat_rad = zeros(0, K);
        tiles.albedo = repmat(co.a0, 1, K);
        tiles.emissivity = repmat(co.e0, 1, K);
    otherwise
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
        tiles.exitanceLw_W_m2 = [];
        if backend == "knocke"
            P2 = 0.5 * (3 * x.^2 - 1);
            tiles.albedo     = co.a0 + co.a1 * (x .* seasonal) + co.a2 * P2;     % MxK
            tiles.emissivity = co.e0 + co.e1 * (x .* seasonal) + co.e2 * P2;
        elseif backend == "ceres"
            tiles = localCeres(unit, tiles, cer, nv.RotationItrsToGcrs, utc, ell, M, K);
        else
            tiles = localClimatology(unit, tiles, clm, nv.RotationItrsToGcrs, utc, M, K);
        end
end
tiles.provenance = NonGravForces.provenance("env.earthRadiance", 'Options', struct('backend', backend, 'grid', grid));
end

function cer = localCheckCeres(unit, cer, R, K)
%LOCALCHECKCERES The fields earthRadiance reads from an io.importCERES struct, and the rotation.
need = ["mjd" "lat_deg" "lon_deg" "rsut" "rlut" "rsdt" "referenceHeight_m" "latitude"];
if ~(isstruct(cer) && isscalar(cer) && all(isfield(cer, need)))
    NonGravForces.fail(unit, "inputCeres", "backend ceres needs 'Ceres', the struct io.importCERES returns (fields %s) -> import the EBAF files first", strjoin(need, " "));
end
T = numel(cer.mjd); sz = [numel(cer.lat_deg) numel(cer.lon_deg) T];
if T < 2 || ~isequal(size(cer.rsut, 1:3), sz) || ~isequal(size(cer.rlut, 1:3), sz) || ~isequal(size(cer.rsdt, 1:3), sz) || cer.latitude ~= "geodetic"
    NonGravForces.fail(unit, "inputCeres", "'Ceres' grids must be lat x lon x month with at least two months and geodetic latitudes -> use io.importCERES with a 'Period' of at least one day");
end
if ~(isnumeric(R) && isequal(size(R, 1:3), [3 3 K]) && all(isfinite(R(:))))
    NonGravForces.fail(unit, "inputCeres", "backend ceres needs 'RotationItrsToGcrs' 3x3x%d -> pass referenceSystem.rotationItrsToGcrs(epoch, 'Eop', eop)", K);
end
end

function tiles = localCeres(unit, tiles, cer, R, utc, ell, M, K)
%LOCALCERES Interpolate the CERES grids at the tile centres: geodetic bilinear in space, linear in time.
nI = pagemtimes(R, 'transpose', tiles.n, 'none');                  % 3xMxK tile directions in the ITRS
x = reshape(nI(1, :, :), M, K); y = reshape(nI(2, :, :), M, K); z = reshape(nI(3, :, :), M, K);
f = ell.f;
latD = atan2d(z, (1 - f)^2 * hypot(x, y));                          % geodetic latitude of the grid [deg]
lon = mod(atan2d(y, x), 360);
tiles.lat_rad = deg2rad(latD);
% time: the two nodes bracketing each epoch
t = double(utc.day(:)') + double(utc.sec(:)') / 86400;               % MJD UTC, 1xK
k0 = discretize(t, cer.mjd);
k0(t == cer.mjd(end)) = numel(cer.mjd) - 1;
if any(isnan(k0))
    bad = find(isnan(k0), 1);
    NonGravForces.fail(unit, "dataRange", "epoch %d (MJD %.4f) lies outside the CERES months held (MJD %.2f .. %.2f) -> import with a 'Period' that covers every epoch", ...
        bad, t(bad), cer.mjd(1), cer.mjd(end));
end
wt = (t - cer.mjd(k0)) ./ (cer.mjd(k0 + 1) - cer.mjd(k0));             % 1xK
% space: cell-centre indices and weights (uniform grids, checked by io.importCERES)
la = cer.lat_deg; lo = cer.lon_deg; nLat = numel(la); nLon = numel(lo);
u = (min(max(latD, la(1)), la(end)) - la(1)) / (la(2) - la(1));      % 0 .. nLat-1
i0 = min(floor(u), nLat - 2); wy = u - i0; i0 = i0 + 1;               % rows i0, i0+1
v = mod(lon - lo(1), 360) / (lo(2) - lo(1));                          % 0 .. nLon
j0 = floor(v); wx = v - j0; j0 = mod(j0, nLon) + 1; j1 = mod(j0, nLon) + 1;
F = struct();
for name = ["rsut" "rlut" "rsdt"]
    G = cer.(name);
    val = zeros(M, K);
    for side = 0:1
        kk = k0 + side;                                                % 1xK month index
        page = (kk - 1) * nLat * nLon;                                 % linear offset of that month
        g00 = double(G(i0 + (j0 - 1) * nLat + page));
        g10 = double(G(i0 + 1 + (j0 - 1) * nLat + page));
        g01 = double(G(i0 + (j1 - 1) * nLat + page));
        g11 = double(G(i0 + 1 + (j1 - 1) * nLat + page));
        b = g00 + wy .* (g10 - g00) + wx .* (g01 - g00) + (wx .* wy) .* ((g11 - g10) - (g01 - g00));
        if side == 0, val = (1 - wt) .* b; else, val = val + wt .* b; end
    end
    F.(name) = val;
end
alb = zeros(M, K);
lit = F.rsdt > 0;
alb(lit) = F.rsut(lit) ./ F.rsdt(lit);
tiles.albedo = alb;
tiles.emissivity = [];
tiles.exitanceLw_W_m2 = F.rlut;
tiles.ceres = struct('product', string(cer.product), 'version', string(cer.version), ...
    'months', datetime(cer.mjd([k0; k0 + 1]) + 2400000.5, 'ConvertFrom', 'juliandate'), 'timeWeight', wt);
end

function clm = localCheckClimatology(unit, clm, R, K, ell)
%LOCALCHECKCLIMATOLOGY The two GROOPS grids the ceresClimatology backend reads, and the rotation.
%   The two files are separate products and nothing in either says it belongs with the other, so
%   every shared quantity is COMPARED here rather than taken from whichever was read first.
need = ["albedo" "longwave"];
if ~(isstruct(clm) && isscalar(clm) && all(isfield(clm, need)))
    NonGravForces.fail(unit, "inputClimatology", "backend ceresClimatology needs 'Climatology', a struct with fields albedo and longwave, each an io.importGroopsGrid result -> import CERES_SYN1deg_albedo_climatology.dat and CERES_SYN1deg_longwaveFlux_climatology.dat");
end
inner = ["lon_deg" "lat_deg" "radius_m" "values" "mjd" "isTimeSeries"];
for f = need
    g = clm.(f);
    if ~(isstruct(g) && isscalar(g) && all(isfield(g, inner)))
        NonGravForces.fail(unit, "inputClimatology", "Climatology.%s must be an io.importGroopsGrid result (fields %s) -> read the file with io.importGroopsGrid", f, strjoin(inner, " "));
    end
    if ~g.isTimeSeries || numel(g.mjd) ~= size(g.values, 3) + 1
        NonGravForces.fail(unit, "inputClimatology", "Climatology.%s must be a griddedDataTimeSeries whose node count is one more than its grid count (got %d nodes, %d grids) -> the monthly climatology is piecewise constant between nodes", f, numel(g.mjd), size(g.values, 3));
    end
end
a = clm.albedo; b = clm.longwave;
if ~isequal(a.lon_deg, b.lon_deg) || ~isequal(a.lat_deg, b.lat_deg) || ~isequal(a.mjd, b.mjd) || ~isequal(size(a.values), size(b.values))
    NonGravForces.fail(unit, "inputClimatology", "the two Climatology grids must share their axes and time nodes -> they are separate files and nothing pairs them but this check");
end
% Longitude is indexed arithmetically below, so its uniformity is asserted rather than assumed --
% but to the tolerance the FILE actually holds, not to float noise. The SYN1deg axis jitters by
% 1.366e-05 degrees between neighbouring centres and starts at -179.499994 rather than -179.5
% (measured 2026-09-17); that is 3.2e-08 relative, which is single-precision, so the axis was
% produced in float32 somewhere upstream. It is uniform in intent and not in its bytes.
% Consequences, both handled: the spacing used below is the AVERAGE over the whole axis rather
% than the first gap, so the jitter cannot accumulate across 360 cells, and the residual
% mispositioning is under 1e-05 degrees of a 1 degree cell. Latitude is NOT indexed this way: the
% real centres are 89.4966, 88.4899, ... with spacings from 0.9933 to 1.0067, and are searched.
d = diff(a.lon_deg);
if numel(a.lon_deg) < 2 || numel(a.lat_deg) < 2 || d(1) <= 0 || max(abs(d - mean(d))) > 1e-3 * mean(d)
    NonGravForces.fail(unit, "inputClimatology", "the Climatology longitudes must be at least two, ascending and uniform to 1e-3 of their spacing (this file: %.3e) -> this backend indexes longitude arithmetically and searches latitude", max(abs(d - mean(d))) / mean(d));
end
if any(diff(a.lat_deg) >= 0)
    NonGravForces.fail(unit, "inputClimatology", "the Climatology latitudes must be strictly decreasing, as io.importGroopsGrid returns them -> do not reorder the grid");
end
if ~(isnumeric(R) && isequal(size(R, 1:3), [3 3 K]) && all(isfinite(R(:))))
    NonGravForces.fail(unit, "inputClimatology", "backend ceresClimatology needs 'RotationItrsToGcrs' 3x3x%d -> pass referenceSystem.rotationItrsToGcrs(epoch, 'Eop', eop)", K);
end
% THE REFERENCE HEIGHT IS MEASURED, AND MEASURING IT IS WHAT CAUGHT THE LATITUDE CONVENTION (R3).
% The file's per-row radii equal the ellipsoid radius at the stored latitude READ AS GEOCENTRIC
% plus a constant -- 30000.1 m for the SYN1deg climatology, to a spread of 0.168 m over all 180
% rows. Read as geodetic the same comparison spreads by 71.6 m, 400 times worse, and that is how
% this backend learned that io.importGroopsGrid's latitudes are geocentric rather than geodetic as
% its help had said. A consumer that converts to geodetic before the lookup misplaces every tile by
% up to 0.19 deg, a fifth of a cell, in a way no smoke test would show.
th = deg2rad(a.lat_deg(:));
rEll = ell.a * ell.b ./ sqrt((ell.a * sin(th)).^2 + (ell.b * cos(th)).^2);   % radius at GEOCENTRIC latitude
h = a.radius_m(:) - rEll;
if any(~isfinite(h)) || max(h) - min(h) > 1
    NonGravForces.fail(unit, "inputClimatology", "the Climatology radii are not one constant height above the WGS84 ellipsoid at the stored (geocentric) latitudes: spread %.3f m -> either the grid is not this product, or its latitudes are geodetic after all", max(h) - min(h));
end
clm.referenceHeight_m = mean(h);
end

function tiles = localClimatology(unit, tiles, clm, R, utc, M, K)
%LOCALCLIMATOLOGY Interpolate the GROOPS monthly climatology at the tile centres.
%   It takes no ellipsoid: the grid's axis is GEOCENTRIC, so the tile latitude needs no conversion.
%   Bilinear in space exactly as the ceres backend, and PIECEWISE CONSTANT in time: the container
%   bounds T grids by T+1 nodes and each grid IS the month's mean, so interpolating between two of
%   them would invent a value the product does not carry.
nI = pagemtimes(R, 'transpose', tiles.n, 'none');                  % 3xMxK tile directions in the ITRS
x = reshape(nI(1, :, :), M, K); y = reshape(nI(2, :, :), M, K); z = reshape(nI(3, :, :), M, K);
% GEOCENTRIC latitude, because that is what the file's axis holds (see localCheckClimatology).
% The ceres backend above converts to geodetic because ITS product is on a geodetic axis; the two
% differ by up to 0.19 deg and the conversion belongs to the grid, not to the tile.
latD = atan2d(z, hypot(x, y));
lon = mod(atan2d(y, x), 360);
tiles.lat_rad = deg2rad(latD);
g = clm.albedo;
la = g.lat_deg(:)'; lo = g.lon_deg(:)'; nLat = numel(la); nLon = numel(lo);
% time: the ONE grid whose node interval contains the epoch
t = double(utc.day(:)') + double(utc.sec(:)') / 86400;               % MJD UTC, 1xK
kk = discretize(t, g.mjd);
kk(t == g.mjd(end)) = numel(g.mjd) - 1;
if any(isnan(kk))
    bad = find(isnan(kk), 1);
    NonGravForces.fail(unit, "dataRange", "epoch %d (MJD %.4f) lies outside the climatology's nodes (MJD %.2f .. %.2f) -> the file covers %s to %s", ...
        bad, t(bad), g.mjd(1), g.mjd(end), string(datetime(g.mjd(1) + 2400000.5, 'ConvertFrom', 'juliandate'), 'yyyy-MM'), ...
        string(datetime(g.mjd(end) + 2400000.5, 'ConvertFrom', 'juliandate'), 'yyyy-MM'));
end
% space: latitude is SEARCHED because the file's centres are not uniformly spaced; longitude is
% uniform and periodic, and was checked to be so.
lam = min(max(latD, la(end)), la(1));                                % clamp poleward of the centres
i0 = discretize(-lam, -la);                                          % -la is increasing
i0(lam == la(1)) = 1;
i0 = min(max(i0, 1), nLat - 1);
% RESHAPE, because la is a ROW and i0 is MxK: for a vector source MATLAB returns the source's
% orientation, so la(i0) comes back 1x(M*K) and subtracting an MxK latitude would broadcast to
% (M*K)x(M*K). With K = 1 that is a square array of the right size in one dimension and the wrong
% one in the other, which is exactly the shape a test compares against and passes for K > 1.
la0 = reshape(la(i0), size(i0));
la1 = reshape(la(i0 + 1), size(i0));
wy = (la0 - lam) ./ (la0 - la1);                                    % 0 at row i0, 1 at row i0+1
dLon = (lo(end) - lo(1)) / (nLon - 1);            % the AVERAGE spacing; see localCheckClimatology
v = mod(lon - lo(1), 360) / dLon;
j0 = floor(v); wx = v - j0; j0 = mod(j0, nLon) + 1; j1 = mod(j0, nLon) + 1;
page = (kk - 1) * nLat * nLon;                                       % linear offset of that month
out = struct();
for name = ["albedo" "longwave"]
    G = clm.(name).values;
    g00 = double(G(i0 + (j0 - 1) * nLat + page));
    g10 = double(G(i0 + 1 + (j0 - 1) * nLat + page));
    g01 = double(G(i0 + (j1 - 1) * nLat + page));
    g11 = double(G(i0 + 1 + (j1 - 1) * nLat + page));
    out.(name) = g00 + wy .* (g10 - g00) + wx .* (g01 - g00) + (wx .* wy) .* ((g11 - g10) - (g01 - g00));
end
tiles.albedo = out.albedo;
tiles.emissivity = [];
tiles.exitanceLw_W_m2 = out.longwave;                                % a MEASURED flux, W m^-2
tiles.climatology = struct('referenceHeight_m', clm.referenceHeight_m, 'gridIndex', kk, ...
    'month', datetime(g.mjd(kk) + 2400000.5, 'ConvertFrom', 'juliandate'), ...
    'nodes', [g.mjd(1) g.mjd(end)], 'gridCount', size(g.values, 3));
end

function co = localCoefficients()
% Knocke, Ries & Tapley (1988) Table 1. The test compares these against
% tools/pyxcheck/knocke_coefficients.csv, an independent restatement.
co = struct('a0', 0.34, 'a1', 0.10, 'a2', 0.29, 'e0', 0.68, 'e1', -0.07, 'e2', -0.18, ...
            'period_days', 365.25, 't0', "1981-12-22");
end
