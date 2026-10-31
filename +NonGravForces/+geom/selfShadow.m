function [f, info] = selfShadow(macro, dir_body, options)
%SELFSHADOW Unobstructed fraction of every plate toward body-frame directions: ray cast or cube-map table.
%
%   [f, info] = NonGravForces.geom.selfShadow(macro, dir_body)                       % ray cast
%   [f, info] = NonGravForces.geom.selfShadow(macro, dir_body, 'Method', "lut", 'LutResolution', 16)
%   lut = NonGravForces.geom.selfShadow(macro, [], 'Method', "lut")                   % table only
%   f = NonGravForces.geom.selfShadow(macro, dir_body, 'Method', "lut", 'Lut', lut)   % reuse a table
%
%   OCCLUSION, AND NOTHING ELSE. For each plate k and each direction d (unit
%   vector FROM the body TOWARD a radiation source, or toward where gas molecules
%   come from), f(k) is the fraction of the plate's area from which a ray along
%   +d leaves the body without striking another plate. Whether the plate FACES d
%   is not this unit's fact: every consumer multiplies by (n_k . d)+ itself
%   (kernel.radiationPanelForce, geom.thermalState, geom.flowShadow). Until
%   2026-09-15 this unit also returned 0 for back-facing plates, and that
%   radiation-beam convention, applied to a thermal gas, deleted the edge-on
%   shear of every plate the stream grazes (HANDOVER LP-C, finding F1).
%
%   For a back-facing direction (n_k . d < 0) the value returned is that of the
%   MIRROR direction d - 2 (n_k . d) n_k, which lies in the plate's front
%   hemisphere. It is a continuous extension across the plate's horizon: no
%   consumer weights it, and it keeps the table's interpolation honest at grazing
%   incidence, where a raw back-facing ray would run into the body and drag the
%   interpolated value toward 0. A direction lying in the plate's own plane is
%   replaced by its front-side limit, 1e-6 rad out of the plane, so both sides of
%   the horizon approach one value. A convex closed body has f = 1 on every plate;
%   GRACE-FO is not convex (the array flanks overhang the nadir plate).
%
%   METHODS
%     raycast   Stratified sampling: an m x m grid of points on each plate
%               polygon (m = 'Samples', default 40; the polygon is mapped
%               bilinearly from the unit square for quadrilaterals, and by fan
%               triangles otherwise); a point is occluded when the ray from it
%               along +d intersects any other plate polygon in front of it.
%               Exact to the sampling: the error of the area fraction is at
%               most ~2/m for a straight shadow edge.
%     lut       Equiangular cube map. The ray cast is evaluated once at the
%               (n+1)^2 nodes of each of the six faces of a cube projected onto
%               the unit sphere, uniform in the two face angles (n =
%               'LutResolution' cells per face edge; 6 n^2 + 2 distinct
%               directions, nodes shared by two or three faces cast once), and
%               read back by bilinear interpolation in those angles. The look-up
%               is O(1) per direction: the face is the direction's largest
%               component and the cell follows from two arctangents. The
%               interpolant is continuous across face edges. The table is
%               returned in info.lut and can be passed back as 'Lut'.
%
% INPUT ARGUMENTS
%   Name       Size  Type                          Units  Description
%   ----       ----  ----                          -----  -----------
%   macro      1x1   NonGravForces.geom.MacroModel  -      With vertices on every plate (hasVertices)
%   dir_body   3xK   double                         -      Unit directions from the body, body frame;
%                                                          [] with Method "lut" returns only the table
%
% NAME-VALUE ARGUMENTS
%   Name            Size  Type    Units  Description
%   ----            ----  ----    -----  -----------
%   Method          1x1   string  -      "raycast" (default) | "lut"
%   Samples         1x1   double  -      grid points per plate edge for the ray cast (default 40)
%   LutResolution   1x1   double  -      cells per cube-face edge, integer in [1, 256]; default:
%                                         the resolution of 'Lut' if given, else 16 (5.6 deg)
%   Lut             1x1   struct  -      a previously returned table for this macro model; its own
%                                         sample count is kept and reported in info.samples (Samples
%                                         is not compared with it)
%   Plate           1x1   double  -      evaluate this one plate only (index into the macro model):
%                                         f is then 1xK; default [] (every plate)
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   f        NxK   double  -      Unobstructed fraction per plate and direction, in [0, 1]; or,
%                                 with dir_body = [] and Method "lut", the table struct itself
%   info     1x1   struct  -      .method, .samples, .lut (struct .f N x 6(n+1)^2, .nodes 3 x 6(n+1)^2
%                                 node directions face by face, .resolution, .samples, .macroName,
%                                 .geometry; only for "lut"), .provenance (built only when info
%                                 is requested)
%
% ACCURACY
%   Ray cast: area-fraction error <= ~2/Samples at a shadow edge (a straight
%   edge crossing an m x m grid misclassifies at most ~2m of m^2 cells): 40
%   samples give ~5e-2 worst case and ~1e-2 typical, measured against the
%   pure-Python ray cast tools/pyxcheck/selfshadow_check.py (tests/tNgfGeometry.m).
%   Table: bilinear interpolation adds the angular discretisation, 90/n deg
%   between nodes (5.6 deg at the default 16); a shadow edge sharper than that is
%   smeared over one cell, and the mean table-vs-ray-cast difference over random
%   directions is recorded per run with its convergence over n. The plate polygons
%   themselves approximate the exterior geometry (geom.graceFO: ~1-4 cm).
%
% ERRORS
%   NonGravForces:selfShadow:inputMacro       not a MacroModel, or without vertices on every plate
%   NonGravForces:selfShadow:inputShape       dir_body not 3xK, not unit, or not finite; 'Plate' not a plate index
%   NonGravForces:selfShadow:inputResolution  LutResolution not an integer in [1, 256]
%   NonGravForces:selfShadow:inputLut         'Lut' built for another geometry or resolution (or not a table)
%
% REFERENCES
%   Wöske et al. (2019) Adv. Space Res. 63:1096 (finite-element ray tracing vs plates);
%   Ronchi C., Iacono R., Paolucci P.S. (1996) J. Comput. Phys. 124:93 (the equiangular
%   cubed sphere); HANDOVER LP-C WP-C1 and WP-C3.
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C");
%   f = NonGravForces.geom.selfShadow(mm, [0; sind(70); cosd(70)])   % low elevation: the port flank shadows the nadir plate
%   lut = NonGravForces.geom.selfShadow(mm, [], 'Method', "lut");
%   fl = NonGravForces.geom.selfShadow(mm, [0; sind(70); cosd(70)], 'Method', "lut", 'Lut', lut);
%
% See also NonGravForces.geom.flowShadow, NonGravForces.geom.MacroModel, NonGravForces.env.radiationSources.
%
% KEYWORDS self-shadowing occlusion ray cast look-up table cube map
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic); occlusion-only semantics and cube map 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    macro
    dir_body double
    options.Method (1,1) string {mustBeMember(options.Method, ["raycast","lut"])} = "raycast"
    options.Samples (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.Samples, 2)} = 40
    options.LutResolution double = []
    options.Lut struct = struct([])
    options.Plate double = []
end
unit = "selfShadow";
if ~isa(macro, 'NonGravForces.geom.MacroModel') || ~macro.hasVertices
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel with vertices_m on every plate -> geom.graceFO and geom.boxWing (bus) provide them");
end
plates = 1:macro.n;
if ~isempty(options.Plate)
    if ~(isscalar(options.Plate) && any(options.Plate == plates))
        NonGravForces.fail(unit, "inputShape", "'Plate' must be one plate index in 1..%d -> got %s", macro.n, mat2str(options.Plate));
    end
    plates = options.Plate;
end
if ~isempty(dir_body)
    if ~(ismatrix(dir_body) && size(dir_body, 1) == 3) || any(~isfinite(dir_body(:)))
        NonGravForces.fail(unit, "inputShape", "dir_body must be a finite 3xK array -> one unit direction per column");
    end
    if any(abs(vecnorm(dir_body) - 1) > 1e-6)
        NonGravForces.fail(unit, "inputShape", "dir_body columns must be unit vectors (norm off by %.2e) -> normalise", max(abs(vecnorm(dir_body) - 1)));
    end
elseif options.Method == "raycast"
    NonGravForces.fail(unit, "inputShape", "dir_body is empty -> pass directions, or Method lut to build the table only");
end
switch options.Method
    case "raycast"
        f = localRaycast(macro, dir_body, options.Samples, plates);
        lut = [];
    case "lut"
        n = options.LutResolution;
        if isempty(n)
            if ~isempty(options.Lut) && isfield(options.Lut, 'resolution'), n = options.Lut.resolution; else, n = 16; end
        end
        if ~(isscalar(n) && isfinite(n) && n == round(n) && n >= 1 && n <= 256)
            NonGravForces.fail(unit, "inputResolution", "LutResolution must be an integer number of cells per cube-face edge in [1, 256] -> got %s; 16 gives 5.6 deg", mat2str(n));
        end
        if ~isempty(options.Lut)
            lut = options.Lut;
            ok = isfield(lut, 'geometry') && isfield(lut, 'resolution') && isfield(lut, 'f') && ...
                 lut.resolution == n && isequal(lut.geometry, localGeometry(macro)) && ...
                 isequal(size(lut.f), [macro.n, 6 * (n + 1)^2]);
            if ~ok
                NonGravForces.fail(unit, "inputLut", "the Lut was not built for this macro model at resolution %g -> build one with selfShadow(macro, [], 'Method', ""lut"")", n);
            end
        else
            lut = localBuildTable(macro, n, options.Samples);
        end
        if isempty(dir_body)
            f = lut;                                             % table only
        else
            f = localLookup(lut, dir_body, plates);
        end
end
if nargout > 1
    % the sample count the numbers were made with: a reused table's own, not the argument (2026-10-31; the
    % argument was reported even when a table built with another count was looked up)
    samples = options.Samples;
    if ~isempty(lut) && isfield(lut, 'samples'), samples = lut.samples; end
    info.method = options.Method;
    info.samples = samples;
    info.lut = lut;
    info.provenance = NonGravForces.provenance("geom.selfShadow", 'Options', struct('method', options.Method, 'samples', samples, 'lutResolution', options.LutResolution));
end
end

% ------------------------------------------------------------------------- ray cast
function f = localRaycast(macro, d, m, plates)
% d: 3xK unit directions. Sample points on each plate, occlusion by every other plate. A
% back-facing direction is replaced, for that plate only, by its mirror in the plate's plane.
N = macro.n; K = size(d, 2);
f = zeros(numel(plates), K);
[P, w] = localSamples(macro, m);                     % P{k} 3xS points, w{k} 1xS weights (sum 1)
for row = 1:numel(plates)
    k = plates(row);
    nk = macro.N(:, k);
    g = nk' * d;                                     % 1xK
    dk = d;
    back = g < 0;
    if any(back)
        dk(:, back) = d(:, back) - 2 * nk .* g(back);    % mirror into the front hemisphere
    end
    % A direction in the plate's own plane is replaced by its front-side limit, 1e-6 rad out of the
    % plane: a ray skimming the plane otherwise meets a neighbouring plate exactly on their shared edge
    % and counts as occluded, a value that belongs to neither side of the horizon.
    gm = abs(g);
    skim = gm < 1e-6;
    if any(skim)
        dk(:, skim) = dk(:, skim) + (1e-6 - gm(skim)) .* nk;
        dk(:, skim) = dk(:, skim) ./ vecnorm(dk(:, skim));
    end
    S = size(P{k}, 2);
    occluded = false(S, K);
    for j = 1:N
        if j == k, continue; end
        nj = macro.N(:, j); Vj = macro.V{j};
        denom = nj' * dk;                                                       % 1xK
        hit = abs(denom) > 1e-12;
        if ~any(hit), continue; end
        % t = n_j . (p_j - P) / (n_j . d), one value per (sample, direction)
        t = (nj' * (Vj(:, 1) - P{k})).' ./ denom;                              % SxK
        ahead = t > 1e-9;
        % intersection points and the point-in-polygon test in the plane of j
        e1 = Vj(:, 2) - Vj(:, 1); e1 = e1 / norm(e1); e2 = cross(nj, e1);
        px = e1' * (Vj - Vj(:, 1)); py = e2' * (Vj - Vj(:, 1));                 % polygon in 2-D
        Xk = reshape(P{k}, 3, S, 1) + reshape(t, 1, S, K) .* reshape(dk, 3, 1, K);    % 3xSxK hit points
        Xr = Xk - Vj(:, 1);
        hx = reshape(sum(Xr .* e1, 1), S, K); hy = reshape(sum(Xr .* e2, 1), S, K);
        inside = localInPolygon(hx, hy, px, py);
        occluded = occluded | (ahead & inside & hit);
    end
    % column-wise sum (not a BLAS matrix product): the result must not depend on how many directions are batched;
    % fully occluded and fully exposed columns are set exactly (the weights sum to 1 only to rounding)
    frac = sum(w{k}' .* double(occluded), 1);
    frac(all(occluded, 1)) = 1;
    frac(~any(occluded, 1)) = 0;
    f(row, :) = 1 - frac;
end
f = min(1, max(0, f));
end

function [P, w] = localSamples(macro, m)
% Stratified sample points on each plate polygon: bilinear map of an m x m grid for quadrilaterals,
% fan triangles (each with m x m / #triangles points) otherwise. Weights are the area shares.
P = cell(1, macro.n); w = cell(1, macro.n);
u = ((1:m) - 0.5) / m; [U, W] = ndgrid(u, u); U = U(:)'; W = W(:)';
for k = 1:macro.n
    V = macro.V{k}; nV = size(V, 2);
    if nV == 4
        pts = (1 - U) .* (1 - W) .* V(:, 1) + U .* (1 - W) .* V(:, 2) + U .* W .* V(:, 3) + (1 - U) .* W .* V(:, 4);
        % bilinear map of a planar quadrilateral: Jacobian varies for non-parallelograms; weight by it
        dU = -(1 - W) .* V(:, 1) + (1 - W) .* V(:, 2) + W .* V(:, 3) - W .* V(:, 4);
        dW = -(1 - U) .* V(:, 1) - U .* V(:, 2) + U .* V(:, 3) + (1 - U) .* V(:, 4);
        jac = vecnorm(cross(dU, dW, 1));
        P{k} = pts; w{k} = jac / sum(jac);
    else
        c = mean(V, 2); pts = []; wts = [];
        for i = 1:nV
            a = V(:, i); b = V(:, mod(i, nV) + 1);
            % uniform points in the triangle (c, a, b) by the square-root map
            r1 = sqrt(U); pt = c + r1 .* ((1 - W) .* (a - c) + W .* (b - c));
            area = 0.5 * norm(cross(a - c, b - c));
            pts = [pts, pt]; wts = [wts, area * ones(1, numel(U)) / numel(U)]; %#ok<AGROW>
        end
        P{k} = pts; w{k} = wts / sum(wts);
    end
end
end

function inside = localInPolygon(hx, hy, px, py)
% Convex polygon (px, py) counter-clockwise: a point is inside when it lies left of every edge.
inside = true(size(hx));
nV = numel(px);
for i = 1:nV
    j = mod(i, nV) + 1;
    cr = (px(j) - px(i)) .* (hy - py(i)) - (py(j) - py(i)) .* (hx - px(i));
    inside = inside & (cr >= -1e-12);
end
end

% ------------------------------------------------------------------------- cube map
function g = localGeometry(macro)
% Everything the ray cast reads, as one column: a table built for other geometry is refused.
g = [macro.N(:); cell2mat(cellfun(@(v) v(:), macro.V(:), 'UniformOutput', false))];
end

function lut = localBuildTable(macro, n, samples)
% Ray cast at the distinct cube-map nodes, then expand to the per-face node layout.
tv = tan(-pi/4 + (0:n) * (pi / 2) / n);
tv(1) = -1; tv(end) = 1;                                  % the face edges exactly, so shared nodes coincide
tv = (tv - flip(tv)) / 2;                                 % exactly antisymmetric about the face centre
UA = [2 3 1]; VA = [3 1 2];
[I, J] = ndgrid(1:n+1, 1:n+1);                            % node (i, j): i along the u axis, j along v
raw = zeros(3, 6 * (n + 1)^2);
for F = 1:6
    m = ceil(F / 2); sg = 1 - 2 * mod(F + 1, 2);          % F = 1 +x, 2 -x, 3 +y, 4 -y, 5 +z, 6 -z
    cols = (F - 1) * (n + 1)^2 + (1:(n + 1)^2);
    raw(m, cols) = sg;
    raw(UA(m), cols) = tv(I(:));
    raw(VA(m), cols) = tv(J(:));
end
[uniq, ~, index] = unique(raw', 'rows');                  % exact: shared edge and corner nodes are identical
nodes = uniq' ./ vecnorm(uniq');
fu = localRaycast(macro, nodes, samples, 1:macro.n);      % N x U
lut = struct('f', fu(:, index), 'nodes', nodes(:, index), 'resolution', n, 'samples', samples, ...
             'macroName', macro.name, 'geometry', localGeometry(macro));
end

function f = localLookup(lut, d, plates)
% O(1) bilinear read of the cube map, element-wise, so a batch equals single calls bitwise.
n = lut.resolution; K = size(d, 2);
[~, m] = max(abs(d), [], 1);                              % major axis, ties to the lower index
col = (0:K-1) * 3;
dm = d(col + m);
face = 2 * m - (dm > 0);                                  % +x 1, -x 2, +y 3, -y 4, +z 5, -z 6
UA = [2 3 1]; VA = [3 1 2];
adm = abs(dm);
a = atan(d(col + UA(m)) ./ adm);                          % face angles in [-pi/4, pi/4]
b = atan(d(col + VA(m)) ./ adm);
x = (a + pi/4) / (pi/2) * n; y = (b + pi/4) / (pi/2) * n;
i0 = min(max(floor(x), 0), n - 1); j0 = min(max(floor(y), 0), n - 1);
wx = min(max(x - i0, 0), 1); wy = min(max(y - j0, 0), 1);
c00 = (face - 1) * (n + 1)^2 + j0 * (n + 1) + i0 + 1;
T = lut.f(plates, :);
% Bilinear interpolation in difference form: four equal corners return that value EXACTLY (the
% weighted sum of four equal values lands one ulp off), so an unobstructed region reads 1 bitwise,
% which geom.flowShadow's limit "nothing in front of the plate gives exactly 1" relies on. The three
% differences are formed once per table node (cell origin c: corners c, c+1, c+n+1, c+n+2), not once
% per direction: the same subtractions of the same values, so the result is bitwise the same, read
% with four gathers instead of five (2026-10-16, interleaved: 109 ns per direction before, 88 after).
L = size(T, 2); s1 = 1:L-n-2;
Dx = zeros(size(T)); Dy = Dx; Dxy = Dx;
Dx(:, s1) = T(:, s1 + 1) - T(:, s1);
Dy(:, s1) = T(:, s1 + n + 1) - T(:, s1);
Dxy(:, s1) = (T(:, s1 + n + 2) - T(:, s1 + 1)) - (T(:, s1 + n + 1) - T(:, s1));
f = T(:, c00) + wx .* Dx(:, c00) + wy .* Dy(:, c00) + (wx .* wy) .* Dxy(:, c00);
% A convex combination of values in [0, 1] can land one ulp outside; the clamp restores what the
% algebra guarantees. A non-finite value cannot arise: every factor is finite by construction.
f = min(1, max(0, f));
end
