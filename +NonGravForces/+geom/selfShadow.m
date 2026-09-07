function [f, info] = selfShadow(macro, dir_body, options)
%SELFSHADOW Exposed (unoccluded) fraction of every plate for light arriving from given body-frame directions.
%
%   [f, info] = NonGravForces.geom.selfShadow(macro, dir_body)                       % ray cast
%   [f, info] = NonGravForces.geom.selfShadow(macro, dir_body, 'Method', "lut", 'LutResolution', 642)
%   lut = NonGravForces.geom.selfShadow(macro, [], 'Method', "lut", 'LutResolution', 642)   % table only
%
%   For each plate k and each direction d (unit vector FROM the body TO the
%   source, the convention of env.radiationSources' src.dir), f(k) is the
%   fraction of the plate's area that a parallel beam from d reaches without
%   hitting another plate. Plates facing away from the beam (n_k . d <= 0)
%   get f = 0 (they are not illuminated; the kernel applies the cosine
%   separately). A convex closed body has f = 1 on every lit plate; GRACE-FO
%   is not convex (the array flanks overhang the nadir plate), so f < 1
%   occurs for low elevations, which is why this exists.
%
%   METHODS
%     raycast   Stratified sampling: an m x m grid of points on each plate
%               polygon (m = 'Samples', default 40; the polygon is mapped
%               bilinearly from the unit square for quadrilaterals, and by
%               fan triangles otherwise); a point is occluded when the ray
%               from it along +d intersects any other plate polygon in
%               front of it. Exact to the sampling: the error of the area
%               fraction is at most ~2/m for a straight shadow edge (measured
%               against an analytic two-plate case and a pure-Python ray cast
%               in tests/tNgfGeom.m).
%     lut       Look-up table over incidence directions: the ray cast is
%               evaluated once at the V vertices of a unit icosphere
%               ('LutResolution' V = 12, 42, 162, 642, 2562, ...) and
%               interpolated barycentrically inside the spherical triangle
%               containing d (C^0). Options.geom.lutResolution is the tier
%               default (642). The table is returned in info.lut and can be
%               passed back as 'Lut' to skip the build.
%
% INPUT ARGUMENTS
%   Name       Size  Type                          Units  Description
%   ----       ----  ----                          -----  -----------
%   macro      1x1   NonGravForces.geom.MacroModel  -      With vertices on every plate (hasVertices)
%   dir_body   3xK   double                         -      Unit directions to the source, body frame;
%                                                          [] with Method "lut" returns only the table
%
% NAME-VALUE ARGUMENTS
%   Name            Size  Type    Units  Description
%   ----            ----  ----    -----  -----------
%   Method          1x1   string  -      "raycast" (default) | "lut"
%   Samples         1x1   double  -      grid points per plate edge for the ray cast (default 40)
%   LutResolution   1x1   double  -      icosphere vertex count 12, 42, 162, 642, 2562, ...;
%                                         default: the resolution of 'Lut' if given, else 642
%   Lut             1x1   struct  -      a previously returned info.lut for this macro model
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   f        NxK   double  -      Exposed fraction per plate and direction, in [0, 1]
%   info     1x1   struct  -      .method, .samples, .lut (struct .verts 3xV, .faces
%                                 3xF, .f NxV, .resolution, .samples, .macroName; only
%                                 for "lut"), .provenance
%
% ACCURACY
%   Ray cast: area-fraction error <= ~2/Samples at a shadow edge (a
%   straight edge crossing an m x m grid misclassifies at most ~2m of m^2
%   cells): 40 samples give ~5e-2 worst case and ~1e-2 typical; measured
%   values per run in tests/last_run.txt, including the convergence over
%   Samples and LutResolution (T4). LUT: the barycentric interpolation adds
%   the table's angular discretisation (~6 deg between vertices at 642,
%   ~3 deg at 2562); the test records the LUT-vs-raycast difference. The
%   plate polygons themselves are approximations of the exterior geometry
%   (geom.graceFO: ~1-4 cm).
%
% ERRORS
%   NonGravForces:selfShadow:inputMacro       not a MacroModel, or without vertices on every plate
%   NonGravForces:selfShadow:inputShape       dir_body not 3xK, not unit, or NaN
%   NonGravForces:selfShadow:inputResolution  LutResolution not of the form 10*4^n + 2
%   NonGravForces:selfShadow:inputLut         'Lut' built for another macro model or resolution
%
% REFERENCES
%   Wöske et al. (2019) Adv. Space Res. 63:1096 (finite-element ray tracing vs plates);
%   HANDOVER D.6 (LUT strategy), C.3.
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C");
%   f = NonGravForces.geom.selfShadow(mm, [0; sind(70); cosd(70)])   % low elevation: the port flank shadows the nadir plate
%
% See also NonGravForces.geom.MacroModel, NonGravForces.geom.graceFO, NonGravForces.env.radiationSources.
%
% KEYWORDS self-shadowing ray cast look-up table icosphere occlusion
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    macro
    dir_body double
    options.Method (1,1) string {mustBeMember(options.Method, ["raycast","lut"])} = "raycast"
    options.Samples (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.Samples, 2)} = 40
    options.LutResolution double = []
    options.Lut struct = struct([])
end
unit = "selfShadow";
if ~isa(macro, 'NonGravForces.geom.MacroModel') || ~macro.hasVertices
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel with vertices_m on every plate -> geom.graceFO and geom.boxWing (bus) provide them");
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
info.method = options.Method;
info.samples = options.Samples;
switch options.Method
    case "raycast"
        f = localRaycast(macro, dir_body, options.Samples);
    case "lut"
        if isempty(options.LutResolution)
            if ~isempty(options.Lut) && isfield(options.Lut, 'resolution'), options.LutResolution = options.Lut.resolution; else, options.LutResolution = 642; end
        end
        if ~isscalar(options.LutResolution)
            NonGravForces.fail(unit, "inputResolution", "LutResolution must be a scalar -> 12, 42, 162, 642, 2562, ...");
        end
        level = log(max(options.LutResolution - 2, 1) / 10) / log(4);
        if abs(level - round(level)) > 1e-9 || level < 0 || level > 7
            NonGravForces.fail(unit, "inputResolution", "LutResolution %g is not 10*4^n + 2 (12, 42, 162, 642, 2562, ...) -> pick one of those", options.LutResolution);
        end
        if ~isempty(options.Lut)
            lut = options.Lut;
            if ~(isfield(lut, 'macroName') && lut.macroName == macro.name && isfield(lut, 'resolution') && lut.resolution == options.LutResolution && isequal(size(lut.f), [macro.n, options.LutResolution]))
                NonGravForces.fail(unit, "inputLut", "the Lut was built for '%s' at resolution %g -> build one for this model and resolution", string(lut.macroName), lut.resolution);
            end
        else
            [verts, faces] = icosphere(round(level));
            lut = struct('verts', verts, 'faces', faces, 'f', localRaycast(macro, verts, options.Samples), ...
                         'resolution', options.LutResolution, 'samples', options.Samples, 'macroName', macro.name);
        end
        info.lut = lut;
        if isempty(dir_body)
            f = lut;                                             % table only
        else
            f = localInterpolate(lut, dir_body);
        end
end
info.provenance = NonGravForces.provenance("geom.selfShadow", 'Options', struct('method', options.Method, 'samples', options.Samples, 'lutResolution', options.LutResolution));
end

% ------------------------------------------------------------------------- ray cast
function f = localRaycast(macro, d, m)
% d: 3xK unit directions to the source. Sample points on each plate, occlusion by every other plate.
N = macro.n; K = size(d, 2);
f = zeros(N, K);
[P, w] = localSamples(macro, m);                     % P{k} 3xS points, w{k} 1xS weights (sum 1)
for k = 1:N
    lit = (macro.N(:, k)' * d) > 0;                  % 1xK plates facing the beam
    if ~any(lit), continue; end
    dl = d(:, lit); Kl = size(dl, 2);
    S = size(P{k}, 2);
    occluded = false(S, Kl);
    for j = 1:N
        if j == k, continue; end
        nj = macro.N(:, j); Vj = macro.V{j};
        denom = nj' * dl;                                                       % 1xKl
        hit = abs(denom) > 1e-12;
        if ~any(hit), continue; end
        % t = n_j . (p_j - P) / (n_j . d), one value per (sample, direction)
        t = (nj' * (Vj(:, 1) - P{k})).' ./ denom;                              % SxKl
        ahead = t > 1e-9;
        % intersection points and the point-in-polygon test in the plane of j
        e1 = Vj(:, 2) - Vj(:, 1); e1 = e1 / norm(e1); e2 = cross(nj, e1);
        px = e1' * (Vj - Vj(:, 1)); py = e2' * (Vj - Vj(:, 1));                 % polygon in 2-D
        Xk = reshape(P{k}, 3, S, 1) + reshape(t, 1, S, Kl) .* reshape(dl, 3, 1, Kl);   % 3xSxKl hit points
        Xr = Xk - Vj(:, 1);
        hx = reshape(sum(Xr .* e1, 1), S, Kl); hy = reshape(sum(Xr .* e2, 1), S, Kl);
        inside = localInPolygon(hx, hy, px, py);
        occluded = occluded | (ahead & inside & hit);
    end
    % column-wise sum (not a BLAS matrix product): the result must not depend on how many directions are batched;
    % fully occluded and fully exposed columns are set exactly (the weights sum to 1 only to rounding)
    frac = sum(w{k}' .* double(occluded), 1);
    frac(all(occluded, 1)) = 1;
    frac(~any(occluded, 1)) = 0;
    f(k, lit) = 1 - frac;
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

% ------------------------------------------------------------------------- LUT interpolation
function f = localInterpolate(lut, d)
% Barycentric interpolation of the table inside the icosphere triangle that contains each direction.
K = size(d, 2);
A = lut.verts(:, lut.faces(1, :)); B = lut.verts(:, lut.faces(2, :)); C = lut.verts(:, lut.faces(3, :));
f = zeros(size(lut.f, 1), K);
for k = 1:K
    dk = d(:, k);
    % barycentric coordinates of the central projection of dk onto each face plane
    nrm = cross(B - A, C - A, 1);                       % face normals (outward, consistent winding)
    s = sum(nrm .* A, 1) ./ max(sum(nrm .* dk, 1), 1e-300);   % ray scale to the plane
    X = s .* dk;                                          % 3xF intersection points
    % barycentric via areas
    w0 = sum(cross(B - X, C - X, 1) .* nrm, 1); w1 = sum(cross(C - X, A - X, 1) .* nrm, 1); w2 = sum(cross(A - X, B - X, 1) .* nrm, 1);
    tot = sum(nrm.^2, 1);
    bc = [w0; w1; w2] ./ tot;
    ok = all(bc >= -1e-9, 1) & s > 0;
    idx = find(ok, 1);
    if isempty(idx), [~, idx] = max(min(bc, [], 1)); end
    b = max(0, bc(:, idx)); b = b / sum(b);
    f(:, k) = lut.f(:, lut.faces(1, idx)) * b(1) + lut.f(:, lut.faces(2, idx)) * b(2) + lut.f(:, lut.faces(3, idx)) * b(3);
end
end
