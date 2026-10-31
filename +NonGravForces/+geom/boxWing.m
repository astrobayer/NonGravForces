function mm = boxWing(dims, optics, options)
%BOXWING Generic box-wing macro model: a rectangular bus of six plates and optional solar-array plates.
%
%   mm = NonGravForces.geom.boxWing([Lx Ly Lz], optics)
%   mm = NonGravForces.geom.boxWing([Lx Ly Lz], optics, 'Array', arr, 'Name', "MySat", 'Mass_kg', 480)
%
%   Builds a NonGravForces.geom.MacroModel for a body-frame box centred at the
%   origin with edge lengths Lx, Ly, Lz (the six plates carry outward normals
%   +-x, +-y, +-z, centroids at the face centres and four-corner polygons for
%   self-shadowing) and, optionally, one or more solar-array plates given by
%   area, unit normal and centroid (each array is one plate whose back side is
%   not modelled; add a second entry with the opposite normal for a two-sided
%   wing). Every optical coefficient is the caller's: nothing here is a
%   spacecraft.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   dims     1x3   double  m      [Lx Ly Lz] edge lengths of the bus, positive
%   optics   1x1   struct  -      .bus and (if 'Array' is given) .array, each with
%                                 fields vis (1x3 [abs spec diff], sums to 1) and
%                                 ir (1x3 [emis spec diff], sums to 1)
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   Array    1xP   struct  -      per array plate: .area_m2, .normal (3x1 unit),
%                                 .centroid_m (3x1), optional .name; default none
%   Name     1x1   string  -      Model name (default "boxWing")
%   Mass_kg  1x1   double  kg     Spacecraft mass; default NaN (unknown: a Model then needs 'Mass')
%   HeatCapacity_J_m2K 1x2 double J/(m^2 K)  [bus array] heat capacity per area
%                                 for thermalState transient; default [NaN NaN]
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   mm       1x1   NonGravForces.geom.MacroModel  6 + P plates, polygons on the bus plates
%
% ACCURACY
%   Exact for a box: areas, normals, centroids and polygons follow from the
%   dimensions. Array plates have no polygon (their shape is not given), so a
%   model with arrays cannot be ray cast (selfShadow refuses); the bus alone can.
%
% ERRORS
%   NonGravForces:boxWing:inputDims     dims not three positive finite lengths
%   NonGravForces:boxWing:inputOptics   optics lacks .bus.vis/.bus.ir (or .array.vis/.array.ir when arrays
%                                       are given), or a vis/ir triple does not have three elements
%   NonGravForces:boxWing:inputArray    an array entry lacks area_m2/normal/centroid_m, its area is not a
%                                       positive scalar, or its normal is not 3-element unit (1e-5)
%   NonGravForces:MacroModel:*          from the constructor: inputEnergy when a vis/ir triple does not
%                                       sum to 1, inputPanels for a non-positive HeatCapacity_J_m2K or Mass_kg
%
% REFERENCES
%   None beyond the caller's numbers: the geometry is a box and every coefficient is an input.
%
% EXAMPLE
%   o = struct('bus', struct('vis', [0.34 0.40 0.26], 'ir', [0.62 0.23 0.15]));
%   mm = NonGravForces.geom.boxWing([3 1.9 0.7], o); mm.n            % 6
%
% See also NonGravForces.geom.MacroModel, NonGravForces.geom.graceFO.
%
% KEYWORDS box-wing macro model factory
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    dims (1,3) double
    optics (1,1) struct
    options.Array struct = struct([])
    options.Name (1,1) string = "boxWing"
    options.Mass_kg (1,1) double = NaN
    options.HeatCapacity_J_m2K (1,2) double = [NaN NaN]
end
unit = "boxWing";
if any(~isfinite(dims)) || any(dims <= 0)
    NonGravForces.fail(unit, "inputDims", "dims must be three positive finite lengths, got %s -> [Lx Ly Lz] in metres", mat2str(dims));
end
if ~(isfield(optics, 'bus') && all(isfield(optics.bus, ["vis","ir"])))
    NonGravForces.fail(unit, "inputOptics", "optics needs .bus.vis and .bus.ir ([abs spec diff], [emis spec diff]) -> see the help");
end
P = numel(options.Array);
if P > 0 && ~(isfield(optics, 'array') && all(isfield(optics.array, ["vis","ir"])))
    NonGravForces.fail(unit, "inputOptics", "arrays are given, so optics needs .array.vis and .array.ir -> add them");
end
Lx = dims(1); Ly = dims(2); Lz = dims(3);
faces = { ...
    "+x", [ 1;0;0], Ly*Lz, [ Lx/2;0;0]; "-x", [-1;0;0], Ly*Lz, [-Lx/2;0;0]; ...
    "+y", [0; 1;0], Lx*Lz, [0; Ly/2;0]; "-y", [0;-1;0], Lx*Lz, [0;-Ly/2;0]; ...
    "+z", [0;0; 1], Lx*Ly, [0;0; Lz/2]; "-z", [0;0;-1], Lx*Ly, [0;0;-Lz/2]};
tmpl = struct('name',"",'kind',"plate",'area_m2',0,'normal',zeros(3,1),'centroid_m',zeros(3,1), ...
    'absVis',0,'specVis',0,'diffVis',0,'emisIR',0,'specIR',0,'diffIR',0,'group',"", ...
    'vertices_m',zeros(3,0),'heatCapacity_J_m2K',NaN);
panels = repmat(tmpl, 1, 6 + P);
half = dims(:) / 2;
for k = 1:6
    nrm = faces{k, 2}; ctr = faces{k, 4};
    panels(k).name = "bus" + faces{k, 1}; panels(k).area_m2 = faces{k, 3}; panels(k).normal = nrm; panels(k).centroid_m = ctr;
    panels(k) = localOptics(panels(k), optics.bus, unit, "bus"); panels(k).group = "bus";
    panels(k).heatCapacity_J_m2K = options.HeatCapacity_J_m2K(1);
    % corners: the two in-plane axes span the face
    ax = find(nrm ~= 0); others = setdiff(1:3, ax);
    u = zeros(3, 1); u(others(1)) = half(others(1)); w = zeros(3, 1); w(others(2)) = half(others(2));
    corners = [ctr - u - w, ctr + u - w, ctr + u + w, ctr - u + w];
    e1 = corners(:, 2) - corners(:, 1); e2 = corners(:, 4) - corners(:, 1);
    if dot(cross(e1, e2), nrm) < 0, corners = corners(:, [1 4 3 2]); end   % counter-clockwise about the outward normal
    panels(k).vertices_m = corners;
end
for p = 1:P
    a = options.Array(p);
    if ~all(isfield(a, ["area_m2","normal","centroid_m"])) || ~(isscalar(a.area_m2) && a.area_m2 > 0) || numel(a.normal) ~= 3 || numel(a.centroid_m) ~= 3 || abs(norm(a.normal) - 1) > 1e-5
        NonGravForces.fail(unit, "inputArray", "Array(%d) needs positive area_m2, a 3x1 unit normal and a 3x1 centroid_m -> check the entry", p);
    end
    k = 6 + p;
    if isfield(a, 'name') && ~isempty(a.name), panels(k).name = string(a.name); else, panels(k).name = "array" + p; end
    panels(k).area_m2 = a.area_m2; panels(k).normal = a.normal(:); panels(k).centroid_m = a.centroid_m(:);
    panels(k) = localOptics(panels(k), optics.array, unit, "array"); panels(k).group = "array";
    panels(k).heatCapacity_J_m2K = options.HeatCapacity_J_m2K(2);
end
mm = NonGravForces.geom.MacroModel(panels, 'Name', options.Name, 'Mass_kg', options.Mass_kg, 'Reference', "boxWing factory: caller dimensions and optics");
end

function p = localOptics(p, o, unit, which)
if numel(o.vis) ~= 3 || numel(o.ir) ~= 3
    NonGravForces.fail(unit, "inputOptics", "optics.%s.vis and .ir must be 1x3 -> [abs spec diff] and [emis spec diff]", which);
end
p.absVis = o.vis(1); p.specVis = o.vis(2); p.diffVis = o.vis(3);
p.emisIR = o.ir(1); p.specIR = o.ir(2); p.diffIR = o.ir(3);
end
