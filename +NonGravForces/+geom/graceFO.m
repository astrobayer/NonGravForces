function mm = graceFO(satellite, options)
%GRACEFO GRACE-FO macro model from JPL D-56935 (Wen et al. 2019), Tables 4 and 5.
%
%   mm = NonGravForces.geom.graceFO("C")
%   mm = NonGravForces.geom.graceFO("D", 'SolarArrayOperating', false)
%
%   Eight plates in the Science Reference Frame (+X to the K-band horn, +Z
%   nadir, +Y completes). Side panels at 50 deg (normals sin/cos 50 exactly
%   as printed). The body is NOT convex: the overhanging flanks expose their
%   undersides, so sum(A n) ~= 0 and self-shadowing is physical.
%
%   The S-band boom (Table 5 note 4: projection area 0.0461901 m^2, a
%   cylinder) is excluded: a flat plate cannot represent it, and its axis is
%   in no source in hand (rule R3). It is 0.27 % of the area.
%
%   PLATE POLYGONS AND CENTROIDS (2026-09-07). Each plate carries four
%   corners built from the exterior dimensions of Figs 2-3 so that the eight
%   polygons close on each other without intersecting: length 3.1225 m,
%   height 0.7210 m, zenith width 0.6950 m, array flanks along the printed
%   50 deg slant from the zenith edge; the flank meets the nadir plane after
%   0.7210 / sin 50 = 0.9412 m, which puts the nadir edge at half-width
%   0.9525 m (the printed nadir width 1.9440 m would let the plate protrude
%   2 cm through the flank, which is what a ray cast then reports as
%   spurious occlusion, measured 2026-09-07); beyond the nadir plane the
%   flank continues by the printed underside width 0.2282913 / 3.1225 =
%   0.0731 m, whose lower face is the inner flank. Polygon areas are within
%   2 % of the printed areas (nadir -2.0 %, front/rear -1.9 %, arrays
%   +0.4 %); THE PRINTED AREAS STAY THE FORCE AUTHORITY, the polygons serve
%   geom.selfShadow only. Centroids are the exact area centroids of these
%   polygons (the LP-0 values placed the flank centroids 0.39 m too high,
%   an error the polygons exposed; HANDOVER Part I 2026-09-07).
%
%   HEAT CAPACITIES (for geom.thermalState "transient") are Wöske (2019)
%   Table 6.1 layer stacks, which that thesis calls educated guesses: solar
%   panel walls (arrays, zenith) 6015, Kapton/insulation/honeycomb walls
%   (front, rear, inner flanks) 3843, teflon radiator (nadir) 219 J/(m^2 K).
%   Since 2026-10-31 the standard preset's thermal model is transient, so these
%   values enter every standard Model. They were FITTED against ACT1B before that
%   switch (tools/heatCapacityFit.m, one factor per group of 0.25-4, 21 days of
%   2020-01 C and D, run 41140282748) and KEPT: the best factors improve rms_y +
%   rms_z by 1 % over Woeske's, C and D prefer different corners of the grid, and
%   the kapton factor is not determined at all. The radiator factor moves the
%   per-day z scale from 1.08 (x1) to 0.96 (x4) at the same residual, so S_zz is
%   only known to that band (HANDOVER, "The checkAlignment lags").
%
% INPUT ARGUMENTS
%   Name       Size Type   Units Description
%   ----       ---- ----   ----- -----------
%   satellite  1x1  string -     "C" (FM1, 601.214 kg) or "D" (FM2, 601.209 kg)
%
% NAME-VALUE ARGUMENTS
%   SolarArrayOperating 1x1 logical  true: absVis 0.65 / diffVis 0.30 (closes).
%                                    false: absVis 0.72 with diffVis 0.23 so
%                                    that the budget closes (the printed 0.30
%                                    would sum to 1.07 - a handbook defect).
%
% OUTPUT ARGUMENTS
%   mm         1x1  NonGravForces.geom.MacroModel with .inertia (Table 4),
%                   .V (plate polygons), .heatCapacity, .mass_kg (the handbook
%                   mass of that satellite, NonGravForces.const.graceFO)
%
% ERRORS
%   MATLAB:validators:mustBeMember   satellite is not "C" or "D" (arguments block)
%   plus the contract errors of geom.MacroModel, which validates the table.
%
% ACCURACY
%   Areas and normals are transcriptions, no approximation. Polygons and
%   centroids follow from the exterior dimensions to ~2 cm (the handbook's
%   own dimensions are inconsistent at that level); they enter torques and
%   self-shadowing only. Wöske et al. (2019) quote ~10 % force error of an
%   8-plate model against finite-element ray-tracing for GRACE. Heat
%   capacities are Wöske's guesses (Table 6.1): the transient thermal model
%   inherits their status.
%
% REFERENCES
%   Wen et al. (2019) JPL D-56935 Tab. 4, Tab. 5, Sec. 2.5, Sec. 3.2.3.
%   Wöske F. (2019) PhD thesis, Univ. Bremen, Table 6.1 (wall layer stacks).
%   Siemes et al. (2023) JSWSC 13:16, tuned array coefficients (opt-in later).
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C"); mm.silhouetteDefect()
%
% See also NonGravForces.geom.MacroModel, NonGravForces.geom.selfShadow, NonGravForces.geom.thermalState.
%
% KEYWORDS GRACE-FO macro model handbook panels
% Created 2026-08-27 by Claude Fable 5 (Anthropic); polygons, centroids and heat capacities 2026-09-07 by Claude Fable 5.1.
arguments
    satellite (1,1) string {mustBeMember(satellite, ["C","D"])}
    options.SolarArrayOperating (1,1) logical = true
end
s50 = 0.766044; c50 = 0.642787;
if options.SolarArrayOperating, aV = 0.65; dV = 0.30; else, aV = 0.72; dV = 0.23; end
Lx = 3.1225; Hz = 0.7210; Wt = 0.6950;
wInner = 0.2282913 / Lx;            % underside of the overhang (printed area / length)
tCross = Hz / s50;                  % slant length from the zenith edge to the nadir plane
yNadir = Wt / 2 + tCross * c50;     % nadir half-width where the flank crosses the nadir plane (0.9525 m)
% Wöske (2019) Table 6.1 layer stacks, cp*rho*h summed [J/(m^2 K)]
hcArray  = 540*2000*0.0025 + 1000*30*0.070 + 900*45*0.030;    % solar panel + insulation + honeycomb
hcKapton = 1095*1400*0.0001 + 1000*33*0.075 + 900*45*0.030;   % Kapton + insulation + honeycomb
hcTeflon = 1095*2000*0.0001;                                   % teflon radiator foil
% name, area, normal, absVis specVis diffVis, emisIR specIR diffIR, group, heat capacity
T = { ...
 "Front",           0.9551567, [ 1; 0; 0],    0.34,0.40,0.26, 0.62,0.23,0.15, "bus",   hcKapton; ...
 "Rear",            0.9551567, [-1; 0; 0],    0.34,0.40,0.26, 0.62,0.23,0.15, "bus",   hcKapton; ...
 "StarboardOuter",  3.1554792, [0; s50;-c50], aV,0.05,dV,     0.81,0.03,0.16, "array", hcArray; ...
 "StarboardInner",  0.2282913, [0;-s50; c50], 0.34,0.40,0.26, 0.62,0.23,0.15, "bus",   hcKapton; ...
 "PortOuter",       3.1554792, [0;-s50;-c50], aV,0.05,dV,     0.81,0.03,0.16, "array", hcArray; ...
 "PortInner",       0.2282913, [0; s50; c50], 0.34,0.40,0.26, 0.62,0.23,0.15, "bus",   hcKapton; ...
 "Nadir",           6.0711120, [0; 0; 1],     0.12,0.68,0.20, 0.75,0.19,0.06, "bus",   hcTeflon; ...
 "Zenith",          2.1673620, [0; 0;-1],     aV,0.05,dV,     0.81,0.03,0.16, "array", hcArray};
% Plate polygons (body frame). Cross-section in the y-z plane: zenith edge at (+-Wt/2, -Hz/2),
% flank direction (sgn*c50, s50), nadir edge at (+-yNadir, Hz/2), overhang end wInner beyond it.
polys = cell(8, 1);
sec = [ -yNadir, Hz/2;  yNadir, Hz/2;  Wt/2, -Hz/2;  -Wt/2, -Hz/2 ]';            % 2x4 (y; z)
polys{1} = [ (Lx/2) * ones(1, 4); sec(1, :); sec(2, :) ];
polys{2} = [-(Lx/2) * ones(1, 4); sec(1, :); sec(2, :) ];
for sgn = [1 -1]
    a = [sgn * Wt / 2; -Hz / 2];                       % zenith edge
    b = [sgn * yNadir; Hz / 2];                        % nadir edge (flank crosses the nadir plane)
    e = b + wInner * [sgn * c50; s50];                 % end of the overhang
    outer = [ -Lx/2, Lx/2, Lx/2, -Lx/2; a(1), a(1), e(1), e(1); a(2), a(2), e(2), e(2) ];
    inner = [ -Lx/2, Lx/2, Lx/2, -Lx/2; b(1), b(1), e(1), e(1); b(2), b(2), e(2), e(2) ];
    if sgn == 1, polys{3} = outer; polys{4} = inner; else, polys{5} = outer; polys{6} = inner; end
end
polys{7} = [ -Lx/2, Lx/2, Lx/2, -Lx/2; -yNadir, -yNadir, yNadir, yNadir; Hz/2 * ones(1, 4) ];   % nadir (+z)
polys{8} = [ -Lx/2, Lx/2, Lx/2, -Lx/2;  Wt/2,  Wt/2, -Wt/2, -Wt/2; -Hz/2 * ones(1, 4) ];       % zenith (-z)
n = size(T, 1);
panels = repmat(struct('name',"",'kind',"plate",'area_m2',0,'normal',zeros(3,1),'centroid_m',zeros(3,1), ...
    'absVis',0,'specVis',0,'diffVis',0,'emisIR',0,'specIR',0,'diffIR',0,'group',"", ...
    'vertices_m',zeros(3,0),'heatCapacity_J_m2K',NaN), 1, n);
for k = 1:n
    panels(k).name = T{k,1}; panels(k).area_m2 = T{k,2}; panels(k).normal = T{k,3};
    panels(k).absVis = T{k,4}; panels(k).specVis = T{k,5}; panels(k).diffVis = T{k,6};
    panels(k).emisIR = T{k,7}; panels(k).specIR = T{k,8}; panels(k).diffIR = T{k,9}; panels(k).group = T{k,10};
    panels(k).heatCapacity_J_m2K = T{k,11};
    v = localWind(polys{k}, panels(k).normal);
    panels(k).vertices_m = v;
    panels(k).centroid_m = localCentroid(v);
end
if satellite == "C"      % Table 4, FM1
    I = [110.491, -1.024,  0.347; -1.024, 580.673, 0.036; 0.347, 0.036, 649.690];
    cg_m = [-0.249; -0.403; 1.266] * 1e-3;   cg2acc_m = [-0.08; 0.00; -0.15] * 1e-3;
else                     % Table 4, FM2
    I = [109.635, -0.497,  0.359; -0.497, 579.759, 0.057; 0.359, 0.057, 648.795];
    cg_m = [-0.285;  1.346; 0.506] * 1e-3;   cg2acc_m = [-0.06; -0.01; -0.04] * 1e-3;
end
c = NonGravForces.const();
mm = NonGravForces.geom.MacroModel(panels, 'Name', "GRACE-" + satellite, 'Inertia', I, 'Mass_kg', c.graceFO.("mass" + satellite), ...
    'CentreOfMass', cg_m, 'ProofMassOffset', cg2acc_m, ...
    'Reference', "Wen et al. 2019, JPL D-56935, Tables 4-5; polygons from Figs 2-3; heat capacities Wöske 2019 Table 6.1");
end

function v = localWind(v, nrm)
% Order the corners counter-clockwise about the outward normal (the construction above
% fixes the geometry; this fixes only the traversal direction, which MacroModel then checks).
ctr = mean(v, 2);
e1 = v(:, 1) - ctr; e1 = e1 / norm(e1); e2 = cross(nrm, e1);
ang = atan2(e2' * (v - ctr), e1' * (v - ctr));
[~, order] = sort(ang);
v = v(:, order);
end

function c = localCentroid(v)
% Area centroid of a planar polygon by a triangle fan from the first corner.
c = zeros(3, 1); a = 0;
for i = 2:size(v, 2) - 1
    tri = [v(:, 1), v(:, i), v(:, i + 1)];
    ai = 0.5 * norm(cross(tri(:, 2) - tri(:, 1), tri(:, 3) - tri(:, 1)));
    c = c + ai * mean(tri, 2); a = a + ai;
end
c = c / a;
end
