function mm = grace(satellite, options)
%GRACE GRACE-A/B macro model from the GRACE Product Specification Document (Bettadpur 2012), Sec. I.5.3.
%
%   mm = NonGravForces.geom.grace("A")
%   mm = NonGravForces.geom.grace("B", 'SolarArrayOperating', false)
%
%   Eight plates in the Satellite Frame (+X to the K-band horn, +Z nadir, +Y
%   completes), which the document says "may be regarded as the same" as the
%   Science Reference Frame for this purpose (Sec. I.5.3). Both satellites are
%   built from ONE table: the document gives the macro model "for only one
%   satellite" because the two are identical (Sec. I.5.3.3).
%
%   SURFACE TABLE (Sec. I.5.3.5, PDF p. 23), transcribed from the document
%   itself. It is the same table, digit for digit, as JPL D-56935 Table 5 that
%   geom.graceFO transcribes, and a test compares the two transcriptions. The
%   note-2 alternative for a non-operating array (absVis 0.72 with the printed
%   diffuse 0.30 and specular 0.05) sums to 1.07; as in geom.graceFO, diffVis
%   becomes 0.23 so that the budget closes.
%
%   The Boom (0.0461901 m^2, SiOx/Kapton) is excluded. Its axis is known here
%   (note 4: the area is its projection along any direction in the X-Y plane,
%   so the cylinder lies along Z; Figs 3-4 draw it below the nadir panel), but
%   no kernel consumes a "cylinder" element and a flat plate cannot stand in
%   for one. It is 0.27 % of the area.
%
%   PLATE POLYGONS (Sec. I.5.3.4, Figs 3-4, PDF pp. 21-22). Built as in
%   geom.graceFO, from R-5's own printed dimensions: length 3123 mm, height
%   720 mm, zenith width 694 mm, arrays along the printed 50 deg slant. The
%   flank meets the nadir plane at half-width 0.9512 m, not at the printed
%   1944/2 = 0.972 m (the handbook's dimensions do not close at the 2 cm level,
%   the same inconsistency as GRACE-FO's), and continues beyond it by the
%   printed underside width 0.2282913 / 3.123 = 0.0731 m. That puts the lower
%   edge of the overhang 0.0731 sin 50 = 56.0 mm below the nadir plane, which
%   is the 56 mm Fig 3 prints: an independent check the tests assert.
%   The origin is the C.o.G. of Fig 3: the nadir plane lies 275 mm below it
%   and the zenith 720 - 275 = 445 mm above it, which the thruster table
%   (Sec. I.5.3.2: A14/A22 at Z = +275 mm, A12/A24 at Z = -444 mm) confirms.
%   Polygon areas are within 2.2 % of the printed ones; THE PRINTED AREAS STAY
%   THE FORCE AUTHORITY, the polygons serve geom.selfShadow and the torque arms.
%
%   HEAT CAPACITIES: the Wöske (2019) Table 6.1 layer stacks geom.graceFO uses,
%   per material; R-5 prints the same material on every plate as D-56935.
%
%   MASS: the launch mass of Sec. I.5.3.3, 487.2 kg for both. The document
%   points to MAS1A/MAS1B for the mass during the mission, and for classic
%   GRACE that product carries the SPACECRAFT mass (476.85 kg for GRACE-A on
%   2008-01-01), unlike GRACE-FO's tank mass. io.massModel serves GRACE-FO only
%   and refuses anything at or below 500 kg (HANDOVER O-46).
%
%   NOT IN THE DOCUMENT, hence not set: the inertia tensor (NaN; the gravity-
%   gradient torque is then refused by kernel.bodyTorques, never guessed) and
%   the CoG -> proof-mass offset (NaN). The centre of mass is the origin.
%
% INPUT ARGUMENTS
%   Name       Size Type   Units Description
%   ----       ---- ----   ----- -----------
%   satellite  1x1  string -     "A" (GRACE-1, FM-1) or "B" (GRACE-2, FM-2)
%
% NAME-VALUE ARGUMENTS
%   SolarArrayOperating 1x1 logical  true: absVis 0.65 / diffVis 0.30 (closes).
%                                    false: absVis 0.72 with diffVis 0.23 so
%                                    that the budget closes (the printed 0.30
%                                    would sum to 1.07 - a handbook defect).
%
% OUTPUT ARGUMENTS
%   mm         1x1  NonGravForces.geom.MacroModel with .V (plate polygons),
%                   .heatCapacity, .mass_kg (the launch mass,
%                   NonGravForces.const.grace), .inertia NaN(3)
%
% ERRORS
%   MATLAB:validators:mustBeMember   satellite is not "A" or "B" (arguments block)
%   plus the contract errors of geom.MacroModel, which validates the table.
%
% ACCURACY
%   Areas, normals and optical coefficients are transcriptions, no
%   approximation. Polygons follow from the exterior dimensions to ~2 cm (the
%   document's own dimensions are inconsistent at that level); they enter
%   torques and self-shadowing only. The side silhouette of the polygons plus
%   the printed boom area is within 0.3 % of the 2.4633921 m^2 Fig 4 prints;
%   the difference is the dome, notch and recesses the figure draws and the
%   plates do not. Wöske et al. (2019) quote ~10 % force error of an 8-plate
%   model against finite-element ray-tracing for GRACE. The launch mass
%   ignores the gas used since launch (476.85 kg on 2008-01-01 for GRACE-A,
%   -2.1 %); pass 'Mass' to Model for the mission value.
%
% REFERENCES
%   Bettadpur S. (2012) GRACE Product Specification Document, GRACE 327-720
%     (CSR-GR-03-02), Rev 4.6 draft, 29 May 2012: Sec. I.5.3.2 (thrusters),
%     I.5.3.3 (mass), I.5.3.4 Figs 3-4 (exterior dimensions), I.5.3.5 (surface
%     properties), PDF pp. 20-23.
%   Wöske F. (2019) PhD thesis, Univ. Bremen, Table 6.1 (wall layer stacks).
%
% EXAMPLE
%   mm = NonGravForces.geom.grace("A"); mm.silhouetteDefect()
%
% See also NonGravForces.geom.graceFO, NonGravForces.geom.MacroModel, NonGravForces.geom.selfShadow.
%
% KEYWORDS GRACE macro model product specification panels classic
% Created 2026-09-22 by Claude (Anthropic) for Matthias Weigelt.
arguments
    satellite (1,1) string {mustBeMember(satellite, ["A","B"])}
    options.SolarArrayOperating (1,1) logical = true
end
s50 = 0.766044; c50 = 0.642787;
if options.SolarArrayOperating, aV = 0.65; dV = 0.30; else, aV = 0.72; dV = 0.23; end
Lx = 3.123; Hz = 0.720; Wt = 0.694;  % Figs 3-4 [m]
zNadir = 0.275; zZenith = zNadir - Hz;  % Fig 3: C.o.G. 275 mm above the nadir plane
wInner = 0.2282913 / Lx;            % underside of the overhang (printed area / length)
tCross = Hz / s50;                  % slant length from the zenith edge to the nadir plane
yNadir = Wt / 2 + tCross * c50;     % nadir half-width where the flank crosses the nadir plane (0.9512 m)
% Wöske (2019) Table 6.1 layer stacks, cp*rho*h summed [J/(m^2 K)], as in geom.graceFO
hcArray  = 540*2000*0.0025 + 1000*30*0.070 + 900*45*0.030;    % solar panel + insulation + honeycomb
hcKapton = 1095*1400*0.0001 + 1000*33*0.075 + 900*45*0.030;   % Kapton + insulation + honeycomb
hcTeflon = 1095*2000*0.0001;                                   % teflon (FEP) radiator foil
% R-5 Sec. I.5.3.5: name, area, normal, absVis specVis diffVis, emisIR specIR diffIR, group, heat capacity
T = { ...
 "Front",           0.9551567, [ 1; 0; 0],    0.34,0.40,0.26, 0.62,0.23,0.15, "bus",   hcKapton; ...
 "Rear",            0.9551567, [-1; 0; 0],    0.34,0.40,0.26, 0.62,0.23,0.15, "bus",   hcKapton; ...
 "StarboardOuter",  3.1554792, [0; s50;-c50], aV,0.05,dV,     0.81,0.03,0.16, "array", hcArray; ...
 "StarboardInner",  0.2282913, [0;-s50; c50], 0.34,0.40,0.26, 0.62,0.23,0.15, "bus",   hcKapton; ...
 "PortOuter",       3.1554792, [0;-s50;-c50], aV,0.05,dV,     0.81,0.03,0.16, "array", hcArray; ...
 "PortInner",       0.2282913, [0; s50; c50], 0.34,0.40,0.26, 0.62,0.23,0.15, "bus",   hcKapton; ...
 "Nadir",           6.0711120, [0; 0; 1],     0.12,0.68,0.20, 0.75,0.19,0.06, "bus",   hcTeflon; ...
 "Zenith",          2.1673620, [0; 0;-1],     aV,0.05,dV,     0.81,0.03,0.16, "array", hcArray};
% Plate polygons (body frame). Cross-section in the y-z plane: zenith edge at (+-Wt/2, zZenith),
% flank direction (sgn*c50, s50), nadir edge at (+-yNadir, zNadir), overhang end wInner beyond it.
polys = cell(8, 1);
sec = [ -yNadir, zNadir;  yNadir, zNadir;  Wt/2, zZenith;  -Wt/2, zZenith ]';     % 2x4 (y; z)
polys{1} = [ (Lx/2) * ones(1, 4); sec(1, :); sec(2, :) ];
polys{2} = [-(Lx/2) * ones(1, 4); sec(1, :); sec(2, :) ];
for sgn = [1 -1]
    a = [sgn * Wt / 2; zZenith];                       % zenith edge
    b = [sgn * yNadir; zNadir];                        % nadir edge (flank crosses the nadir plane)
    e = b + wInner * [sgn * c50; s50];                 % end of the overhang
    outer = [ -Lx/2, Lx/2, Lx/2, -Lx/2; a(1), a(1), e(1), e(1); a(2), a(2), e(2), e(2) ];
    inner = [ -Lx/2, Lx/2, Lx/2, -Lx/2; b(1), b(1), e(1), e(1); b(2), b(2), e(2), e(2) ];
    if sgn == 1, polys{3} = outer; polys{4} = inner; else, polys{5} = outer; polys{6} = inner; end
end
polys{7} = [ -Lx/2, Lx/2, Lx/2, -Lx/2; -yNadir, -yNadir, yNadir, yNadir; zNadir * ones(1, 4) ];    % nadir (+z)
polys{8} = [ -Lx/2, Lx/2, Lx/2, -Lx/2;  Wt/2,  Wt/2, -Wt/2, -Wt/2; zZenith * ones(1, 4) ];       % zenith (-z)
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
c = NonGravForces.const();
mm = NonGravForces.geom.MacroModel(panels, 'Name', "GRACE-" + satellite, 'Mass_kg', c.grace.("mass" + satellite), ...
    'CentreOfMass', zeros(3, 1), 'ProofMassOffset', nan(3, 1), ...
    'Reference', "Bettadpur 2012, GRACE 327-720 Rev 4.6, Sec. I.5.3.3-I.5.3.5; polygons from Figs 3-4; heat capacities Wöske 2019 Table 6.1");
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
