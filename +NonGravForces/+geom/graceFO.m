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
%   cylinder) is excluded: a flat plate cannot represent it. It is 0.27 % of
%   the area; LP-4 adds it as kind "cylinder".
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
%   Alpha      1x1  double -     energy accommodation default, 0.85
%                                (March et al. 2019, 400-500 km)
%
% OUTPUT ARGUMENTS
%   mm         1x1  NonGravForces.geom.MacroModel with .inertia (Table 4)
%
% ACCURACY
%   Areas and normals are transcriptions, no approximation. Panel centroids
%   are derived from the exterior dimensions (Figs 2-3) and are accurate to
%   ~1 cm; they enter torques and self-shadowing only. Wöske et al. (2019)
%   quote ~10 % force error of an 8-plate model against finite-element
%   ray-tracing for GRACE.
%
% REFERENCES
%   Wen et al. (2019) JPL D-56935 Tab. 4, Tab. 5, Sec. 2.5, Sec. 3.2.3.
%   Siemes et al. (2023) JSWSC 13:16, tuned array coefficients (opt-in later).
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C"); mm.silhouetteDefect()
%
% See also NonGravForces.geom.MacroModel.
%
% KEYWORDS GRACE-FO macro model handbook panels
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    satellite (1,1) string {mustBeMember(satellite, ["C","D"])}
    options.SolarArrayOperating (1,1) logical = true
    options.Alpha (1,1) double {mustBeInRange(options.Alpha, 0, 1)} = 0.85
end
s50 = 0.766044; c50 = 0.642787;
if options.SolarArrayOperating, aV = 0.65; dV = 0.30; else, aV = 0.72; dV = 0.23; end
Lx = 3.1225; Wy = 1.9440; Hz = 0.7210; Wt = 0.6950;
% name, area, normal, centroid, absVis specVis diffVis, emisIR specIR diffIR, group
T = { ...
 "Front",           0.9551567, [ 1; 0; 0],   [ Lx/2; 0; 0],           0.34,0.40,0.26, 0.62,0.23,0.15, "bus"; ...
 "Rear",            0.9551567, [-1; 0; 0],   [-Lx/2; 0; 0],           0.34,0.40,0.26, 0.62,0.23,0.15, "bus"; ...
 "StarboardOuter",  3.1554792, [0; s50;-c50],[0; (Wy+Wt)/4; -Hz/2],   aV,0.05,dV,     0.81,0.03,0.16, "array"; ...
 "StarboardInner",  0.2282913, [0;-s50; c50],[0; (Wy+Wt)/4;  0],      0.34,0.40,0.26, 0.62,0.23,0.15, "bus"; ...
 "PortOuter",       3.1554792, [0;-s50;-c50],[0;-(Wy+Wt)/4; -Hz/2],   aV,0.05,dV,     0.81,0.03,0.16, "array"; ...
 "PortInner",       0.2282913, [0; s50; c50],[0;-(Wy+Wt)/4;  0],      0.34,0.40,0.26, 0.62,0.23,0.15, "bus"; ...
 "Nadir",           6.0711120, [0; 0; 1],    [0; 0;  Hz/2],           0.12,0.68,0.20, 0.75,0.19,0.06, "bus"; ...
 "Zenith",          2.1673620, [0; 0;-1],    [0; 0; -Hz/2],           aV,0.05,dV,     0.81,0.03,0.16, "array"};
n = size(T,1);
panels = repmat(struct('name',"",'kind',"plate",'area_m2',0,'normal',zeros(3,1),'centroid_m',zeros(3,1), ...
    'absVis',0,'specVis',0,'diffVis',0,'emisIR',0,'specIR',0,'diffIR',0,'alpha',options.Alpha,'group',""), 1, n);
for k = 1:n
    panels(k).name = T{k,1}; panels(k).area_m2 = T{k,2}; panels(k).normal = T{k,3}; panels(k).centroid_m = T{k,4};
    panels(k).absVis = T{k,5}; panels(k).specVis = T{k,6}; panels(k).diffVis = T{k,7};
    panels(k).emisIR = T{k,8}; panels(k).specIR = T{k,9}; panels(k).diffIR = T{k,10}; panels(k).group = T{k,11};
end
if satellite == "C"      % Table 4, FM1
    I = [110.491, -1.024,  0.347; -1.024, 580.673, 0.036; 0.347, 0.036, 649.690];
    cg_m = [-0.249; -0.403; 1.266] * 1e-3;   cg2acc_m = [-0.08; 0.00; -0.15] * 1e-3;
else                     % Table 4, FM2
    I = [109.635, -0.497,  0.359; -0.497, 579.759, 0.057; 0.359, 0.057, 648.795];
    cg_m = [-0.285;  1.346; 0.506] * 1e-3;   cg2acc_m = [-0.06; -0.01; -0.04] * 1e-3;
end
mm = NonGravForces.geom.MacroModel(panels, 'Name', "GRACE-" + satellite, 'Inertia', I, ...
    'CentreOfMass', cg_m, 'ProofMassOffset', cg2acc_m, ...
    'Reference', "Wen et al. 2019, JPL D-56935, Tables 4-5");
end
