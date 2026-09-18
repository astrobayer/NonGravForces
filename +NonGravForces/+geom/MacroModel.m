classdef MacroModel
%MACROMODEL Panel ("macro") geometry of a spacecraft, frozen and validated once.
%
%   mm = NonGravForces.geom.MacroModel(panels, 'Name', "GRACE-C", 'Inertia', I)
%   s  = mm.silhouetteDefect()     % 3x1; zero for a closed convex body
%
% METHODS
%   s = silhouetteDefect(obj)  sum_k A_k n_k, the area-weighted sum of the outward normals, 3x1
%                              in the body frame. It is ZERO for a closed convex body, so it is
%                              the one-line check that a panel set closes: a non-zero value means
%                              the panels do not enclose a volume, usually a missing face or a
%                              normal pointing inward. It is a DIAGNOSTIC and no force reads it;
%                              geom.graceFO's own tests assert it against the handbook geometry.
%
%   PANELS is a struct array, one element per surface element, with fields
%     name        string     label
%     kind        string     "plate" (the only kind so far; "cylinder" and
%                            "sphere" are reserved until a kernel consumes them)
%     area_m2     1x1 double area (the authority for every force)
%     normal      3x1 double outward unit normal, body frame
%     centroid_m  3x1 double application point, body frame (torques,
%                            self-shadowing); NaN if unknown
%     absVis, specVis, diffVis   visible band, must sum to 1
%     emisIR, specIR, diffIR     infrared band, must sum to 1
%     group       string     "bus" | "array" | "boom" (informational)
%   (a per-panel accommodation coefficient `alpha` was required until 2026-09-15 and read
%   by nothing; it is no longer required and, when present, ignored: HANDOVER O-23 closed
%   by removal. The accommodation is the environment's: kernel.accommodation.)
%   and OPTIONALLY (added 2026-09-07, LP-4)
%     vertices_m  3xV double corners of the plate polygon, body frame, wound
%                            counter-clockwise about the outward normal, V >= 3;
%                            used ONLY by geom.selfShadow (occlusion). Absent or
%                            empty: no ray casting possible for this model.
%     heatCapacity_J_m2K 1x1 double heat capacity per unit area of the wall
%                            behind the surface; used ONLY by geom.thermalState
%                            "transient". NaN if unknown.
%     internalShare 1x1 double share of the internal dissipation radiated through
%                            this plate; given on every plate or on none, >= 0,
%                            summing to 1 (geom.thermalState 'InternalShare').
%
%   NAME-VALUE 'Mass_kg' (added 2026-09-15, LP-C F10): the spacecraft mass a Model
%   divides by when it is given no 'Mass' of its own; NaN (the default) means
%   unknown, and a Model built on such a macro model refuses to guess.
%
%   Validation at construction is the ONLY validation: kernels trust the
%   packed arrays. Energy closure |a+s+d-1| <= 1e-6 in both bands, unit
%   normals to 1e-5 (the handbook prints sin/cos 50 deg to six digits, whose
%   norm is 1 - 7.3e-7, measured): the packed N holds the normals RENORMALISED to
%   unit length (decided 2026-09-07, HANDOVER O-2: the printed digits are a
%   truncation of a unit vector, so one renormalisation here restores the
%   authority's intent and no kernel renormalises), while panels(k).normal
%   keeps the value as given; positive areas, finite
%   centroids or all-NaN. Vertices, when given: planar to 1 mm about the
%   normal, polygon area within 5 % of area_m2 (the polygons are derived
%   from exterior dimensions, the areas are the printed authority), and
%   wound counter-clockwise (a clockwise polygon is refused, never flipped), and
%   a finite centroid_m within 5 cm of the polygon's area centroid (R1).
%
% PROPERTIES (read-only)
%   name        string
%   n           double     number of elements
%   panels      struct     as given, validated
%   A           1xN double areas
%   N           3xN double normals
%   C           3xN double centroids (NaN allowed)
%   optVis      3xN double [abs; spec; diff] visible
%   optIR       3xN double [emis; spec; diff] infrared
%   V           1xN cell   vertices (3xV each) or empty cells
%   hasVertices logical    every panel carries vertices
%   heatCapacity 1xN double J/(m^2 K), NaN where unknown
%   internalShare 1xN double internal dissipation shares, NaN where not given
%   mass_kg     1x1 double spacecraft mass, NaN if not given
%   inertia     3x3 double kg m^2 (NaN if not given)
%   centreOfMass    3x1 double m, body frame (origin = ACC proof mass for GRACE-FO)
%   proofMassOffset 3x1 double m, CG -> proof mass (CoM-coupling term, LP-7)
%   reference   string     source of the numbers
%
% NAME-VALUE ARGUMENTS
%   The properties above are read-only once built; these are how they are SET at construction.
%   The names are PascalCase where the properties they fill are camelCase, and three of them were
%   described nowhere in this help until 2026-09-18 -- the property was documented and the way to
%   fill it was not, which is a gap a reader cannot close by guessing.
%
%   Name             Size  Type    Units      Description
%   ----             ----  ----    -----      -----------
%   Name             1x1   string  -          Model label, e.g. "GRACE-C"; default ""
%   Inertia          3x3   double  kg m^2     Inertia tensor about the CoM, body frame; default NaN(3).
%                                             Required by kernel.bodyTorques' gravity-gradient term.
%   CentreOfMass     3x1   double  m          CoM in the body frame, whose ORIGIN for GRACE-FO is the
%                                             accelerometer proof mass; default zeros. Read by the
%                                             torque reduction as the point moments are taken about.
%   ProofMassOffset  3x1   double  m          CG -> proof mass, body frame; default zeros. This is the
%                                             `dr` of sensor.Electrostatic's CoM-coupling term, and it
%                                             is NOT the negative of CentreOfMass in general: the two
%                                             are independent statements about the same geometry, and
%                                             a model that carries one need not carry the other.
%   Reference        1x1   string  -          Provenance of the numbers, e.g. "JPL D-56935 Tab.4";
%                                             default "". Carried into provenance, never computed with.
%   Mass_kg          1x1   double  kg         Spacecraft mass; default NaN (a Model then needs 'Mass').
%                                             MAS1B is COLD-GAS TANK mass, not this: see io.massModel.
%
% ERRORS
%   NonGravForces:MacroModel:inputPanels     missing field / wrong shape; internalShare or Mass_kg invalid
%   NonGravForces:MacroModel:inputEnergy     a+s+d ~= 1 in a band
%   NonGravForces:MacroModel:inputNormal     non-unit normal
%   NonGravForces:MacroModel:inputVertices   vertices not 3xV, not planar, area
%                                            off by > 5 %, wound clockwise, or
%                                            centroid_m > 5 cm from the polygon centroid
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C"); sum(mm.A .* mm.N, 2)  % not zero: non-convex
%
% See also NonGravForces.geom.graceFO, NonGravForces.geom.boxWing, NonGravForces.geom.selfShadow.
%
% KEYWORDS macro model panels geometry optical properties vertices
% Created 2026-08-27 by Claude Fable 5 (Anthropic); vertices and heat capacity 2026-09-07 by Claude Fable 5.1; mass and internal shares, alpha removed 2026-09-15 by Claude Opus 5.
    properties (SetAccess = private)
        name (1,1) string = ""
        n (1,1) double = 0
        panels struct
        A double
        N double
        C double
        optVis double
        optIR double
        V cell
        hasVertices (1,1) logical = false
        heatCapacity double
        internalShare double
        mass_kg (1,1) double = NaN
        inertia (3,3) double = nan(3)
        centreOfMass (3,1) double = zeros(3,1)
        proofMassOffset (3,1) double = zeros(3,1)
        reference (1,1) string = ""
    end
    methods
        function obj = MacroModel(panels, nv)
            arguments
                panels (1,:) struct
                nv.Name (1,1) string = ""
                nv.Inertia (3,3) double = nan(3)
                nv.CentreOfMass (3,1) double = zeros(3,1)
                nv.ProofMassOffset (3,1) double = zeros(3,1)
                nv.Reference (1,1) string = ""
                nv.Mass_kg (1,1) double = NaN
            end
            req = ["name","kind","area_m2","normal","centroid_m","absVis","specVis","diffVis","emisIR","specIR","diffIR","group"];
            missing = setdiff(req, string(fieldnames(panels)));
            if ~isempty(missing)
                NonGravForces.fail("MacroModel", "inputPanels", "panel struct lacks %s -> add the fields", strjoin(missing, ", "));
            end
            n = numel(panels);
            A = [panels.area_m2]; N = reshape([panels.normal], 3, n); C = reshape([panels.centroid_m], 3, n);
            if any(~isfinite(A)) || any(A <= 0)
                NonGravForces.fail("MacroModel", "inputPanels", "areas must be finite and positive -> check panel %d", find(~(A > 0), 1));
            end
            nn = sqrt(sum(N.^2, 1));
            if any(abs(nn - 1) > 1e-5)
                NonGravForces.fail("MacroModel", "inputNormal", "normal of panel %d has norm %.3e -> normalise it", find(abs(nn-1) > 1e-5, 1), nn(find(abs(nn-1)>1e-5,1)));
            end
            N = N ./ nn;                                   % packed normals are exactly unit (O-2); panels keep the printed value
            oV = [[panels.absVis]; [panels.specVis]; [panels.diffVis]];
            oI = [[panels.emisIR]; [panels.specIR]; [panels.diffIR]];
            if any(abs(sum(oV,1) - 1) > 1e-6)
                k = find(abs(sum(oV,1)-1) > 1e-6, 1);
                NonGravForces.fail("MacroModel", "inputEnergy", "visible a+s+d = %.4f on panel '%s' -> the transcription is wrong", sum(oV(:,k)), panels(k).name);
            end
            if any(abs(sum(oI,1) - 1) > 1e-6)
                k = find(abs(sum(oI,1)-1) > 1e-6, 1);
                NonGravForces.fail("MacroModel", "inputEnergy", "infrared e+s+d = %.4f on panel '%s' -> the transcription is wrong", sum(oI(:,k)), panels(k).name);
            end
            % optional vertices: planar, counter-clockwise, area within 5 % of the printed area
            V = cell(1, n);
            if isfield(panels, 'vertices_m')
                for k = 1:n
                    vk = panels(k).vertices_m;
                    if isempty(vk), continue; end
                    if ~(ismatrix(vk) && size(vk, 1) == 3 && size(vk, 2) >= 3 && all(isfinite(vk(:))))
                        NonGravForces.fail("MacroModel", "inputVertices", "vertices_m of panel '%s' must be a finite 3xV array with V >= 3 -> one corner per column", panels(k).name);
                    end
                    ctr = mean(vk, 2);
                    off = N(:, k)' * (vk - ctr);
                    if any(abs(off) > 1e-3)
                        NonGravForces.fail("MacroModel", "inputVertices", "vertices_m of panel '%s' are %.1f mm off the plane of its normal -> fix the polygon or the normal", panels(k).name, 1e3 * max(abs(off)));
                    end
                    e1 = vk(:, 2) - vk(:, 1); e1 = e1 / norm(e1); e2 = cross(N(:, k), e1);
                    x = e1' * (vk - vk(:, 1)); y = e2' * (vk - vk(:, 1));
                    signedArea = 0.5 * sum(x .* circshift(y, -1) - circshift(x, -1) .* y);
                    if signedArea <= 0
                        NonGravForces.fail("MacroModel", "inputVertices", "vertices_m of panel '%s' are wound clockwise about the outward normal -> reverse the column order", panels(k).name);
                    end
                    if all(isfinite(C(:, k)))
                        cPoly = localCentroid(vk);
                        if norm(cPoly - C(:, k)) > 0.05
                            NonGravForces.fail("MacroModel", "inputVertices", "centroid_m of panel '%s' is %.2f m from the centroid of its polygon -> one of the two is wrong (rule R1: derive the centroid from the polygon)", panels(k).name, norm(cPoly - C(:, k)));
                        end
                    end
                    if abs(signedArea - A(k)) / A(k) > 0.05
                        NonGravForces.fail("MacroModel", "inputVertices", "polygon area %.4f m^2 of panel '%s' differs from area_m2 %.4f by more than 5 %% -> the polygon does not describe this plate", signedArea, panels(k).name, A(k));
                    end
                    V{k} = vk;
                end
            end
            hc = nan(1, n);
            if isfield(panels, 'heatCapacity_J_m2K')
                for k = 1:n
                    if ~isempty(panels(k).heatCapacity_J_m2K), hc(k) = panels(k).heatCapacity_J_m2K; end
                end
                if any(hc <= 0)
                    NonGravForces.fail("MacroModel", "inputPanels", "heatCapacity_J_m2K must be positive or NaN -> check panel %d", find(hc <= 0, 1));
                end
            end
            share = nan(1, n);
            if isfield(panels, 'internalShare')
                for k = 1:n
                    if ~isempty(panels(k).internalShare), share(k) = panels(k).internalShare; end
                end
                if any(isfinite(share)) && ~(all(isfinite(share)) && all(share >= 0) && abs(sum(share) - 1) <= 1e-9)
                    NonGravForces.fail("MacroModel", "inputPanels", "internalShare must be given on every panel or on none, non-negative, summing to 1 -> got a sum of %.6g", sum(share(isfinite(share))));
                end
            end
            if ~(isnan(nv.Mass_kg) || (isfinite(nv.Mass_kg) && nv.Mass_kg > 0))
                NonGravForces.fail("MacroModel", "inputPanels", "'Mass_kg' must be a positive finite mass or NaN (unknown) -> got %g", nv.Mass_kg);
            end
            obj.name = nv.Name; obj.n = n; obj.panels = panels; obj.A = A; obj.N = N; obj.C = C;
            obj.optVis = oV; obj.optIR = oI; obj.V = V; obj.hasVertices = all(~cellfun(@isempty, V));
            obj.internalShare = share; obj.mass_kg = nv.Mass_kg;
            obj.heatCapacity = hc;
            obj.inertia = nv.Inertia; obj.centreOfMass = nv.CentreOfMass; obj.proofMassOffset = nv.ProofMassOffset; obj.reference = nv.Reference;
        end
        function s = silhouetteDefect(obj)
            %SILHOUETTEDEFECT sum_k A_k n_k; zero for a closed convex body.
            s = sum(obj.A .* obj.N, 2);
        end
    end
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
