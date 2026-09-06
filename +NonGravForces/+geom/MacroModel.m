classdef MacroModel
%MACROMODEL Panel ("macro") geometry of a spacecraft, frozen and validated once.
%
%   mm = NonGravForces.geom.MacroModel(panels, 'Name', "GRACE-C", 'Inertia', I)
%
%   PANELS is a struct array, one element per surface element, with fields
%     name        string     label
%     kind        string     "plate" (the only kind so far; "cylinder" and
%                            "sphere" are reserved, see LP-4)
%     area_m2     1x1 double area
%     normal      3x1 double outward unit normal, body frame
%     centroid_m  3x1 double application point, body frame (torques,
%                            self-shadowing); NaN if unknown
%     absVis, specVis, diffVis   visible band, must sum to 1
%     emisIR, specIR, diffIR     infrared band, must sum to 1
%     alpha       1x1 double energy accommodation coefficient default (0..1)
%     group       string     "bus" | "array" | "boom" (informational)
%
%   Validation at construction is the ONLY validation: kernels trust the
%   packed arrays. Energy closure |a+s+d-1| <= 1e-6 in both bands, unit
%   normals to 1e-5 (the handbook prints sin/cos 50 deg to six digits, whose
%   squares sum to 1 - 5e-7; kernels renormalise nothing, the 5e-7 is
%   carried as a documented transcription defect), positive areas, finite centroids or all-NaN.
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
%   alpha       1xN double accommodation defaults
%   inertia     3x3 double kg m^2 (NaN if not given)
%   centreOfMass    3x1 double m, body frame (origin = ACC proof mass for GRACE-FO)
%   proofMassOffset 3x1 double m, CG -> proof mass (CoM-coupling term, LP-7)
%   reference   string     source of the numbers
%
% ERRORS
%   NonGravForces:MacroModel:inputPanels     missing field / wrong shape
%   NonGravForces:MacroModel:inputEnergy     a+s+d ~= 1 in a band
%   NonGravForces:MacroModel:inputNormal     non-unit normal
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C"); sum(mm.A .* mm.N, 2)  % not zero: non-convex
%
% See also NonGravForces.geom.graceFO, NonGravForces.geom.boxWing.
%
% KEYWORDS macro model panels geometry optical properties
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    properties (SetAccess = private)
        name (1,1) string = ""
        n (1,1) double = 0
        panels struct
        A double
        N double
        C double
        optVis double
        optIR double
        alpha double
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
            end
            req = ["name","kind","area_m2","normal","centroid_m","absVis","specVis","diffVis","emisIR","specIR","diffIR","alpha","group"];
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
            al = [panels.alpha];
            if any(al < 0 | al > 1)
                NonGravForces.fail("MacroModel", "inputPanels", "alpha must be in [0,1] -> check panel %d", find(al<0|al>1,1));
            end
            obj.name = nv.Name; obj.n = n; obj.panels = panels; obj.A = A; obj.N = N; obj.C = C;
            obj.optVis = oV; obj.optIR = oI; obj.alpha = al; obj.inertia = nv.Inertia; obj.centreOfMass = nv.CentreOfMass; obj.proofMassOffset = nv.ProofMassOffset; obj.reference = nv.Reference;
        end
        function s = silhouetteDefect(obj)
            %SILHOUETTEDEFECT sum_k A_k n_k; zero for a closed convex body.
            s = sum(obj.A .* obj.N, 2);
        end
    end
end
