function [F, tau, info] = assemble(Fk, macro, options)
%ASSEMBLE Reduce per-panel forces to the total force and the torque about the centre of mass.
%
%   [F, tau] = NonGravForces.kernel.assemble(Fk, macro)
%   [F, tau, info] = NonGravForces.kernel.assemble(Fk, macro, 'CentreOfMass', r_cm)
%
%   F = sum_k F_k and tau = sum_k (c_k - r_cm) x F_k, with c_k the plate
%   centroids of the macro model (body frame). Force and torque are two
%   reductions over ONE per-panel array (principle P2), so the torque path can
%   never lose a term the force path has (the predecessor's separate torque
%   implementation dropped the vis/IR split, 8 %). Several effects are summed
%   by concatenating their per-panel arrays along the panel dimension, or by
%   adding the arrays before the call; the result is the same to rounding.
%
% INPUT ARGUMENTS
%   Name     Size    Type    Units  Description
%   ----     ----    ----    -----  -----------
%   Fk       3xNxK   double  N      Per-panel forces, body frame (kernel.dragPanelForce,
%                                   radiationPanelForce, thermalRecoil, or their sum)
%   macro    1x1     MacroModel     Panel model whose centroids C (3xN) are the lever arms
%
% NAME-VALUE ARGUMENTS
%   Name          Size  Type    Units  Description
%   ----          ----  ----    -----  -----------
%   CentreOfMass  3x1   double  m      Reference point of the torque, default macro.centreOfMass
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   F        3xK   double  N      Total force, body frame
%   tau      3xK   double  N m    Torque about CentreOfMass, body frame
%   info     1x1   struct  -      .lever_m (3xN, c_k - r_cm), .N, .K, .provenance
%
% ACCURACY
%   Sums in a fixed panel order (column-wise, not through a BLAS product), so a
%   batch of K epochs is bitwise the same as K single calls (vectorisation
%   category). Rounding 1e-16 relative to the largest panel force.
%
% ERRORS
%   NonGravForces:assemble:inputMacro   macro is not a MacroModel
%   NonGravForces:assemble:inputShape   Fk not 3xNxK with N = macro.n, or CentreOfMass not 3x1
%   NonGravForces:assemble:numericNaN   NaN/Inf in Fk or CentreOfMass, or non-finite panel centroids
%
% REFERENCES
%   HANDOVER D.2 principle P2 and Part C defect D8; theory §7.
%
% EXAMPLE
%   Fk = NonGravForces.kernel.radiationPanelForce(mm, src) + NonGravForces.kernel.thermalRecoil(mm, T);
%   [F, tau] = NonGravForces.kernel.assemble(Fk, mm);
%
% See also NonGravForces.kernel.dragPanelForce, NonGravForces.kernel.radiationPanelForce, NonGravForces.kernel.bodyTorques.
%
% KEYWORDS assemble force torque reduction per-panel
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    Fk double
    macro
    options.CentreOfMass double = []
end
unit = "assemble";
if ~isa(macro, 'NonGravForces.geom.MacroModel')
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel -> the same model the per-panel forces were computed with");
end
N = macro.n;
if ndims(Fk) > 3 || size(Fk, 1) ~= 3 || size(Fk, 2) ~= N
    NonGravForces.fail(unit, "inputShape", "Fk must be 3xNxK with N = %d -> got %s", N, mat2str(size(Fk)));
end
if isempty(options.CentreOfMass)
    rcm = macro.centreOfMass;
else
    rcm = options.CentreOfMass;
    if ~isequal(size(rcm), [3 1])
        NonGravForces.fail(unit, "inputShape", "CentreOfMass must be 3x1 -> got %s", mat2str(size(rcm)));
    end
end
if any(~isfinite(Fk(:))) || any(~isfinite(rcm))
    NonGravForces.fail(unit, "numericNaN", "Fk or CentreOfMass contain NaN/Inf -> screen the per-panel forces");
end
K = size(Fk, 3);
lever = macro.C - rcm;                              % 3xN
if any(~isfinite(lever(:)))
    NonGravForces.fail(unit, "numericNaN", "the panel centroids are not all finite, so no torque can be formed -> build the MacroModel with real centroid_m (MacroModel permits all-NaN centroids for force-only models)");
end
F = zeros(3, K); tau = zeros(3, K);
for k = 1:N
    f = reshape(Fk(:, k, :), 3, K);
    F = F + f;
    r = lever(:, k);
    tau = tau + [r(2) * f(3, :) - r(3) * f(2, :); r(3) * f(1, :) - r(1) * f(3, :); r(1) * f(2, :) - r(2) * f(1, :)];
end
info.lever_m = lever;
info.N = N; info.K = K;
info.provenance = NonGravForces.provenance("kernel.assemble", 'Options', struct('CentreOfMass', rcm));
end
