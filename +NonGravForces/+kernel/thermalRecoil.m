function [F, info] = thermalRecoil(macro, T_K)
%THERMALRECOIL Per-panel recoil force of Lambertian thermal emission in the body frame.
%
%   [F, info] = NonGravForces.kernel.thermalRecoil(macro, T_K)
%
%   A plate of area A_k, infrared emissivity eps_k and temperature T_k emits
%   eps_k sigma T_k^4 A_k watts as a Lambertian (cosine) source and recoils
%   along its inward normal:
%
%       F_k = - (2/3) eps_k sigma T_k^4 A_k / c  n_k          (Wöske 2019 eq. 5.14)
%
%   The factor 2/3 is the hemispherical integral of the cosine law (only the
%   normal components add up). The temperatures come from geom.thermalState
%   (steady state or transient); this unit does no energy balance itself, so a
%   caller can also feed measured or otherwise modelled temperatures.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   macro    1x1   MacroModel     Panel model (A, N, optIR(1,:) = emissivity)
%   T_K      NxK   double  K      Surface temperature per plate and epoch, >= 0 (an unlit steady-state
%                                  plate without internal power sits at 0 K and emits nothing)
%
% OUTPUT ARGUMENTS
%   Name     Size    Type    Units  Description
%   ----     ----    ----    -----  -----------
%   F        3xNxK   double  N      Recoil force on each plate, body frame
%   info     1x1     struct  -      .emitted_W_m2 (NxK, eps sigma T^4), .emittedPower_W (1xK, total),
%                                   .provenance
%
% ACCURACY
%   Exact for a Lambertian emitter at a uniform plate temperature; rounding
%   1e-16. The model error is the temperature (geom.thermalState: steady
%   state ignores the wall heat capacity, the transient model the conduction
%   between plates) and the Lambertian assumption; for GRACE, TRP is the
%   second-largest radiation term (Wöske et al. 2019; Vielberg & Kusche 2020;
%   Siemes et al. 2023: 6.4 nm/s^2 cross-track RMS), so the temperature error
%   dominates, not this formula. Measured here (steadyState, 2026-09-17) on GRACE-FO 2020-01-01: 5.263
%   nm/s^2 cross-track RMS, between Hladczuk et al. (2024)'s 4.22 for GRACE-FO
%   and Siemes' 6.4 for GRACE-A.
%
%   THE LITERATURE DISAGREES ABOUT WHETHER MODELLING THIS INSTANTANEOUSLY HELPS,
%   and this block cited only the agreeing side until 2026-09-17.
%     Hladczuk et al. (2024), GRACE-FO, Table 3: adding instantaneous reradiation
%       to a panel model takes the cross-track residual from 4.22 to 3.04 nm/s^2,
%       and thermal inertia takes it further to 1.70. Modelling it helps.
%     Vielberg & Kusche (2020), GRACE-A, one year from 2010: instantaneous thermal
%       reradiation fits "more than 39% worse" to calibrated cross-track
%       accelerations than the standard model, and the conclusions state it "did
%       not lead to improvements". Modelling it hurts.
%   The standard and precise presets follow the first (standard with Hladczuk's
%   further step, the transient walls of geom.thermalState, since 2026-10-31, both;
%   fast sets thermal.model "none"), for reasons that
%   are about the evidence and not about preference: it is the right satellite, it is four years later, and its
%   whole Section 3 is devoted to repairing the cross-track calibration that
%   Vielberg & Kusche name as a possible cause of their own result (they also
%   name neglected horizontal winds in the drag model, which this toolbox models
%   through env.wind). A user who wants the other position sets thermal.model
%   "none". It is a genuine open question in the field, not a settled one.
%
% ERRORS
%   NonGravForces:thermalRecoil:inputMacro   macro is not a MacroModel
%   NonGravForces:thermalRecoil:inputShape   T_K is not NxK
%   NonGravForces:thermalRecoil:inputRange   a negative temperature
%   NonGravForces:thermalRecoil:numericNaN   NaN/Inf in T_K
%
% REFERENCES
%   Wöske F. (2019) PhD thesis, Univ. Bremen, Sect. 5.1.4 eq. 5.14; Rievers B.,
%   Lämmerzahl C. (2011) Ann. Phys. 523:439; Vielberg K., Kusche J. (2020)
%   J. Geodesy 94:71; Siemes C. et al. (2023) Adv. Space Res.
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C");
%   T = NonGravForces.geom.thermalState(mm, src);          % NxK steady state
%   F = NonGravForces.kernel.thermalRecoil(mm, T);          % (2/3) of the emitted power over c
%
% See also NonGravForces.geom.thermalState, NonGravForces.kernel.radiationPanelForce, NonGravForces.kernel.assemble.
%
% KEYWORDS thermal radiation pressure recoil emission Lambert
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    macro
    T_K double
end
unit = "thermalRecoil";
if ~isa(macro, 'NonGravForces.geom.MacroModel')
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel -> build one with geom.graceFO or geom.boxWing");
end
N = macro.n;
if ~(ismatrix(T_K) && size(T_K, 1) == N)
    NonGravForces.fail(unit, "inputShape", "T_K must be NxK with N = %d plates -> got %s (geom.thermalState returns it)", N, mat2str(size(T_K)));
end
if any(~isfinite(T_K(:)))
    NonGravForces.fail(unit, "numericNaN", "T_K contains NaN/Inf -> check thermalState's input fluxes");
end
if any(T_K(:) < 0)
    NonGravForces.fail(unit, "inputRange", "T_K must be >= 0 K -> got %g", min(T_K(:)));
end
K = size(T_K, 2);
c = NonGravForces.const();
cLight = referenceSystem.iersConstants().speedOfLight;
eps = macro.optIR(1, :)';                           % Nx1 emissivity
emitted = eps .* c.sigmaSB .* T_K.^4;               % NxK  W/m^2
mag = (2 / 3) * emitted .* macro.A(:) / cLight;     % NxK  N
F = zeros(3, N, K);
for k = 1:N
    F(:, k, :) = reshape(-macro.N(:, k) .* mag(k, :), 3, 1, K);
end
info.emitted_W_m2 = emitted;
info.emittedPower_W = sum(emitted .* macro.A(:), 1);
info.provenance = NonGravForces.provenance("kernel.thermalRecoil");
end
