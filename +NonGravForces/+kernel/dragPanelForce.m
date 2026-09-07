function [F, info] = dragPanelForce(macro, rho, species, T_K, vrel_body, options)
%DRAGPANELFORCE Per-panel free-molecular gas force in the body frame, species by species.
%
%   [F, info] = NonGravForces.kernel.dragPanelForce(macro, rho, species, T_K, vrel, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300)
%   [F, info] = NonGravForces.kernel.dragPanelForce(macro, rho, [], [], vrel, 'Model', "constantCd", 'Cd', 2.2)
%   [F, info] = NonGravForces.kernel.dragPanelForce(___, 'Shadow', f)
%
%   For every plate k and epoch, with v the satellite velocity relative to the
%   gas (body frame), gamma_k = n_k . v/|v| and the tangential stream direction
%   t_k = (-v/|v| + gamma_k n_k)/l_k,
%
%       F_k = (1/2) |v|^2 A_k f_k  sum_j rho_j [ Cp_jk (-n_k) + Ctau_jk t_k ],
%
%   where the sum runs over the MSIS species j with mass density rho_j = n_j m_j
%   and its own speed ratio s_j = |v| / sqrt(2 k T / m_j) (Doornbos 2011 eq. 3.57:
%   the mixture coefficient is the mass-weighted mean of the per-species
%   coefficients; treating the gas as air with R = 287 J/kg/K biases s by
%   sqrt(29/16) = 1.35 at 500 km, defect D1b of the predecessor). The
%   coefficients come from gsiCoefficients; f_k in [0, 1] is the self-shadow
%   fraction of the plate for the stream direction (geom.selfShadow), 1 when
%   'Shadow' is omitted (no self-shadowing, a documented choice of the caller).
%   The model constantCd uses rho directly and needs neither species nor T.
%   Forces are per plate (principle P2): assemble sums them and forms torques.
%
% INPUT ARGUMENTS
%   Name        Size  Type    Units   Description
%   ----        ----  ----    -----   -----------
%   macro       1x1   MacroModel      Panel model (areas A, normals N in body frame)
%   rho         1xK   double  kg/m^3  Total mass density (env.density); used by constantCd
%                                     and reported for the species-sum check
%   species     8xK   double  1/m^3   Number densities in const.speciesName order (env.density);
%                                     may be [] for constantCd
%   T_K         1xK   double  K       Gas kinetic temperature at the satellite (density's T(2,:));
%                                     may be [] for constantCd
%   vrel_body   3xK   double  m/s     Satellite velocity relative to the gas, body frame, |v| > 0
%
% NAME-VALUE ARGUMENTS
%   Name              Size      Type     Units  Description
%   ----              ----      ----     -----  -----------
%   Model             1x1       string   -      "sentman" | "maxwell" | "constantCd" (required; no default)
%   Alpha             1x1|1xK   double   -      Accommodation coefficient (kernel.accommodation), required by sentman/maxwell
%   Tw_K              1x1|1xK   double   K      Wall temperature, required by sentman/maxwell
%   Cd                1x1       double   -      Drag coefficient, required by constantCd
%   SpecularFraction  1x1       double   -      Required by maxwell
%   Shadow            NxK       double   -      Self-shadow fraction per plate and epoch, default all 1
%
% OUTPUT ARGUMENTS
%   Name     Size    Type    Units  Description
%   ----     ----    ----    -----  -----------
%   F        3xNxK   double  N      Force on each plate, body frame
%   info     1x1     struct  -      .gamma (NxK), .speed_m_s (1xK), .speedRatio (8xK, NaN for
%                                   constantCd), .Cp .Ctau (NxK, mass-weighted mixture
%                                   coefficients per unit rho), .rhoSpecies_kg_m3 (1xK, sum of
%                                   n_j m_j), .model, .provenance
%
% ACCURACY
%   Exact assembly of the gsiCoefficients closed forms; rounding 1e-15. The
%   species sum does NOT reproduce the density backend's total exactly: measured
%   3.7e-4 to 5.7e-4 relative against msis00 (2026-09-07, three epochs at 480-500 km),
%   because NRLMSISE-00 forms its total from its own molecular masses and this unit
%   re-forms it from const.speciesMass. The forces use the species sum, so the
%   difference is a statement about the reported total, not an error in the force.
%   (This block claimed 1e-6 until 2026-09-07; the number had been asserted, not
%   measured - rule R3.) The drag force error is
%   dominated by the density (10-30 %, HANDOVER E) and by the gas-surface
%   model (5-10 %), not by this assembly. A plate exactly edge-on (l = 0 or
%   gamma = 0) gets its exact limit, no NaN.
%
% ERRORS
%   NonGravForces:dragPanelForce:inputMacro    macro is not a NonGravForces.geom.MacroModel
%   NonGravForces:dragPanelForce:inputModel    'Model' missing or not in the set
%   NonGravForces:dragPanelForce:inputShape    rho, species, T_K, vrel_body or Shadow with the wrong shape
%   NonGravForces:dragPanelForce:inputRange    |v| = 0, rho < 0, T_K <= 0, or Shadow outside [0, 1]
%   NonGravForces:dragPanelForce:numericNaN    NaN/Inf in an input
%   NonGravForces:gsiCoefficients:*            propagated from the coefficient model (e.g. inputCd, inputSpecular)
%
% REFERENCES
%   Sentman (1961); Doornbos E. (2011) Thermospheric density and wind
%   determination from satellite dynamics, Ch. 3; Moe K., Moe M.M. (2005)
%   Planet. Space Sci. 53:793-801; HANDOVER C.1 defects D1, D1b.
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C");
%   [rho, sp, T] = NonGravForces.env.density("msis00", ep, lla, sw);
%   F = NonGravForces.kernel.dragPanelForce(mm, rho, sp, T(2,:), vrel, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
%   a = squeeze(sum(F, 2)) / 601.214;        % 3xK acceleration, ~5e-8 m/s^2 along -x at solar minimum
%
% See also NonGravForces.kernel.gsiCoefficients, NonGravForces.kernel.assemble, NonGravForces.env.density, NonGravForces.geom.selfShadow.
%
% KEYWORDS drag per-panel force free molecular species Sentman
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    macro
    rho double
    species double
    T_K double
    vrel_body double
    options.Model string = string.empty
    options.Alpha double = []
    options.Tw_K double = []
    options.Cd double = []
    options.SpecularFraction double = []
    options.Shadow double = []
end
unit = "dragPanelForce";
if ~isa(macro, 'NonGravForces.geom.MacroModel')
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel -> build one with geom.graceFO or geom.boxWing");
end
if isempty(options.Model) || ~isscalar(options.Model) || ~any(options.Model == ["sentman","maxwell","constantCd"])
    NonGravForces.fail(unit, "inputModel", "'Model' must be sentman, maxwell or constantCd -> state the gas-surface model explicitly");
end
model = options.Model;
K = size(vrel_body, 2);
if ~(ismatrix(vrel_body) && size(vrel_body, 1) == 3 && isequal(size(rho), [1 K]))
    NonGravForces.fail(unit, "inputShape", "vrel_body must be 3xK and rho 1xK -> got %s and %s", mat2str(size(vrel_body)), mat2str(size(rho)));
end
if model ~= "constantCd" && ~(isequal(size(species), [8 K]) && isequal(size(T_K), [1 K]))
    NonGravForces.fail(unit, "inputShape", "species must be 8xK and T_K 1xK for %s -> got %s and %s (env.density returns them)", model, mat2str(size(species)), mat2str(size(T_K)));
end
N = macro.n;
if isempty(options.Shadow)
    fShadow = ones(N, K);
else
    fShadow = options.Shadow;
    if ~isequal(size(fShadow), [N K])
        NonGravForces.fail(unit, "inputShape", "Shadow must be NxK = %dx%d -> got %s (geom.selfShadow for the stream direction)", N, K, mat2str(size(fShadow)));
    end
    if any(fShadow(:) < 0 | fShadow(:) > 1)
        NonGravForces.fail(unit, "inputRange", "Shadow fractions must lie in [0, 1] -> check the selfShadow output");
    end
end
if any(~isfinite(vrel_body(:))) || any(~isfinite(rho)) || any(~isfinite(species(:))) || any(~isfinite(T_K(:))) || any(~isfinite(fShadow(:)))
    NonGravForces.fail(unit, "numericNaN", "an input contains NaN/Inf -> screen the density, temperature and velocity");
end
U = vecnorm(vrel_body);
if any(U <= 0) || any(rho < 0) || (model ~= "constantCd" && any(T_K <= 0))
    NonGravForces.fail(unit, "inputRange", "|vrel| must be > 0, rho >= 0 and T_K > 0 -> check the inputs");
end
v = vrel_body ./ U;                                 % 3xK, direction of the satellite through the gas
Nrm = macro.N;                                      % 3xN outward normals
A = macro.A(:);                                     % Nx1
c = NonGravForces.const();
% incidence per plate: gamma = n.v, l = |n x v| (atan2 form, well conditioned at 0 and 180 deg)
gamma = Nrm' * v;                                   % NxK
l = zeros(N, K); tdir = zeros(3, N, K);
for k = 1:N
    nk = Nrm(:, k);
    cr = [nk(2) * v(3, :) - nk(3) * v(2, :); nk(3) * v(1, :) - nk(1) * v(3, :); nk(1) * v(2, :) - nk(2) * v(1, :)];
    l(k, :) = vecnorm(cr);
    t = (-v + gamma(k, :) .* nk);                   % tangential part of the gas stream g = -v (unnormalised, norm l)
    ok = l(k, :) > 0;
    t(:, ok) = t(:, ok) ./ l(k, ok);
    t(:, ~ok) = 0;                                  % edge-on or head-on: no tangential direction, Ctau -> 0 there
    tdir(:, k, :) = reshape(t, 3, 1, K);
end
theta = atan2(l, gamma);
F = zeros(3, N, K);
switch model
    case "constantCd"
        [Cp, Ctau, ginfo] = NonGravForces.kernel.gsiCoefficients("constantCd", [], theta, [], [], [], 'Cd', options.Cd);
        Cpm = rho .* Cp; Ctm = rho .* Ctau;         % NxK, mass-weighted (single "species")
        sRatio = nan(8, K); rhoSp = rho;
    otherwise
        if isempty(options.Alpha) || isempty(options.Tw_K)
            NonGravForces.fail(unit, "inputModel", "%s needs 'Alpha' and 'Tw_K' -> pass the accommodation coefficient (kernel.accommodation) and the wall temperature", model);
        end
        alpha = NonGravForces.kernel.accommodation("fixed", species, T_K, 'Alpha', options.Alpha);
        Tw = options.Tw_K + zeros(1, K);
        if any(Tw <= 0)
            NonGravForces.fail(unit, "inputRange", "Tw_K must be > 0 -> got %g", min(Tw));
        end
        mass = c.speciesMass / c.NA;                % kg per molecule, 1x8
        sRatio = U ./ sqrt(2 * c.kB * T_K ./ mass(:));   % 8xK
        Cpm = zeros(N, K); Ctm = zeros(N, K); rhoSp = zeros(1, K);
        extra = {};
        if model == "maxwell", extra = {'SpecularFraction', options.SpecularFraction}; end
        for j = 1:8
            rhoj = species(j, :) * mass(j);         % 1xK
            rhoSp = rhoSp + rhoj;
            if ~any(rhoj > 0), continue; end
            [Cpj, Ctj] = NonGravForces.kernel.gsiCoefficients(model, sRatio(j, :), theta, alpha, Tw, T_K, extra{:});
            Cpm = Cpm + rhoj .* Cpj;
            Ctm = Ctm + rhoj .* Ctj;
        end
        ginfo = struct('model', model);
end
scale = 0.5 * U.^2 .* A .* fShadow;                 % NxK
fn = -scale .* Cpm;                                 % along -n
ft = scale .* Ctm;                                  % along t
for k = 1:N
    F(:, k, :) = reshape(Nrm(:, k) .* fn(k, :), 3, 1, K) + tdir(:, k, :) .* reshape(ft(k, :), 1, 1, K);
end
info.model = ginfo.model;
info.gamma = gamma;
info.theta_rad = theta;
info.speed_m_s = U;
info.speedRatio = sRatio;
info.rhoSpecies_kg_m3 = rhoSp;
info.Cp = Cpm ./ max(rhoSp, realmin);
info.Ctau = Ctm ./ max(rhoSp, realmin);
info.shadow = fShadow;
info.provenance = NonGravForces.provenance("kernel.dragPanelForce", 'Options', struct('Model', model, 'Alpha', options.Alpha, 'Tw_K', options.Tw_K, 'Cd', options.Cd, 'SpecularFraction', options.SpecularFraction, 'shadowGiven', ~isempty(options.Shadow)));
end
