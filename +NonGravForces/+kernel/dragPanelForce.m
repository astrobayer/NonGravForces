function [F, info] = dragPanelForce(macro, rho, species, T_K, vrel_body, options)
%DRAGPANELFORCE Per-panel free-molecular gas force in the body frame, species by species.
%
%   [F, info] = NonGravForces.kernel.dragPanelForce(macro, rho, species, T_K, vrel, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300)
%   [F, info] = NonGravForces.kernel.dragPanelForce(macro, rho, [], [], vrel, 'Model', "constantCd", 'Cd', 2.2)
%   [F, info] = NonGravForces.kernel.dragPanelForce(___, 'Shadow', geom.flowShadow(...))
%   [F, info] = NonGravForces.kernel.dragPanelForce(___, 'Alpha', kernel.accommodation("walker2014", ...))
%
%   For every plate k and epoch, with v the satellite velocity relative to the
%   gas (body frame), gamma_k = n_k . v/|v| and the tangential stream direction
%   t_k = (-v/|v| + gamma_k n_k)/l_k,
%
%       F_k = (1/2) |v|^2 A_k  sum_j rho_j [ (fP Cp_inc + fN Cp_re)_jk (-n_k) + (fT Ctau)_jk t_k ],
%
%   where the sum runs over the MSIS species j with mass density rho_j = n_j m_j
%   and its own speed ratio s_j = |v| / sqrt(2 k T / m_j) (Doornbos 2011 eq. 3.57:
%   the mixture coefficient is the mass-weighted mean of the per-species
%   coefficients; treating the gas as air with R = 287 J/kg/K biases s by
%   sqrt(29/16) = 1.35 at 500 km, defect D1b of the predecessor). Cp_inc and
%   Cp_re are gsiCoefficients' incident and re-emitted normal coefficients.
%
%   SELF-SHADOWING ('Shadow'). A geom.flowShadow struct supplies, per plate,
%   epoch and species, the occluded fraction of the incident number flux (fN,
%   scaling the re-emitted momentum), of the normal momentum (fP) and of the
%   tangential momentum (fT); constantCd reads its parallel-beam field .beam. A
%   plain NxK array is one fraction applied to every term of every species (the
%   hyperthermal approximation). Omitted, there is no self-shadowing. Until
%   2026-09-15 Model passed a radiation-beam fraction that was 0 on every plate
%   the stream grazes, which deleted their Sentman shear (HANDOVER LP-C, F1).
%
%   ACCOMMODATION ('Alpha'). A scalar or 1xK coefficient, or a
%   kernel.accommodation struct. A struct with .coverage describes a surface of
%   two states (walker2014): the re-emitted momentum is
%   theta Cp_re(alpha_a) + (1 - theta) Cp_re(alpha_s), which is the drag-coefficient
%   mixture of that model exactly, because alpha enters nothing else.
%
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
%   Name              Size      Type           Units  Description
%   ----              ----      ----           -----  -----------
%   Model             1x1       string         -      "sentman" | "maxwell" | "constantCd" (required; no default)
%   Alpha             1x1|1xK   double|struct  -      Accommodation coefficient, or a kernel.accommodation struct;
%                                                     required by sentman/maxwell
%   Tw_K              1x1|1xK   double         K      Wall temperature, required by sentman/maxwell
%   Cd                1x1       double         -      Drag coefficient, required by constantCd
%   SpecularFraction  1x1       double         -      Required by maxwell
%   Shadow            NxK       double|struct  -      Self-shadowing: a geom.flowShadow struct, or one fraction per
%                                                     plate and epoch; default none
%
% OUTPUT ARGUMENTS
%   Name     Size    Type    Units  Description
%   ----     ----    ----    -----  -----------
%   F        3xNxK   double  N      Force on each plate, body frame
%   info     1x1     struct  -      .gamma (NxK), .theta_rad (NxK), .speed_m_s (1xK), .speedRatio
%                                   (8xK, NaN for constantCd), .Cp .Ctau (NxK, mass-weighted mixture
%                                   coefficients per unit rho, shadowing included), .rhoSpecies_kg_m3
%                                   (1xK, sum of n_j m_j), .model, .shadow ("none" | "fraction" |
%                                   "flow"), .provenance
%
% ACCURACY
%   Exact assembly of the gsiCoefficients closed forms; rounding 1e-15. The
%   species sum does NOT reproduce the density backend's total exactly: measured
%   3.7e-4 to 5.7e-4 relative against msis00 (2026-09-07, three epochs at 480-500 km),
%   because NRLMSISE-00 forms its total from its own molecular masses and this unit
%   re-forms it from const.speciesMass. The forces use the species sum, so the
%   difference is a statement about the reported total, not an error in the force.
%   The drag force error is dominated by the density (10-30 %, HANDOVER E) and by
%   the gas-surface model (5-10 %), not by this assembly. A plate exactly edge-on
%   (l = 0 or gamma = 0) gets its exact limit, no NaN. The self-shadow ratios
%   carry their own accuracy statement (geom.flowShadow).
%
% ERRORS
%   NonGravForces:dragPanelForce:inputMacro    macro is not a NonGravForces.geom.MacroModel
%   NonGravForces:dragPanelForce:inputModel    'Model' missing or not in the set, or its required options absent
%   NonGravForces:dragPanelForce:inputShape    rho, species, T_K, vrel_body, Alpha or Shadow with the wrong shape
%   NonGravForces:dragPanelForce:inputRange    |v| = 0, rho < 0, T_K <= 0, Tw <= 0, a fraction or a number/normal ratio outside [0, 1]
%                                              (the tangential ratio is signed and need only be finite)
%   NonGravForces:dragPanelForce:numericNaN    NaN/Inf in an input
%   NonGravForces:gsiCoefficients:*            propagated from the coefficient model (e.g. inputCd, inputSpecular)
%   NonGravForces:accommodation:*              propagated when a numeric Alpha is validated
%
% REFERENCES
%   Sentman (1961); Doornbos E. (2011) Thermospheric density and wind
%   determination from satellite dynamics, Ch. 3; Moe K., Moe M.M. (2005)
%   Planet. Space Sci. 53:793-801; Walker, Mehta & Koller (2014) J. Spacecraft
%   Rockets 51(5); HANDOVER C.1 defects D1, D1b; HANDOVER LP-C WP-C1, WP-C7.
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C");
%   [rho, sp, T] = NonGravForces.env.density("msis00", ep, lla, sw);
%   F = NonGravForces.kernel.dragPanelForce(mm, rho, sp, T(2,:), vrel, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
%   a = squeeze(sum(F, 2)) / 601.214;        % 3xK acceleration, ~5e-8 m/s^2 along -x at solar minimum
%
% See also NonGravForces.kernel.gsiCoefficients, NonGravForces.kernel.accommodation, NonGravForces.geom.flowShadow, NonGravForces.kernel.assemble.
%
% KEYWORDS drag per-panel force free molecular species Sentman self-shadowing accommodation
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic); flow shadowing and two-state accommodation 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    macro
    rho double
    species double
    T_K double
    vrel_body double
    options.Model string = string.empty
    options.Alpha = []
    options.Tw_K double = []
    options.Cd double = []
    options.SpecularFraction double = []
    options.Shadow = []
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
[shadowKind, fBeam, fN, fP, fT] = localShadow(unit, options.Shadow, N, K, model);
if any(~isfinite(vrel_body(:))) || any(~isfinite(rho)) || any(~isfinite(species(:))) || any(~isfinite(T_K(:)))
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
        [Cp, Ctau] = NonGravForces.kernel.gsiCoefficients("constantCd", [], theta, [], [], [], 'Cd', options.Cd);
        Cpm = rho .* Cp; Ctm = rho .* Ctau;         % NxK, mass-weighted (single "species")
        sRatio = nan(8, K); rhoSp = rho;
    otherwise
        if isempty(options.Alpha) || isempty(options.Tw_K)
            NonGravForces.fail(unit, "inputModel", "%s needs 'Alpha' and 'Tw_K' -> pass the accommodation (kernel.accommodation) and the wall temperature", model);
        end
        [alphaA, coverage, alphaS] = localAccommodation(unit, options.Alpha, species, T_K, K);
        Tw = options.Tw_K + zeros(1, K);
        if any(Tw <= 0)
            NonGravForces.fail(unit, "inputRange", "Tw_K must be > 0 -> got %g", min(Tw));
        end
        mass = c.speciesMass / c.NA;                % kg per molecule, 1x8
        sRatio = U ./ sqrt(2 * c.kB * T_K ./ mass(:));   % 8xK
        Cpm = zeros(N, K); Ctm = zeros(N, K); rhoSp = zeros(1, K);
        extra = {'Provenance', false};
        if model == "maxwell", extra = [extra, {'SpecularFraction', options.SpecularFraction}]; end
        for j = 1:8
            rhoj = species(j, :) * mass(j);         % 1xK
            rhoSp = rhoSp + rhoj;
            if ~any(rhoj > 0), continue; end
            [~, Ctj, gi] = NonGravForces.kernel.gsiCoefficients(model, sRatio(j, :), theta, alphaA, Tw, T_K, extra{:});
            CpRe = gi.CpReemitted;
            if ~isempty(coverage)
                [~, ~, gs] = NonGravForces.kernel.gsiCoefficients(model, sRatio(j, :), theta, alphaS, Tw, T_K, extra{:});
                CpRe = coverage .* gi.CpReemitted + (1 - coverage) .* gs.CpReemitted;
            end
            switch shadowKind
                case "flow"
                    Cpm = Cpm + rhoj .* (fP(:, :, j) .* gi.CpIncident + fN(:, :, j) .* CpRe);
                    Ctm = Ctm + rhoj .* (fT(:, :, j) .* Ctj);
                otherwise
                    Cpm = Cpm + rhoj .* (gi.CpIncident + CpRe);
                    Ctm = Ctm + rhoj .* Ctj;
            end
        end
end
scale = 0.5 * U.^2 .* A;                            % NxK
if shadowKind == "fraction" || (shadowKind == "flow" && model == "constantCd")
    scale = scale .* fBeam;
end
fn = -scale .* Cpm;                                 % along -n
ft = scale .* Ctm;                                  % along t
for k = 1:N
    F(:, k, :) = reshape(Nrm(:, k) .* fn(k, :), 3, 1, K) + tdir(:, k, :) .* reshape(ft(k, :), 1, 1, K);
end
if nargout > 1
    info.model = model;
    info.gamma = gamma;
    info.theta_rad = theta;
    info.speed_m_s = U;
    info.speedRatio = sRatio;
    info.rhoSpecies_kg_m3 = rhoSp;
    info.Cp = Cpm ./ max(rhoSp, realmin);
    info.Ctau = Ctm ./ max(rhoSp, realmin);
    info.shadow = shadowKind;
    info.provenance = NonGravForces.provenance("kernel.dragPanelForce", 'Options', struct('Model', model, 'Tw_K', options.Tw_K, 'Cd', options.Cd, 'SpecularFraction', options.SpecularFraction, 'shadow', shadowKind));
end
end

function [kind, fBeam, fN, fP, fT] = localShadow(unit, S, N, K, model)
% The three accepted forms of 'Shadow', validated once.
fBeam = []; fN = []; fP = []; fT = [];
if isempty(S)
    kind = "none";
    return
end
if isstruct(S)
    if ~(isscalar(S) && isfield(S, 'kind') && S.kind == "flow" && all(isfield(S, {'beam', 'number', 'normal', 'tangential'})))
        NonGravForces.fail(unit, "inputShape", "a struct 'Shadow' must be a geom.flowShadow result -> build it with geom.flowShadow(macro, vrel_body, T_K, ...)");
    end
    kind = "flow";
    fBeam = S.beam;
    if ~isequal(size(fBeam), [N K])
        NonGravForces.fail(unit, "inputShape", "Shadow.beam must be NxK = %dx%d -> got %s (the flowShadow was built for other epochs or another macro model)", N, K, mat2str(size(fBeam)));
    end
    if model ~= "constantCd"
        fN = S.number; fP = S.normal; fT = S.tangential;
        if ~(isequal(size(fN), [N K 8]) && isequal(size(fP), [N K 8]) && isequal(size(fT), [N K 8]))
            NonGravForces.fail(unit, "inputShape", "Shadow.number/.normal/.tangential must be NxKx8 = %dx%dx8 for %s -> build the flowShadow WITH the gas temperature", N, K, model);
        end
        % the tangential ratio is a signed quantity (geom.flowShadow, RANGES): finite is all it must be
        if any(~isfinite(fT(:)))
            NonGravForces.fail(unit, "numericNaN", "a tangential shadow ratio is NaN/Inf -> check the flowShadow output");
        end
    end
    vals = [fBeam(:); fN(:); fP(:)];
elseif isnumeric(S)
    kind = "fraction";
    fBeam = S;
    if ~isequal(size(fBeam), [N K])
        NonGravForces.fail(unit, "inputShape", "Shadow must be NxK = %dx%d -> got %s", N, K, mat2str(size(fBeam)));
    end
    vals = fBeam(:);
else
    NonGravForces.fail(unit, "inputShape", "'Shadow' must be a geom.flowShadow struct or an NxK fraction array -> got %s", class(S));
end
if any(~isfinite(vals))
    NonGravForces.fail(unit, "numericNaN", "a shadow fraction is NaN/Inf -> check the selfShadow/flowShadow output");
end
if any(vals < 0 | vals > 1)
    NonGravForces.fail(unit, "inputRange", "shadow fractions must lie in [0, 1] -> check the selfShadow/flowShadow output");
end
end

function [alphaA, coverage, alphaS] = localAccommodation(unit, a, species, T_K, K)
% A numeric coefficient is validated through accommodation("fixed"); a struct is taken as built.
if isstruct(a)
    if ~(isscalar(a) && all(isfield(a, {'alpha', 'coverage', 'alphaSubstrate'})) && isequal(size(a.alpha), [1 K]))
        NonGravForces.fail(unit, "inputShape", "a struct 'Alpha' must be a kernel.accommodation result for these %d epochs -> build it with kernel.accommodation", K);
    end
    alphaA = a.alpha; coverage = a.coverage; alphaS = a.alphaSubstrate;
    if ~isempty(coverage) && ~(isequal(size(coverage), [1 K]) && isequal(size(alphaS), [1 K]))
        NonGravForces.fail(unit, "inputShape", "the accommodation struct's coverage and alphaSubstrate must be 1x%d -> build it with kernel.accommodation", K);
    end
    v = [alphaA(:); coverage(:); alphaS(:)];
    if any(~isfinite(v)) || any(v < 0 | v > 1)
        NonGravForces.fail(unit, "inputRange", "accommodation coefficients and coverage must lie in [0, 1] -> check the accommodation struct");
    end
else
    acc = NonGravForces.kernel.accommodation("fixed", species, T_K, 'Alpha', a);
    alphaA = acc.alpha; coverage = []; alphaS = [];
end
end
