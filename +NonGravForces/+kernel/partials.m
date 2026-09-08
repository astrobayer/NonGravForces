function [dF, info] = partials(what, wrt, macro, options)
%PARTIALS Analytic derivatives of the per-panel kernel forces with respect to one scalar parameter.
%
%   dF = NonGravForces.kernel.partials("drag", "alpha", macro, 'Rho', rho, 'Species', sp, ...
%                                      'T_K', T, 'Vrel_body', v, 'Model', "sentman", 'Alpha', a, 'Tw_K', tw)
%   dF = NonGravForces.kernel.partials("drag", "speed", macro, ___)
%   dF = NonGravForces.kernel.partials("radiation", "pressureScale", macro, 'Sources', src)
%   dF = NonGravForces.kernel.partials("thermal", "temperature", macro, 'Temperature_K', T)
%
%   Returns d(F_k)/d(parameter) for every panel and epoch, in the same 3xNxK
%   layout and body frame as the force itself, so a caller differentiates the
%   assembled force by running kernel.assemble on the derivative (both are
%   linear reductions). This unit differentiates only what has a closed form;
%   the state derivatives d/dr and d/dv run through the atmosphere model and a
%   MEX gateway and are taken numerically by Model.partials.
%
%   PARAMETERS
%     "drag" / "alpha"          the energy accommodation coefficient. Only the
%                               re-emitted normal momentum depends on it, so the
%                               tangential derivative is identically zero.
%     "drag" / "densityScale"   a multiplicative scale on every species number
%                               density. The force is linear in it, so the
%                               derivative is the force (this is the sensitivity
%                               app.densityInversion inverts).
%     "drag" / "speed"          the magnitude of the relative velocity with its
%                               direction held fixed. Both the dynamic pressure
%                               (quadratic) and every species speed ratio
%                               (linear) depend on it.
%     "radiation" / "pressureScale"  a multiplicative scale on every source
%                               pressure; the force is linear in it.
%     "thermal" / "temperature" each panel's own temperature. The force goes as
%                               T^4, so dF_k/dT_k = 4 F_k / T_k and the
%                               derivative with respect to another panel is zero;
%                               only the diagonal is returned.
%
%   The coefficient derivatives come from kernel.gsiCoefficients with
%   'Derivatives', true, where they sit beside the closed form they
%   differentiate, so the two cannot drift apart (rule R1).
%
% INPUT ARGUMENTS
%   Name    Size  Type    Units  Description
%   ----    ----  ----    -----  -----------
%   what    1x1   string  -      "drag" | "radiation" | "thermal"
%   wrt     1x1   string  -      the parameter, from the list above
%   macro   1x1   MacroModel     the panel model the force was computed with
%
% NAME-VALUE ARGUMENTS
%   Name              Size    Type    Units   Description
%   ----              ----    ----    -----   -----------
%   Rho               1xK     double  kg/m^3  total mass density (drag)
%   Species           8xK     double  1/m^3   number densities (drag; not for constantCd)
%   T_K               1xK     double  K       gas kinetic temperature (drag)
%   Vrel_body         3xK     double  m/s     relative velocity, body frame (drag)
%   Model             1x1     string  -       "sentman" | "maxwell" | "constantCd" (drag)
%   Alpha             1x1|1xK double  -       accommodation coefficient (drag, kinetic models)
%   Tw_K              1x1|1xK double  K       wall temperature (drag, kinetic models)
%   Cd                1x1     double  -       drag coefficient (drag, constantCd)
%   SpecularFraction  1x1     double  -       specular fraction (drag, maxwell)
%   Shadow            NxK     double  -       self-shadow fractions (drag), default all 1
%   Sources           1x1     struct  -       source list (radiation), as env.radiationSources
%   RadShadow         NxKxM   double  -       self-shadow fractions (radiation), default all 1
%   Temperature_K     NxK     double  K       panel temperatures (thermal)
%
% OUTPUT ARGUMENTS
%   Name    Size    Type    Units          Description
%   ----    ----    ----    -----          -----------
%   dF      3xNxK   double  N per unit of  derivative of the per-panel force
%                           the parameter
%   info    1x1     struct  -              .what, .wrt, .unit (the parameter's unit as a
%                                          string), .K, .N, .provenance
%
% ACCURACY
%   Closed forms differentiated term by term; rounding only, 1e-15 relative.
%   Verified against Richardson-extrapolated central differences of the force
%   itself: measured 2026-09-08, worst 3e-10 relative over the test grid. A
%   complex-step verifier is NOT possible on this path and the toolbox does not
%   pretend otherwise: MATLAB's erf and erfc raise MATLAB:erf:notFullReal on a
%   complex argument and abs/vecnorm discard the imaginary part, so no complex
%   step survives Sentman's coefficients or the speed (HANDOVER Part I,
%   2026-09-08).
%
% ERRORS
%   NonGravForces:partials:inputMacro     macro is not a MacroModel
%   NonGravForces:partials:unsupportedWrt the parameter is not defined for this force
%   NonGravForces:partials:inputMissing   a name-value this combination needs was not given
%   NonGravForces:partials:inputShape     an input has the wrong shape
%   NonGravForces:gsiCoefficients:*       propagated from the coefficient model
%
% REFERENCES
%   Sentman (1961); Doornbos (2011) Ch. 3; HANDOVER WP-16 and the Part I decision
%   of 2026-09-08 on the verifier.
%
% EXAMPLE
%   dFa = NonGravForces.kernel.partials("drag", "alpha", mm, 'Rho', rho, 'Species', sp, ...
%           'T_K', T, 'Vrel_body', vb, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
%   da = NonGravForces.kernel.assemble(dFa, mm) / mass;      % d(acceleration)/d(alpha)
%
% See also NonGravForces.kernel.gsiCoefficients, NonGravForces.kernel.dragPanelForce, NonGravForces.Model.
%
% KEYWORDS partials derivatives sensitivity accommodation density estimation
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
arguments
    what (1,1) string {mustBeMember(what, ["drag","radiation","thermal"])}
    wrt (1,1) string
    macro
    options.Rho double = []
    options.Species double = []
    options.T_K double = []
    options.Vrel_body double = []
    options.Model string = string.empty
    options.Alpha double = []
    options.Tw_K double = []
    options.Cd double = []
    options.SpecularFraction double = []
    options.Shadow double = []
    options.Sources = []
    options.RadShadow double = []
    options.Temperature_K double = []
end
unit = "partials";
if ~isa(macro, 'NonGravForces.geom.MacroModel')
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel -> the same model the force was computed with");
end
allowed = struct('drag', ["alpha","densityScale","speed"], 'radiation', "pressureScale", 'thermal', "temperature");
if ~any(wrt == allowed.(what))
    NonGravForces.fail(unit, "unsupportedWrt", "%s has no closed-form derivative with respect to '%s' -> use one of %s (d/dr and d/dv are numerical, in Model.partials)", ...
        what, wrt, strjoin(allowed.(what), ", "));
end
N = macro.n;
info.what = what; info.wrt = wrt; info.N = N;

switch what
    case "thermal"
        T = options.Temperature_K;
        if isempty(T)
            NonGravForces.fail(unit, "inputMissing", "thermal needs 'Temperature_K' (NxK) -> pass the panel temperatures the force was computed with");
        end
        F = NonGravForces.kernel.thermalRecoil(macro, T);
        K = size(T, 2);
        dF = 4 * F ./ reshape(max(T, realmin), 1, N, K);     % F ~ T^4, panel by panel
        info.unit = "N/K";
        info.K = K;

    case "radiation"
        src = options.Sources;
        if isempty(src)
            NonGravForces.fail(unit, "inputMissing", "radiation needs 'Sources' -> pass the source list the force was computed with (env.radiationSources)");
        end
        extra = {};
        if ~isempty(options.RadShadow), extra = {'Shadow', options.RadShadow}; end
        dF = NonGravForces.kernel.radiationPanelForce(macro, src, extra{:});   % linear in every source pressure
        info.unit = "N per unit scale";
        info.K = size(dF, 3);

    case "drag"
        for f = ["Rho", "Vrel_body"]
            if isempty(options.(f))
                NonGravForces.fail(unit, "inputMissing", "drag needs '%s' -> pass the inputs the force was computed with", f);
            end
        end
        if isempty(options.Model) || ~isscalar(options.Model)
            NonGravForces.fail(unit, "inputMissing", "drag needs 'Model' -> sentman, maxwell or constantCd");
        end
        model = options.Model;
        vrel = options.Vrel_body; rho = options.Rho;
        K = size(vrel, 2);
        if ~isequal(size(vrel), [3 K]) || ~isequal(size(rho), [1 K])
            NonGravForces.fail(unit, "inputShape", "Vrel_body must be 3xK and Rho 1xK -> got %s and %s", mat2str(size(vrel)), mat2str(size(options.Rho)));
        end
        shadow = options.Shadow;
        if isempty(shadow), shadow = ones(N, K); end
        if ~isequal(size(shadow), [N K])
            NonGravForces.fail(unit, "inputShape", "Shadow must be NxK = %dx%d -> got %s", N, K, mat2str(size(shadow)));
        end
        % the force itself, for the derivatives that are proportional to it
        dragArgs = {'Model', model, 'Alpha', options.Alpha, 'Tw_K', options.Tw_K, 'Cd', options.Cd, ...
                    'SpecularFraction', options.SpecularFraction, 'Shadow', shadow};
        F = NonGravForces.kernel.dragPanelForce(macro, rho, options.Species, options.T_K, vrel, dragArgs{:});
        switch wrt
            case "densityScale"
                dF = F;                                  % linear in the density
                info.unit = "N per unit scale";
            otherwise
                [dF, info] = localCoefficientDerivative(unit, info, wrt, macro, F, options, model, vrel, K, N, shadow);
        end
        info.K = K;
end
info.provenance = NonGravForces.provenance("kernel.partials", 'Options', struct('what', what, 'wrt', wrt));
end

% ---------------------------------------------------------------------------------------
function [dF, info] = localCoefficientDerivative(unit, info, wrt, macro, F, options, model, vrel, K, N, shadow)
%LOCALCOEFFICIENTDERIVATIVE d/dalpha and d/dspeed, which run through the gas-surface coefficients.
if model == "constantCd"
    switch wrt
        case "alpha"
            dF = zeros(size(F));                          % constantCd has no accommodation coefficient
            info.unit = "N";
        case "speed"
            dF = 2 * F ./ reshape(vecnorm(vrel), 1, 1, K); % only the dynamic pressure depends on U
            info.unit = "N/(m/s)";
    end
    return
end
c = NonGravForces.const();
species = options.Species; T_K = options.T_K;
if isempty(species) || isempty(T_K)
    NonGravForces.fail(unit, "inputMissing", "%s needs 'Species' and 'T_K' for a kinetic model -> pass what env.density returned", model);
end
U = vecnorm(vrel);
v = vrel ./ U;
Nrm = macro.N; A = macro.A(:);
gamma = Nrm' * v;                                          % NxK
l = zeros(N, K); tdir = zeros(3, N, K);
for k = 1:N
    nk = Nrm(:, k);
    cr = [nk(2) * v(3, :) - nk(3) * v(2, :); nk(3) * v(1, :) - nk(1) * v(3, :); nk(1) * v(2, :) - nk(2) * v(1, :)];
    l(k, :) = vecnorm(cr);
    t = -v + gamma(k, :) .* nk;
    ok = l(k, :) > 0;
    t(:, ok) = t(:, ok) ./ l(k, ok);
    t(:, ~ok) = 0;
    tdir(:, k, :) = reshape(t, 3, 1, K);
end
theta = atan2(l, gamma);
alpha = NonGravForces.kernel.accommodation("fixed", species, T_K, 'Alpha', options.Alpha);
Tw = options.Tw_K + zeros(1, K);
mass = c.speciesMass / c.NA;
sRatio = U ./ sqrt(2 * c.kB * T_K ./ mass(:));             % 8xK
extra = {};
if model == "maxwell", extra = {'SpecularFraction', options.SpecularFraction}; end
dCpm = zeros(N, K); dCtm = zeros(N, K);
for j = 1:8
    rhoj = species(j, :) * mass(j);
    if ~any(rhoj > 0), continue; end
    [~, ~, gi] = NonGravForces.kernel.gsiCoefficients(model, sRatio(j, :), theta, alpha, Tw, T_K, extra{:}, 'Derivatives', true);
    switch wrt
        case "alpha"
            dCpm = dCpm + rhoj .* gi.dCp_dalpha;
            dCtm = dCtm + rhoj .* gi.dCtau_dalpha;
        case "speed"
            % s_j is proportional to U at fixed temperature, so ds_j/dU = s_j/U
            dCpm = dCpm + rhoj .* gi.dCp_ds .* (sRatio(j, :) ./ U);
            dCtm = dCtm + rhoj .* gi.dCtau_ds .* (sRatio(j, :) ./ U);
    end
end
scale = 0.5 * U.^2 .* A .* shadow;                         % NxK, the same factor dragPanelForce uses
dF = zeros(3, N, K);
for k = 1:N
    normalPart = Nrm(:, k) .* (-scale(k, :) .* dCpm(k, :));
    tangPart = tdir(:, k, :) .* reshape(scale(k, :) .* dCtm(k, :), 1, 1, K);
    dF(:, k, :) = reshape(normalPart, 3, 1, K) + tangPart;
end
if wrt == "speed"
    dF = dF + 2 * F ./ reshape(U, 1, 1, K);                % the dynamic pressure's own U dependence
    info.unit = "N/(m/s)";
else
    info.unit = "N";
end
end
