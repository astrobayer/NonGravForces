function [T, info] = thermalState(macro, src, options)
%THERMALSTATE Per-panel surface temperatures from the radiation source list: steady state or transient lumped wall.
%
%   [T, info] = NonGravForces.geom.thermalState(macro, src)                          % steady state
%   [T, info] = NonGravForces.geom.thermalState(macro, src, 'Model', "transient", 'Epoch', ep, 'T0_K', T0)
%   [T, info] = NonGravForces.geom.thermalState(macro, src, 'InternalPower_W', 250, 'InternalShare', s)
%
%   Absorbed flux per panel and epoch from the body-frame source list of
%   env.radiationSources (Wöske 2019 eq. 5.15):
%       q_abs,k = sum_m  alpha_k(band_m) * c * P_m * max(0, n_k . dir_m) * o_km
%   with alpha the visible absorptivity for "vis" sources and the infrared
%   emissivity for "ir" sources (Kirchhoff), o_km the unobstructed fraction of
%   plate k toward source m from geom.selfShadow ('Exposed', one value PER
%   SOURCE), plus an internal power share q_int,k = InternalPower_W * InternalShare_k / A_k.
%   The occlusion is applied source by source, as the radiation force applies it:
%   until 2026-09-15 one fraction per panel was applied to the flux of every
%   source, and the Sun's fraction zeroed the Earth heating of every panel facing
%   away from the Sun (HANDOVER LP-C, finding F2).
%
%   steadyState  (Wöske eq. 5.17) eps_k sigma T_k^4 = q_abs,k + q_int,k, so
%                T_k = ((q_abs + q_int) / (eps sigma))^(1/4); a panel with no
%                absorbed power is at 0 K (there is no background term:
%                the 2.7 K sky is 1e-5 of the fluxes in play).
%   transient    (Wöske eq. 5.18) C_k dT_k/dt = q_abs,k + q_int,k - eps_k sigma T_k^4
%                with the per-area heat capacity C_k of the macro model
%                (geom.graceFO carries Wöske Table 6.1 stacks), integrated
%                across the epochs from 'T0_K' by classical RK4 with substeps
%                chosen from the linearised time constant (dt <= tau/10); the
%                sources are held constant over each epoch interval.
%
% INPUT ARGUMENTS
%   Name     Size  Type                          Units  Description
%   ----     ----  ----                          -----  -----------
%   macro    1x1   NonGravForces.geom.MacroModel  -      Optics (optVis(1,:), optIR(1,:)), areas, heat capacities
%   src      1x1   struct                         -      env.radiationSources output: .dir 3xMxK, .P MxK, .band 1xM
%
% NAME-VALUE ARGUMENTS
%   Name             Size  Type    Units  Description
%   ----             ----  ----    -----  -----------
%   Model            1x1   string  -      "steadyState" (default) | "transient"
%   Exposed          NxKxM double  -      Unobstructed fraction of every plate toward every
%                                         source, from geom.selfShadow (default: none = 1);
%                                         N x K when the list has one source
%   InternalPower_W  1x1   double  W      Total internal dissipation radiated through the
%                                         panels (default 0)
%   InternalShare    1xN   double  -      Share of the internal power per panel (sums to 1);
%                                         REQUIRED when InternalPower_W > 0 (rule R5)
%   Epoch            1x1   struct  -      referenceSystem epoch struct, K entries (transient)
%   T0_K             1xN   double  K      Initial temperatures (transient)
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   T        NxK   double  K      Surface temperature per panel and epoch
%   info     1x1   struct  -      .model, .absorbed_W_m2 (NxK q_abs), .internal_W_m2 (NxK),
%                                 .emitted_W_m2 (NxK eps sigma T^4), .heatCapacity_J_m2K (1xN),
%                                 .substeps (1x(K-1), transient), .timeConstant_s (NxK,
%                                 transient, linearised C/(4 eps sigma T^3)), .provenance
%
% ACCURACY
%   Steady state is exact for its own assumption (no conduction between or
%   through walls, no thermal inertia); Wöske (2019) shows the transient
%   models change the TRP acceleration at the 1e-10 m/s^2 level across
%   eclipses. Transient: RK4 with dt <= tau/10 has a relative step error
%   below 1e-5 per interval (tested against the linearised exponential
%   relaxation to 1e-3 for a 1 K excursion and against the steady state for
%   t >> tau); the heat capacities are educated guesses (graceFO help), which
%   is the dominant uncertainty, not the integration.
%
% ERRORS
%   NonGravForces:thermalState:inputMacro      not a MacroModel; transient with NaN heat capacity
%   NonGravForces:thermalState:inputSources    src lacks dir/P/band or its sizes disagree
%   NonGravForces:thermalState:inputShape      Exposed not NxKxM in [0,1], T0_K not 1xN, Epoch not K entries
%   NonGravForces:thermalState:inputInternal   InternalPower_W > 0 without a share summing to 1
%   NonGravForces:thermalState:numericNaN      NaN/Inf in the sources or the result
%
% REFERENCES
%   Wöske F. (2019) PhD thesis, Univ. Bremen, sect. 5.1.4 eqs 5.14-5.18, Table 6.1.
%   Rievers B., Lämmerzahl C. (2011) Ann. Phys. 523:439 (thermal recoil modelling).
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C"); [src] = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard"));
%   T = NonGravForces.geom.thermalState(mm, src);            % 8xK, ~250-330 K in sunlight
%
% See also NonGravForces.env.radiationSources, NonGravForces.geom.selfShadow, NonGravForces.kernel.thermalRecoil.
%
% KEYWORDS thermal temperature steady state transient TRP
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    macro
    src (1,1) struct
    options.Model (1,1) string {mustBeMember(options.Model, ["steadyState","transient"])} = "steadyState"
    options.Exposed double = []
    options.InternalPower_W (1,1) double {mustBeNonnegative} = 0
    options.InternalShare double = []
    options.Epoch = []
    options.T0_K double = []
end
unit = "thermalState";
if ~isa(macro, 'NonGravForces.geom.MacroModel')
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel -> geom.graceFO or geom.boxWing");
end
if ~all(isfield(src, ["dir","P","band"])) || ndims(src.dir) > 3 || size(src.dir, 1) ~= 3 || ~isequal(size(src.P), [size(src.dir, 2), size(src.dir, 3)]) || numel(src.band) ~= size(src.dir, 2)
    NonGravForces.fail(unit, "inputSources", "src must be env.radiationSources output (.dir 3xMxK, .P MxK, .band 1xM) -> build it there");
end
if any(~isfinite(src.dir(:))) || any(~isfinite(src.P(:)))
    NonGravForces.fail(unit, "numericNaN", "src contains NaN/Inf -> screen the source list");
end
N = macro.n; M = size(src.dir, 2); K = size(src.dir, 3);
c = referenceSystem.iersConstants().speedOfLight;
sigma = NonGravForces.const().sigmaSB;
isVis = reshape(string(src.band) == "vis", M, 1);
alphaVis = macro.optVis(1, :)'; epsIR = macro.optIR(1, :)';                  % Nx1
% cosines n_k . dir_m for every k, m, epoch: NxMxK
cosI = max(0, reshape(pagemtimes(macro.N', src.dir), N, M, K));
E = reshape(c * src.P, 1, M, K);                                           % irradiance W/m^2, 1xMxK
absorptivity = alphaVis .* isVis' + epsIR .* (~isVis');                    % NxM
w = absorptivity .* cosI .* E;                                             % NxMxK
if ~isempty(options.Exposed)
    X = options.Exposed;
    if ~isequal(size(X, 1:3), [N K M]) || any(~isfinite(X(:))) || any(X(:) < 0 | X(:) > 1)
        NonGravForces.fail(unit, "inputShape", "Exposed must be %dx%dx%d (plates x epochs x sources) in [0,1] -> geom.selfShadow for every source direction", N, K, M);
    end
    w = w .* permute(X, [1 3 2]);                                          % per source, before the sum
end
qAbs = reshape(sum(w, 2), N, K);                                           % NxK
if options.InternalPower_W > 0
    s = options.InternalShare;
    if ~(numel(s) == N && all(s >= 0) && abs(sum(s) - 1) < 1e-9)
        NonGravForces.fail(unit, "inputInternal", "InternalPower_W > 0 needs InternalShare (1x%d, sums to 1) -> e.g. 0.9 on the radiator and the rest spread (Wöske 2019)", N);
    end
    qInt = repmat(options.InternalPower_W * s(:) ./ macro.A(:), 1, K);
else
    qInt = zeros(N, K);
end
q = qAbs + qInt;
info.model = options.Model;
info.absorbed_W_m2 = qAbs;
info.internal_W_m2 = qInt;
info.heatCapacity_J_m2K = macro.heatCapacity;
switch options.Model
    case "steadyState"
        T = (q ./ (epsIR * sigma)).^(1/4);
    case "transient"
        Cw = macro.heatCapacity(:);
        if any(~isfinite(Cw)) || any(Cw <= 0)
            NonGravForces.fail(unit, "inputMacro", "transient needs a positive heatCapacity_J_m2K on every panel -> geom.graceFO carries Wöske's; boxWing takes 'HeatCapacity_J_m2K'");
        end
        ep = options.Epoch;
        if ~(isstruct(ep) && all(isfield(ep, ["day","sec","scale"])) && numel(ep.day) == K)
            NonGravForces.fail(unit, "inputShape", "transient needs 'Epoch' with %d entries matching the source list -> the epochs of radiationSources", K);
        end
        T0 = options.T0_K(:);
        if numel(T0) ~= N || any(~isfinite(T0)) || any(T0 < 0)
            NonGravForces.fail(unit, "inputShape", "transient needs 'T0_K' (1x%d, finite, >= 0) -> e.g. the steady state of the first epoch", N);
        end
        t = double(ep.day(:)') * 86400 + double(ep.sec(:)');
        if any(diff(t) <= 0)
            NonGravForces.fail(unit, "inputShape", "Epoch must be strictly increasing for the transient integration -> sort the epochs");
        end
        T = zeros(N, K); T(:, 1) = T0;
        substeps = zeros(1, max(K - 1, 0));
        for k = 1:K-1
            dt = t(k + 1) - t(k);
            qk = q(:, k);                                          % sources held over the interval
            Tk = T(:, k);
            % linearised time constant tau = C / (4 eps sigma T^3) at the larger of the current and the steady temperature
            Tref = max(Tk, (max(qk, 0) ./ (epsIR * sigma)).^(1/4));
            tau = Cw ./ (4 * epsIR * sigma .* max(Tref, 1).^3);
            n = max(1, ceil(dt / (min(tau) / 10)));
            substeps(k) = n;
            h = dt / n;
            f = @(Tv) (qk - epsIR * sigma .* Tv.^4) ./ Cw;
            for i = 1:n
                k1 = f(Tk); k2 = f(Tk + h/2 * k1); k3 = f(Tk + h/2 * k2); k4 = f(Tk + h * k3);
                Tk = Tk + h / 6 * (k1 + 2*k2 + 2*k3 + k4);
            end
            T(:, k + 1) = max(0, Tk);
        end
        info.substeps = substeps;
        info.timeConstant_s = Cw ./ (4 * epsIR * sigma .* max(T, 1).^3);
end
if any(~isfinite(T(:)))
    NonGravForces.fail(unit, "numericNaN", "a temperature is NaN/Inf -> check the optics (eps > 0) and the sources");
end
info.emitted_W_m2 = epsIR * sigma .* T.^4;
info.provenance = NonGravForces.provenance("geom.thermalState", 'Options', struct('model', options.Model, 'internalPower_W', options.InternalPower_W));
end
