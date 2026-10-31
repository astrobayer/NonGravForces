function [F, info] = radiationPanelForce(macro, src, options)
%RADIATIONPANELFORCE Per-panel radiation force for a source list (Sun, Earth tiles) in the body frame.
%
%   [F, info] = NonGravForces.kernel.radiationPanelForce(macro, src)
%   [F, info] = NonGravForces.kernel.radiationPanelForce(macro, src, 'Shadow', f)
%
%   For every plate k, source m (unit direction d from the satellite to the
%   source, pressure P = E/c in Pa, band vis or ir) and epoch, with
%   gamma = n_k . d > 0 (a lit plate) and the band's triplet [a, s, delta]
%   (absorbed, specular, diffuse; a + s + delta = 1),
%
%       F_km = - P_m A_k f_km gamma [ (a + delta) d + ( 2 s gamma + (2/3) delta ) n_k ]
%
%   (theory eq. 6 = Wöske 2019 eq. 5.4; Montenbruck & Gill 2000 eq. 3.73), and
%   F_k = sum_m F_km. The same expression serves solar, albedo and infrared
%   sources through the band (principle P3: sources, not effects); the eclipse
%   factor is already inside the Sun's P (env.radiationSources). Plates facing
%   away from a source (gamma <= 0) receive nothing from it. Thermal re-emission
%   of the absorbed part is NOT included here: it is kernel.thermalRecoil with
%   the panel temperatures of geom.thermalState.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    ----   -----------
%   macro    1x1   MacroModel     Panel model (A, N, optVis [abs;spec;diff], optIR [emis;spec;diff])
%   src      1x1   struct  -      Source list: .dir 3xMxK unit vectors satellite -> source (body frame),
%                                 .P MxK pressure [Pa], .band 1xM string "vis"|"ir" (env.radiationSources)
%
% NAME-VALUE ARGUMENTS
%   Name     Size    Type    Units  Description
%   ----     ----    ----    -----  -----------
%   Shadow   NxKxM   double  -      Self-shadow fraction per plate, epoch and source (geom.selfShadow
%                                   for each source direction), default all 1 (no self-shadowing)
%
% OUTPUT ARGUMENTS
%   Name     Size    Type    Units  Description
%   ----     ----    ----    -----  -----------
%   F        3xNxK   double  N      Force on each plate, body frame, summed over the sources
%   info     1x1     struct  -      .absorbed_W_m2 (NxK, a P c gamma f summed over sources, f the Shadow
%                                   fraction: the flux thermalState needs), .litCount (NxK, sources with
%                                   gamma > 0, counted before the Shadow),
%                                   .M, .K, .bands, .provenance
%
% ACCURACY
%   Exact evaluation of the flat-plate optical model; rounding 1e-16 relative
%   per source. The model error is the flat-plate assumption: Lambertian
%   diffuse reflection, no wavelength dependence inside a band, no
%   inter-reflection between plates (Wöske 2019 quotes ~10 % against a
%   ray-traced finite-element model). Sources at exactly gamma = 0 contribute
%   0 (no NaN).
%
% ERRORS
%   NonGravForces:radiationPanelForce:inputMacro    macro is not a MacroModel
%   NonGravForces:radiationPanelForce:inputSources  src lacks .dir/.P/.band, shapes disagree, band not vis|ir,
%                                                   or a direction is not a unit vector (1e-9)
%   NonGravForces:radiationPanelForce:inputShape    Shadow not NxKxM
%   NonGravForces:radiationPanelForce:inputRange    P < 0 or Shadow outside [0, 1]
%   NonGravForces:radiationPanelForce:numericNaN    NaN/Inf in src or Shadow
%
% REFERENCES
%   Wöske F. (2019) PhD thesis, Univ. Bremen, eqs 5.2-5.8; Montenbruck O., Gill E.
%   (2000) Satellite Orbits, Sect. 3.4 eq. 3.73; Rievers B., Lämmerzahl C. (2011)
%   Ann. Phys. 523:439.
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C"); opt = NonGravForces.Options("standard");
%   src = NonGravForces.env.radiationSources(ep, r, q, opt);
%   [F, info] = NonGravForces.kernel.radiationPanelForce(mm, src);
%   a = squeeze(sum(F, 2)) / 601.214;        % 3xK, ~2e-8 m/s^2 in sunlight
%
% See also NonGravForces.env.radiationSources, NonGravForces.kernel.thermalRecoil, NonGravForces.kernel.assemble, NonGravForces.geom.selfShadow.
%
% KEYWORDS radiation pressure per-panel force albedo infrared flat plate
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    macro
    src (1,1) struct
    options.Shadow double = []
end
unit = "radiationPanelForce";
if ~isa(macro, 'NonGravForces.geom.MacroModel')
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel -> build one with geom.graceFO or geom.boxWing");
end
if ~all(isfield(src, {'dir','P','band'}))
    NonGravForces.fail(unit, "inputSources", "src needs fields dir (3xMxK), P (MxK), band (1xM) -> use env.radiationSources");
end
d = src.dir; P = src.P; band = string(src.band);
if ndims(d) > 3 || size(d, 1) ~= 3
    NonGravForces.fail(unit, "inputSources", "src.dir must be 3xMxK -> got %s", mat2str(size(d)));
end
M = size(d, 2); K = size(d, 3);
if ~isequal(size(P), [M K]) || ~isequal(size(band), [1 M]) || ~all(band == "vis" | band == "ir")
    NonGravForces.fail(unit, "inputSources", "src.P must be MxK = %dx%d and src.band 1xM of vis|ir -> got P %s, band %s", M, K, mat2str(size(P)), mat2str(size(band)));
end
if any(~isfinite(d(:))) || any(~isfinite(P(:)))
    NonGravForces.fail(unit, "numericNaN", "src.dir or src.P contain NaN/Inf -> screen the source list");
end
if any(P(:) < 0)
    NonGravForces.fail(unit, "inputRange", "src.P must be >= 0 -> check the source list");
end
if any(abs(squeeze(vecnorm(d, 2, 1)) - 1) > 1e-9, 'all')
    NonGravForces.fail(unit, "inputSources", "src.dir must hold unit vectors (1e-9) -> normalise the directions");
end
N = macro.n;
if isempty(options.Shadow)
    fShadow = ones(N, K, M);
else
    fShadow = options.Shadow;
    if ~isequal(size(fShadow, 1:3), [N K M])
        NonGravForces.fail(unit, "inputShape", "Shadow must be NxKxM = %dx%dx%d -> got %s", N, K, M, mat2str(size(fShadow)));
    end
    if any(~isfinite(fShadow(:)))
        NonGravForces.fail(unit, "numericNaN", "Shadow contains NaN/Inf -> check the selfShadow output");
    end
    if any(fShadow(:) < 0 | fShadow(:) > 1)
        NonGravForces.fail(unit, "inputRange", "Shadow fractions must lie in [0, 1] -> check the selfShadow output");
    end
end
Nrm = macro.N; A = macro.A(:);
Fx = zeros(N, K); Fy = Fx; Fz = Fx;
absorbed = zeros(N, K); lit = zeros(N, K);
cLight = referenceSystem.iersConstants().speedOfLight;
for m = 1:M
    dm = reshape(d(:, m, :), 3, K);                 % 3xK
    if band(m) == "vis", opt = macro.optVis; else, opt = macro.optIR; end   % 3xN [abs; spec; diff]
    gamma = Nrm' * dm;                              % NxK
    gp = max(gamma, 0);
    w = P(m, :) .* A .* fShadow(:, :, m) .* gp;     % NxK, P A f gamma
    lit = lit + (gamma > 0);
    absorbed = absorbed + opt(1, :)' .* w * cLight ./ A;   % a E gamma f  [W/m^2]
    % Every plate at once, one NxK accumulator per component: element by element the same products
    % and sums as one plate at a time, so F is bitwise what it was. The loop over plates writing
    % strided slices of a 3xNxK F took 40 % longer on 577 sources (HANDOVER, 2026-10-16).
    ad = opt(1, :) + opt(3, :);                                      % 1xN
    cn = (2 * opt(2, :)') .* gp + (2 / 3) * opt(3, :)';              % NxK
    Fx = Fx - w .* (ad' .* dm(1, :) + cn .* Nrm(1, :)');
    Fy = Fy - w .* (ad' .* dm(2, :) + cn .* Nrm(2, :)');
    Fz = Fz - w .* (ad' .* dm(3, :) + cn .* Nrm(3, :)');
end
F = permute(cat(3, Fx, Fy, Fz), [3 1 2]);                          % 3xNxK
info.absorbed_W_m2 = absorbed;
info.litCount = lit;
info.M = M; info.K = K; info.bands = band;
info.provenance = NonGravForces.provenance("kernel.radiationPanelForce", 'Options', struct('shadowGiven', ~isempty(options.Shadow), 'M', M));
end
