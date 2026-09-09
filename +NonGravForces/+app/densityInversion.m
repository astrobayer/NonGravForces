function [T, info] = densityInversion(model, epoch, r, v, q, aMeasured, options)
%DENSITYINVERSION Thermosphere density from calibrated non-gravitational accelerations.
%
%   T = NonGravForces.app.densityInversion(model, epoch, r, v, q, aMeasured)
%   [T, info] = NonGravForces.app.densityInversion(___, 'Frame', "body", 'Chunk', 200)
%
%   THE INVERSION IS A DIVISION, NOT A SEARCH, and that is a property of the
%   physics rather than a simplification. The gas force is exactly linear in the
%   number densities (every species term carries rho_j as a factor), so one
%   evaluation of the model at its own density gives the whole sensitivity:
%
%       a_drag(rho) = (rho / rho_model) a_drag(rho_model)
%       s           = <a_measured - a_radiation - a_thermal, u> / <a_drag(rho_model), u>
%       rho         = s rho_model
%
%   with u the UNIT VECTOR ALONG THE MODELLED DRAG. Projecting is not a way of
%   discarding information; it is the recognition that only one direction of the
%   measurement carries density at all. The component perpendicular to u is
%   returned as `residual_perp` and is the honest quality indicator: it is what
%   the density cannot explain, and it grows when the attitude, the radiation
%   model or the calibration is wrong.
%
%   WHAT THIS DOES NOT DO. It does not calibrate: `aMeasured` must already be
%   bias- and scale-corrected (app.calibrate, or sensor.Electrostatic inverted by
%   the caller), because a bias is indistinguishable from a density offset over a
%   short arc and pretending otherwise would launder one into the other. It does
%   not separate drag from radiation: it SUBTRACTS the modelled radiation and
%   thermal recoil, so an error in those goes straight into rho, which is why the
%   perpendicular residual is reported beside every value.
%
%   THE SCALE MAY COME OUT NEGATIVE, and it is returned rather than clipped. A
%   negative scale means the measurement points against the modelled drag, which
%   happens in deep solar minimum when the drag falls below the radiation
%   modelling error, and it is a statement about the data. Clipping at zero would
%   turn a visible failure into a plausible-looking density floor.
%
% INPUT ARGUMENTS
%   Name        Size    Type    Units   Description
%   ----        ----    ----    -----   -----------
%   model       1x1     Model   -       A configured NonGravForces.Model
%   epoch       1x1     struct  -       referenceSystem epoch struct with K entries
%   r           3xK     double  m       Position, GCRS
%   v           3xK     double  m/s     Velocity, GCRS
%   q           4xK     double  -       BODY->GCRS quaternions, or [] for the model's provider
%   aMeasured   3xK     double  m/s^2   Calibrated non-gravitational acceleration, in 'Frame'
%
% NAME-VALUE ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   Frame     1x1   string  -      "body" (default) or "gcrs"; aMeasured is read in it
%   Chunk     1x1   double  -      Epochs per evaluation, >= 1; default the model's
%                                  Options.series.chunk. Does not change the result
%   MinDrag   1x1   double  m/s^2  Below this modelled drag magnitude the epoch is marked
%                                  unusable instead of dividing by it; default 1e-12
%
% OUTPUT ARGUMENTS
%   Name   Size    Type       Units    Description
%   ----   ----    ----       -----    -----------
%   T      Kx7     timetable  -        .rho_kg_m3 inverted density, .rho_model_kg_m3 the
%                                      backend's own, .scale their ratio, .dragProjected_m_s2
%                                      the measurement along u, .dragModel_m_s2 the model along
%                                      u, .residualPerp_m_s2 the norm of what is left across u,
%                                      .usable logical (false where the modelled drag is below
%                                      MinDrag, where rho_model is zero, or where the epoch was
%                                      not finite). rho is NaN wherever usable is false
%   info   1x1     struct     -        .K .chunks .frame .usable (count) .provenance
%
% ACCURACY
%   Exact by construction where the model is the truth: feeding a model's own
%   acceleration back in returns its own density, and the round trip is the
%   reference test. Against real data the accuracy is that of everything
%   subtracted first - the radiation model, the macro model, the attitude and the
%   calibration - and no bound here can speak for those.
%
% ERRORS
%   NonGravForces:densityInversion:inputModel   model is not a NonGravForces.Model
%   NonGravForces:densityInversion:inputShape   r, v, q or aMeasured do not match the epoch count
%   NonGravForces:densityInversion:inputRange   Chunk < 1, or MinDrag < 0, or not finite
%   NonGravForces:densityInversion:numericNaN   aMeasured contains NaN or Inf
%   NonGravForces:Model:*                       propagated from the evaluation
%
% REFERENCES
%   Doornbos E. (2011) Thermospheric density and wind determination from
%   satellite dynamics, PhD, TU Delft, ch. 4 (the iterative form; the direct
%   division used here follows from the exact linearity in rho);
%   Bruinsma S. et al. (2004) J. Geophys. Res. 109:A02307; Sutton E. (2008)
%   Effects of solar disturbances on the thermosphere densities, PhD, U Colorado.
%   HANDOVER.md Part G LP-7 (WP-32).
%
% EXAMPLE
%   [T, info] = NonGravForces.app.densityInversion(M, ep, r, v, [], aCalibrated);
%   plot(T.Time(T.usable), T.rho_kg_m3(T.usable))
%
% See also NonGravForces.Model, NonGravForces.app.forceBudget, NonGravForces.env.density, NonGravForces.app.calibrate.
%
% KEYWORDS density inversion thermosphere drag accelerometer Doornbos
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
arguments
    model
    epoch (1,1) struct
    r double
    v double
    q double = []
    aMeasured double = []
    options.Frame (1,1) string {mustBeMember(options.Frame, ["body","gcrs"])} = "body"
    options.Chunk double = []
    options.MinDrag (1,1) double = 1e-12
end
unit = "densityInversion";
if ~isa(model, 'NonGravForces.Model')
    NonGravForces.fail(unit, "inputModel", "model must be a NonGravForces.Model -> build one with NonGravForces.Model(macro, options, ...)");
end
K = numel(epoch.sec);
if ~isequal(size(r), [3 K]) || ~isequal(size(v), [3 K])
    NonGravForces.fail(unit, "inputShape", "r and v must be 3x%d for this epoch -> got %s and %s", K, mat2str(size(r)), mat2str(size(v)));
end
if ~isempty(q) && ~isequal(size(q), [4 K])
    NonGravForces.fail(unit, "inputShape", "q must be 4x%d or [] -> got %s", K, mat2str(size(q)));
end
if ~isequal(size(aMeasured), [3 K])
    NonGravForces.fail(unit, "inputShape", "aMeasured must be 3x%d, one calibrated acceleration per epoch -> got %s", K, mat2str(size(aMeasured)));
end
if any(~isfinite(aMeasured(:)))
    NonGravForces.fail(unit, "numericNaN", "aMeasured contains NaN/Inf -> screen or gap-fill the calibrated series before inverting it");
end
if ~isfinite(options.MinDrag) || options.MinDrag < 0
    NonGravForces.fail(unit, "inputRange", "'MinDrag' must be finite and non-negative -> it is the modelled drag below which an epoch is refused rather than divided by");
end
if isempty(options.Chunk)
    chunk = model.options.series.chunk;
else
    chunk = options.Chunk;
end
if ~isscalar(chunk) || ~isfinite(chunk) || chunk < 1
    NonGravForces.fail(unit, "inputRange", "'Chunk' must be a finite scalar >= 1 -> got %s", mat2str(chunk));
end
chunk = max(1, round(chunk));

first = 1:chunk:K;
last = min(first + chunk - 1, K);
nC = numel(first);
rhoModel = zeros(1, K); dragProj = zeros(1, K); dragMod = zeros(1, K);
perp = zeros(1, K); usable = false(1, K);
for c = 1:nC
    idx = first(c):last(c);
    ep = referenceSystem.daySecondToEpoch([epoch.day(idx) epoch.sec(idx)], epoch.scale);
    qc = []; if ~isempty(q), qc = q(:, idx); end
    [~, ai] = model.accel(ep, r(:, idx), v(:, idx), qc, 'Frame', options.Frame);
    % What the measurement must explain once everything that is not gas is removed.
    residual = aMeasured(:, idx) - ai.radiation - ai.thermal;
    dm = ai.drag;
    nrm = vecnorm(dm);
    ok = nrm > options.MinDrag & isfinite(ai.rho) & ai.rho > 0;
    u = dm ./ max(nrm, realmin);                 % the direction the density lives on
    dragProj(idx) = sum(residual .* u, 1);
    dragMod(idx) = nrm;
    perp(idx) = vecnorm(residual - u .* sum(residual .* u, 1));
    rhoModel(idx) = ai.rho;
    usable(idx) = ok;
end
scale = dragProj ./ max(dragMod, realmin);       % nrm cancels: <res,u>/|dm| = <res,dm>/|dm|^2
rho = scale .* rhoModel;
scale(~usable) = NaN; rho(~usable) = NaN;

t = referenceSystem.epochTo(epoch, "datetime");
T = timetable(t(:), rho(:), rhoModel(:), scale(:), dragProj(:), dragMod(:), perp(:), usable(:), ...
    'VariableNames', {'rho_kg_m3', 'rho_model_kg_m3', 'scale', 'dragProjected_m_s2', ...
                      'dragModel_m_s2', 'residualPerp_m_s2', 'usable'});
T.Properties.VariableUnits = ["kg/m^3", "kg/m^3", "", "m/s^2", "m/s^2", "m/s^2", ""];
T.Properties.Description = sprintf('NonGravForces.app.densityInversion, %s frame, %d epochs in %d chunk(s), %d usable', ...
    options.Frame, K, nC, nnz(usable));
info = struct('K', K, 'chunks', nC, 'chunkSize', chunk, 'frame', options.Frame, ...
    'usable', nnz(usable), 'minDrag', options.MinDrag, ...
    'provenance', NonGravForces.provenance("app.densityInversion", 'Options', ...
        struct('Frame', options.Frame, 'Chunk', chunk, 'MinDrag', options.MinDrag)));
end
