%% Tutorial 12 - Driving an orbit integrator
% NonGravForces computes forces; it does not integrate orbits. An integrator - here gjToolbox, the
% Gauss-Jackson integrator of odeSuite - needs the non-gravitational acceleration as a right-hand
% side it can call at any time and state, thousands of times. |app.DynamicsSource| is that seam:
% it implements gjToolbox's frozen |DynamicsSource| interface, and it has two evaluation paths, the
% normative MATLAB one (|Model.accel| at every call) and a NATIVE one in which everything that
% depends on time only is tabulated once and the whole force chain is evaluated in C. This
% tutorial builds both, shows what the native path refuses, checks the two against each other and
% against |Model.accel|, and propagates a five-minute arc with both bindings.
%
% *What you will learn*
%
% * what |app.DynamicsSource| needs (a complete Model, a TAI run epoch) and its three-phase
%   life: construct, |prepare| per run, |accel| per step, |release|;
% * which options the native path refuses by name (dtm2020, ceresClimatology, transient walls
%   without 'ThermalOrbit' or on the MATLAB path) and why;
% * how 'ThermalOrbit' serves the standard preset's transient wall temperatures as a table;
% * how close the native evaluator is to |Model.accel| at a table node and between nodes;
% * how to propagate with |gjToolbox.propagate| on the "vectorized" and the "native" binding and
%   compare the accelerations along the propagated arc with |Model.accel|.
%
% *Requirements:* none - every file read is committed in |data/|. gjToolbox must be installed
% (|gjToolboxSetup|) with its native stepping kernel built (|buildtool mex| in odeSuite), and
% NonGravForces' MEX target "rhs" must be built (|NonGravForces.buildMex|); the gate's |ciGate|
% builds both. Run time: under a minute.
%
% Requirements: none
%
% NonGravForces v0.20 | 12-Nov-2026 | Claude Opus 5.5

%% Step 1 - Check the integrator and the native evaluator
% |app.DynamicsSource| subclasses |gjToolbox.DynamicsSource|, so the class does not even load
% without gjToolbox. The native path needs the compiled evaluator |mex/ngfRhs.c|, the
% |NonGravForces.buildMex| target "rhs"; |'Check', true| only reports whether it is built and
% current, it compiles nothing. The "native" binding of the integrator also needs gjToolbox's own
% compiled stepping kernel.
assert(~isempty(which('gjToolbox.DynamicsSource')), 'gjToolbox is not on the path: run gjToolboxSetup from the odeSuite checkout');
rep = NonGravForces.buildMex('Check', true, 'Target', "rhs");
fprintf('rhs evaluator %s: up to date %d\n', rep.file, rep.upToDate);
fprintf('gjToolbox native stepping kernel available: %d\n', gjToolbox.nativeStepper("available"));
assert(rep.upToDate, 'build the native evaluator first: NonGravForces.buildMex(''Target'', "rhs")');

%% Step 2 - The data and the run epoch
% The reference day's ten minutes: the GNV1B orbit (rotated into the GCRS), the SCA1B attitude,
% the space weather and the EOP. An integrator counts time in SECONDS FROM A RUN EPOCH; gjToolbox
% takes that epoch as an unzoned datetime in TAI, and the source must be given the same one. Here
% it is the first GNV1B epoch, so t = 0 is 2020-01-01 00:00:00 GPS and t = k - 1 is record k.
% |orbAll| packs the whole excerpt as a struct (.epoch, .r, .v) - the form 'ThermalOrbit' takes.
root = fileparts(which('setupNonGravForces'));
f = @(varargin) fullfile(root, 'data', varargin{:});
eop = referenceSystem.readEop(f('eop', 'finals2000A_2019-09_2020-03.all'));
sw = NonGravForces.io.spaceWeather(f('spaceweather', 'SW-All_2019-09_2020-03.csv'));
cls = NonGravForces.io.spaceWeather(f('spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
swarm = string(f('samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
[ep, r, v] = NonGravForces.io.readGNV1B(f('samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
[eq, q] = NonGravForces.io.readSCA1B(f('samples', 'SCA1B_2020-01-01_C_04.txt'));
att = NonGravForces.geom.AttitudeProvider(eq, q);
macro = NonGravForces.geom.graceFO("C");
e1 = referenceSystem.daySecondToEpoch([ep.day(1) ep.sec(1)], ep.scale);
tai = referenceSystem.epochTo(referenceSystem.convertTimeScale(e1, "TAI"), "datetime");
orbAll = struct('epoch', referenceSystem.daySecondToEpoch([ep.day(:) ep.sec(:)], ep.scale), 'r', r, 'v', v);
fprintf('run epoch (TAI) %s; %d orbit records, t = 0 .. %d s\n', string(tai), numel(ep.sec), round(ep.sec(end) - ep.sec(1)));

%% Step 3 - What the native path refuses
% The native evaluator ports the whole chain to C, and an option it does not port is REFUSED by
% name rather than evaluated as something else. Four cases, each caught and its identifier
% printed:
%
% * |drag.density| "dtm2020" and |erp.backend| "ceresClimatology" have no native evaluator
%   (|unsupportedOption|, at construction) - use them through |Model.accel| or |series|;
% * the standard preset's TRANSIENT wall temperatures depend on their history, while an
%   integrator calls |accel| at times out of order and repeatedly at one t. The native path
%   therefore integrates the walls ONCE along a reference orbit, 'ThermalOrbit', and reads them by
%   time; without it, or on the MATLAB path ('Native' false), where every call would restart the
%   walls from the steady state, construction is refused (|inputModel|).
%
% The first two Models use the standard preset unscaled with steady walls, so that the refusal
% shown is the one about the option named.
base = NonGravForces.Options("standard", "drag.densityScale", "none", "thermal.model", "steadyState");
S = NonGravForces.env.densityScale(swarm, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));
Mstd = NonGravForces.Model(macro, NonGravForces.Options("standard"), 'Attitude', att, 'Eop', eop, ...
    'SpaceWeather', sw, 'DensityScale', S);
cases = {
    "dtm2020",             NonGravForces.Model(macro, base.set("drag.density", "dtm2020"), 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw), {}
    "ceresClimatology",    NonGravForces.Model(macro, base.set("erp.backend", "ceresClimatology"), 'Attitude', att, 'Eop', eop, 'SpaceWeather', sw), {}
    "transient, no orbit", Mstd, {}
    "transient, MATLAB",   Mstd, {'Native', false, 'ThermalOrbit', orbAll}};
for c = 1:size(cases, 1)
    try
        NonGravForces.app.DynamicsSource(cases{c, 2}, 'Epoch', tai, cases{c, 3}{:});
        fprintf('%-22s accepted\n', cases{c, 1});
    catch err
        fprintf('%-22s refused: %s\n', cases{c, 1}, err.identifier);
    end
end

%% Step 4 - The standard preset on the native path, with 'ThermalOrbit'
% The standard preset as shipped, with the reference orbit for its walls. |prepare(context)| binds
% the source to one run: the context's 'Span' (and the integrator's 'StencilSpan') in seconds from
% the run epoch, plus 'Margin_s' on either side, is what the tables cover - the attitude and the
% thermal orbit must cover it too. 'TableStep_s' is the table spacing. The tables are exact at their
% nodes (every multiple of 60 s here), so at a node the native acceleration and |Model.accel| are
% two implementations of the same expressions. |Model.accel| is evaluated along the WHOLE orbit,
% as the table was: the transient walls start from the steady state at the first epoch of a call.
src = NonGravForces.app.DynamicsSource(Mstd, 'Epoch', tai, 'Margin_s', 30, 'TableStep_s', 60, 'ThermalOrbit', orbAll);
src = src.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
aG = Mstd.accel(orbAll.epoch, orbAll.r, orbAll.v, [], 'Frame', "gcrs");
fprintf('%6s %10s %14s %12s\n', 't [s]', 'node', '|a| [m/s^2]', 'rel. diff');
for k = [121 151 181 241 301 361 391 421]
    t = double(ep.sec(k) - ep.sec(1));
    a1 = src.accel(t, r(:, k), v(:, k));
    fprintf('%6d %10d %14.4e %12.2e\n', t, mod(t, 60) == 0, norm(aG(:, k)), norm(a1 - aG(:, k)) / norm(aG(:, k)));
end
d = src.diagnostics(180, r(:, 181), v(:, 181));
fprintf('native diagnostics at t = 180 s: rho %.3e kg/m^3 (scale %.3f), shadow factor %.2f, mass %.3f kg\n', ...
    d.rho, d.densityScale, d.nu, d.mass_kg);
src.release();

%% Step 5 - The two paths of one source
% With steady-state walls ('thermal.model' "steadyState") the MATLAB path can serve the same Model:
% 'Native' false makes |accel| call |Model.accel| itself. The two paths of one source agree at a
% node to rounding; the native one exists to be cheap in an integrator's inner loop, which the
% timing of single-state calls shows.
Ms = NonGravForces.Model(macro, NonGravForces.Options("standard").set("thermal.model", "steadyState"), 'Attitude', att, ...
    'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
srcN = NonGravForces.app.DynamicsSource(Ms, 'Epoch', tai, 'Margin_s', 30, 'TableStep_s', 60);
srcN = srcN.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
srcM = NonGravForces.app.DynamicsSource(Ms, 'Epoch', tai, 'Margin_s', 30, 'TableStep_s', 60, 'Native', false);
srcM = srcM.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
aN = srcN.accel(180, r(:, 181), v(:, 181));
aM = srcM.accel(180, r(:, 181), v(:, 181));
fprintf('native vs MATLAB path at t = 180 s: relative difference %.2e\n', norm(aN - aM) / norm(aM));
n = 20;
t0 = tic; for k = 1:n, srcN.accel(180, r(:, 181), v(:, 181)); end, tN = toc(t0) / n;
t0 = tic; for k = 1:n, srcM.accel(180, r(:, 181), v(:, 181)); end, tM = toc(t0) / n;
fprintf('one state: native %.3f ms, MATLAB %.3f ms (ratio %.3f, this host)\n', 1e3 * tN, 1e3 * tM, tN / tM);
srcN.release(); srcM.release();

%% Step 6 - Propagate five minutes with both bindings
% |gjToolbox.Dynamics(source, 'Binding', b)| tells the integrator how to call the source:
% "vectorized" calls |accel| from MATLAB, "native" hands the source's compiled entry to the C
% stepper, which then steps AND evaluates in C. |gjToolbox.Options| sets the run epoch (the same
% TAI datetime), the step (s), the order and the start-up.
% |gjToolbox.propagate(dyn, [t0 t1], r0, v0, opt)| returns a Solution whose |T| (s from the
% epoch), |R| and |V| (N x 3, GCRS) are the nodes.
%
% The arc has NO GRAVITY: this toolbox computes none, and a gravity source is another toolbox's
% (the 'Gravity' option). The state therefore flies a straight line and climbs tens of kilometres
% in five minutes - out of the Swarm height bracket the observed density scale is defined over,
% where the standard scale refuses. So this run uses the DRIVER scale (one factor per day, no
% height dependence) and steady walls, as the toolbox's own integration test does. The point is
% that the source is the same on both bindings, not that the arc is an orbit.
Ma = NonGravForces.Model(macro, NonGravForces.Options("standard").set("drag.densityScale", "driver").set("thermal.model", "steadyState"), ...
    'Attitude', att, 'Eop', eop, 'SpaceWeather', sw, 'DensityScale', NonGravForces.env.densityDriverScale(sw, cls));
warnState = warning('off', 'gjToolbox:dynamics:fdGradient');
opt = gjToolbox.Options('Epoch', tai, 'Step', 5, 'Order', 8, 'Startup', "onesided");
srcV = NonGravForces.app.DynamicsSource(Ma, 'Epoch', tai, 'Margin_s', 30, 'Native', false);
srcC = NonGravForces.app.DynamicsSource(Ma, 'Epoch', tai, 'Margin_s', 30);
dynV = gjToolbox.Dynamics(srcV, 'Binding', "vectorized", 'ValidateGradient', false, 'FiniteDifferenceFallback', false);
dynC = gjToolbox.Dynamics(srcC, 'Binding', "native", 'ValidateGradient', false, 'FiniteDifferenceFallback', false);
r0 = r(:, 121); v0 = v(:, 121);
t1 = tic; solV = gjToolbox.propagate(dynV, [120 420], r0, v0, opt); secV = toc(t1);
t1 = tic; solC = gjToolbox.propagate(dynC, [120 420], r0, v0, opt); secC = toc(t1);
warning(warnState);
fprintf('vectorized: %d nodes, %d force evaluations, %.2f s\n', numel(solV.T), solV.Diagnostics.ForceEvals.total, secV);
fprintf('native:     %d nodes, %d force evaluations, %.2f s\n', numel(solC.T), solC.Diagnostics.ForceEvals.total, secC);
fprintf('final positions differ by %.3e m\n', norm(solC.R(end, :) - solV.R(end, :)));
fprintf('height above the start: %+.1f km after %d s (no gravity)\n', (norm(solC.R(end, :)) - norm(r0)) / 1e3, round(solC.T(end) - solC.T(1)));

%% Step 7 - The accelerations along the propagated arc
% The last check closes the loop: at every node of the propagated arc, the native source's
% acceleration against |Model.accel| at the same epoch and state (|Frame| "gcrs"). The epoch of
% node time t is the first GNV1B epoch plus t seconds, in GPS. At the 60 s table nodes the two
% agree to rounding; between them what is left is the cubic interpolation of the time tables.
srcP = NonGravForces.app.DynamicsSource(Ma, 'Epoch', tai, 'Margin_s', 30);
srcP = srcP.prepare(struct('Span', [120 420], 'StencilSpan', [120 420], 'MaxColumns', 1));
T = solC.T(:); R = solC.R.'; V = solC.V.';
eT = referenceSystem.daySecondToEpoch([repmat(ep.day(1), numel(T), 1), ep.sec(1) + T], ep.scale);
aModel = Ma.accel(eT, R, V, [], 'Frame', "gcrs");
aSrc = zeros(3, numel(T));
for k = 1:numel(T)
    aSrc(:, k) = srcP.accel(T(k), R(:, k), V(:, k));
end
srcP.release();
rel = vecnorm(aSrc - aModel) ./ vecnorm(aModel);
node = mod(round(T.'), 60) == 0;
fprintf('native vs Model.accel along the arc: worst %.2e at table nodes, %.2e between them\n', max(rel(node)), max(rel(~node)));
figure
subplot(2, 1, 1)
plot(T, 1e9 * vecnorm(aModel), T, 1e9 * vecnorm(aSrc), '--'), grid on
ylabel('|a| (nm/s^2)'), legend(["Model.accel" "native source"]), title('along the propagated arc')
subplot(2, 1, 2)
semilogy(T, max(rel, eps), '.-'), grid on
xlabel('t (s from the run epoch)'), ylabel('relative difference')
