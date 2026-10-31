classdef DynamicsSource < gjToolbox.DynamicsSource
%DYNAMICSSOURCE The non-gravitational acceleration as gjToolbox's frozen integrator interface, with a native C entry.
%
%   src = NonGravForces.app.DynamicsSource(model, 'Epoch', tai)
%   src = NonGravForces.app.DynamicsSource(model, 'Epoch', tai, 'Gravity', gravitySource)
%   src = NonGravForces.app.DynamicsSource(model, 'Epoch', tai, 'Parameters', ["cdScale" "crScale"])
%   src = NonGravForces.app.DynamicsSource(model, 'Epoch', tai, 'ThermalOrbit', orb)   % thermal.model transient
%   sol = gjToolbox.propagate(src, [0 5400], r0, v0, gjToolbox.Options('Epoch', tai, 'Step', 5, 'Binding', "native"))
%
%   WHAT IT IS. NonGravForces computes forces; gjToolbox integrates. This class is
%   the seam, and it is the interface gjToolbox froze (gjToolbox.DynamicsSource):
%   a three-phase lifecycle (construct once, prepare per run, evaluate per step),
%   state 3 (GCRS), t in seconds from the run epoch, no parameter values crossing
%   the boundary. Model.rhs, whose handle matched no signature the integrator
%   takes, is gone (HANDOVER O-4, issue #35).
%
%   TWO EVALUATION PATHS, ONE MODEL. With 'Native' true (the default) prepare()
%   builds every time-only quantity into tables -- the ITRS->GCRS rotation and its
%   rate, Sun and Moon, the attitude nodes, the space-weather drivers, the UTC
%   calendar, the mass, the CERES months, the self-shadow table and, under
%   geom.flowShadow "table", geom.flowShadow's table of the gas-stream ratios -- hands them to
%   the compiled evaluator (mex/ngfRhs.c, NonGravForces.buildMex target "rhs") and
%   publishes a gj_native_entry, so a run bound "native" steps in C and evaluates
%   in C. accel() then calls THAT entry, so the MATLAB side of such a run
%   evaluates one model and not two. With 'Native' false accel() calls
%   Model.accel, which stays normative: the parity of the two is asserted in
%   tests/tNgfIntegration.m in every preset (for the standard preset, whose
%   transient walls 'Native' false cannot serve, through its steady-state
%   variant, and the transient table separately against Model.accel along the
%   reference orbit).
%
%   GRAVITY IS NOT HERE and cannot be: this toolbox computes non-gravitational
%   forces only (P9). A run needs the gravity field of another toolbox, and the
%   native stepper binds ONE entry, so 'Gravity' takes another
%   gjToolbox.DynamicsSource: it is prepared first, its acceleration is added
%   before this one's, and if it publishes a native entry the C evaluator calls
%   that entry itself.
%
%   WHAT THE NATIVE PATH SERVES. Every option Model evaluates except
%   geom.selfShadow "raycast" (use "lut"), an attitude provider built with
%   "squad" (prepare() refuses both), and erp.backend "ceresClimatology" and
%   drag.density "dtm2020", which have no native evaluator and are refused at
%   construction (HANDOVER O-44, O-52). Each is refused by name rather than
%   degraded to something else. thermal.model "transient" (2026-10-30, the owner's ruling
%   that the standard preset goes transient and stays integrable) is served by
%   a TABLE: the walls' temperatures depend on their history, and the integrator
%   calls accel non-monotonically and repeatedly at one t, so prepare()
%   integrates them ONCE, with Model.accel along a reference orbit
%   ('ThermalOrbit', e.g. the arc's GNV1B), and the C evaluator reads them by
%   time like every other table and forms the recoil from them. Only the native
%   path can do that: with 'Native' false Model.accel would restart the walls
%   from the steady state at every call, so that combination is refused. drag.densityScale
%   (swarm, driver, combined) is served since 2026-10-29: prepare() turns the
%   Model's scale into one row per UTC day -- Swarm's bracket, interpolated in
%   height in C as env.densityScale does, or the day's single factor -- and it
%   switches at the aligned driver-bin edge of UTC midnight, which is already a
%   declared breakpoint. A UTC day with no scale (drag.densityScale "swarm" on a
%   day Swarm does not cover), or a height further outside Swarm's bracket than
%   the bracket is wide, is refused by the native evaluator when it is reached
%   (NonGravForces:rhsGateway:dataRange), not at prepare().
%
%   THE TABLES ARE EXACT AT THEIR NODES. The rotation is a cubic Hermite of R with
%   Rdot, everything else a cubic Lagrange, on a 'TableStep_s' grid (default 60 s,
%   about 1e-12 relative on the rotation between nodes). The drivers are piecewise
%   constant over the 3-hour UTC bins the space-weather file resolves, and the
%   source declares those switches as Breakpoints, aligned to the run's grid so
%   that the declaration and the model change at the same instant.
%
% CONSTRUCTOR ARGUMENTS
%   Name     Size  Type                  Units  Description
%   ----     ----  ----                  -----  -----------
%   model    1x1   NonGravForces.Model   -      Configured and complete (attitude, EOP, space weather,
%                                               CERES, mass) for the arc to be integrated
%
% NAME-VALUE ARGUMENTS
%   Name          Size  Type      Units  Description
%   ----          ----  ----      -----  -----------
%   Epoch         1x1   datetime  TAI    The run epoch, unzoned, as gjToolbox.Options carries it.
%                                        Required (default NaT, refused); a context Epoch must equal it.
%   Attitude      1x1   string    -      "provider" (default, the Model's) or "nominal"
%   AlongVelocity 1x1   logical   -      Nominal attitude: +X along the velocity (default true)
%   Parameters    1xP   string    -      Estimated parameters: any of "cdScale", "crScale"; default none
%   Native        1x1   logical   -      Publish the compiled entry (default true)
%   Gravity       1x1   object    -      Another gjToolbox.DynamicsSource added before this one; default []
%   Margin_s      1x1   double    s      Table margin beyond the stencil span (default 3600)
%   TableStep_s   1x1   double    s      Table spacing (default 60)
%   ThermalOrbit  1x1   struct    -      thermal.model transient only, and required by it: the reference
%                                        orbit (.epoch, a referenceSystem epoch struct with K increasing
%                                        entries; .r and .v, 3xK in GCRS, m and m/s; e.g. the three outputs
%                                        of io.readGNV1B packed into a struct) along which
%                                        prepare() integrates the wall temperatures; it must cover the
%                                        tables' span and should start well before it (the walls start
%                                        from the steady state at its first epoch). Default struct()
%
% PROPERTIES (beyond the interface's)
%   model, epoch, attitudeSource, alongVelocity, parameters, native, gravity, margin_s,
%   tableStep_s, thermalOrbit, tables (the prepared tables: .nT, .nb, .breakpoints, .handle,
%   .tPanel, .scale, ...), StateSize (3), ParameterRule, Vectorized (true)
%
% METHODS
%   prepare(context)   Bind to one run: tables, the native entry, the breakpoints, the ParameterSpec
%   accel(t, R, V)     3xm accelerations in GCRS, from the entry (native) or from Model.accel
%   hasGradient()      true only when Parameters were declared
%   gradient(t, r, v)  Gr, Gv by central differences of this source; Gp analytic from Model.partials
%   release()          Free the native context (called under onCleanup by gjToolbox.Dynamics)
%
% ACCURACY
%   The native path agrees with Model.accel to the rounding of two implementations
%   of the same expressions at a table node, and to the tables' interpolation
%   error between nodes; both are measured per run in tests/tNgfIntegration.m
%   (2.5e-14 and 8.1e-11 relative at a 60 s table step, measured 2026-09-16).
%   The time tables are exact at their nodes by construction.
%   thermal.model transient: the temperatures are those of the REFERENCE orbit, not
%   of the integrated state, so the thermal recoil ignores the state's departure from
%   it (a few km of along-track error moves the eclipse edges by under a second); at a
%   node and on the reference orbit the native path equals Model.accel of the same orbit
%   (tNgfIntegration/DynamicsSourceThermalTableReference). Between nodes the cubic
%   reads the temperatures as smooth, which they are not at an eclipse edge (the
%   radiator's time constant is ~44 s): shrink 'TableStep_s' for an arc that needs it.
%
% METHODS BEYOND THE gjToolbox INTERFACE
%   d  = diagnostics(obj, t, r, v)  What the NATIVE evaluator computed on the way to one
%                                   acceleration: .drag .radiation .thermal (3x1 each, m/s^2),
%                                   .rho, .species (8x1), .Tgas_K, .vrel_body (3x1), .lla (3x1),
%                                   .nu, .mass_kg, .rotationBodyToGcrs (3x3), .densityScale (1
%                                   without drag.densityScale; .rho and .species are scaled). It exists to pin
%                                   WHERE the native path and Model.accel differ rather than only
%                                   that they do, which is what the parity tests need when a bound
%                                   moves. Requires 'Native', true and a prepared source; raises
%                                   inputModel otherwise.
%   ep = epochAt(obj, t, m)         The referenceSystem epoch (TAI) of the integrator offsets t,
%                                   repeated m times, so a caller can label a propagated arc with
%                                   the same time scale the tables were built on. t is seconds
%                                   from the source's 'Epoch'.
%
% ERRORS
%   NonGravForces:DynamicsSource:inputModel     model is not a scalar NonGravForces.Model; its options
%                                               cannot be served (raycast; a non-slerp attitude provider;
%                                               transient thermal without 'ThermalOrbit' or with 'Native'
%                                               false); the Model lacks what its options need (Eop, the
%                                               attitude provider, SpaceWeather, Ceres, DensityScale,
%                                               Exponential); 'Parameters' outside cdScale/crScale;
%                                               'Margin_s' < 0 or 'TableStep_s' <= 0; accel() given a
%                                               non-scalar t or mismatched R, V; diagnostics() without Native
%   NonGravForces:DynamicsSource:unsupportedOption  erp.backend ceresClimatology or drag.density dtm2020,
%                                               which have no native evaluator; use Model.accel or series instead
%   NonGravForces:DynamicsSource:inputEpoch     no 'Epoch' (or a zoned one), or a context Epoch that disagrees
%   NonGravForces:DynamicsSource:inputGravity   'Gravity' is not a gjToolbox.DynamicsSource, or (Native)
%                                               publishes no NativeEntry
%   NonGravForces:DynamicsSource:depMissing     the rhs MEX is not built (Native), or sun.source rs.de440
%                                               without the Model's 'Ephemeris'. Without gjToolbox the
%                                               class does not load at all (MATLAB's own error)
%   NonGravForces:DynamicsSource:dataRange      the attitude nodes or 'ThermalOrbit' do not cover the
%                                               tables' span (span plus the margin); a space-weather file
%                                               that does not cover it raises io.spaceWeather's own error
%   NonGravForces:rhsGateway:dataRange          native evaluation outside its domain: t outside the tables,
%                                               a UTC day without a density scale, a height too far
%                                               outside Swarm's bracket, hwm14 below the ellipsoid
%   NonGravForces:rhsGateway:numericModel       native evaluation failed otherwise (e.g. the gravity entry)
%   gjToolbox:dynamics:notPrepared              accel() before prepare()
%
% REFERENCES
%   gjToolbox.DynamicsSource (the frozen interface) and csrc/gj_native.h (ABI 1);
%   HANDOVER.md Part G LP-C, WP-C9.
%
% EXAMPLE
%   % M: a standard-preset Model (thermal.model transient, so 'ThermalOrbit' is required;
%   % built with 'DensityScale', as the preset requires) whose attitude and space weather
%   % cover 2020-01-01; the day's GNV1B is the reference orbit. The run starts at 02:00 UTC
%   % so that the tables' default 3600 s margin lies inside the day.
%   [ep, r, v] = NonGravForces.io.readGNV1B(gnvFile, 'Eop', eop);
%   orb = struct('epoch', ep, 'r', r, 'v', v);
%   src = NonGravForces.app.DynamicsSource(M, 'Epoch', datetime(2020, 1, 1, 2, 0, 37), 'ThermalOrbit', orb);
%   src = src.prepare(struct('Span', [0 600], 'StencilSpan', [0 600], 'MaxColumns', 1));
%   k = find(ep.sec >= 7218, 1);                    % 02:00 UTC is 02:00:18 GPS
%   a = src.accel(0, r(:, k), v(:, k));
%
% See also NonGravForces.Model, gjToolbox.DynamicsSource, NonGravForces.buildMex.
%
% KEYWORDS integrator dynamics source native entry gauss-jackson orbit integration
% Created 2026-09-16 by Claude Opus 5 (Anthropic).
    properties (SetAccess = protected)
        StateSize = 3
        ParameterRule = "none"
    end
    properties (SetAccess = private)
        model
        epoch datetime = NaT
        attitudeSource (1,1) string = "provider"
        alongVelocity (1,1) logical = true
        parameters string = strings(1, 0)
        native (1,1) logical = true
        gravity = []
        margin_s (1,1) double = 3600
        tableStep_s (1,1) double = 60
        thermalOrbit struct = struct()
        Vectorized (1,1) logical = true
        tables struct = struct('prepared', false)
    end
    methods
        function obj = DynamicsSource(model, nv)
            arguments
                model
                nv.Epoch datetime = NaT
                nv.Attitude (1,1) string {mustBeMember(nv.Attitude, ["provider","nominal"])} = "provider"
                nv.AlongVelocity (1,1) logical = true
                nv.Parameters string = strings(1, 0)
                nv.Native (1,1) logical = true
                nv.Gravity = []
                nv.Margin_s (1,1) double = 3600
                nv.TableStep_s (1,1) double = 60
                nv.ThermalOrbit struct = struct()
            end
            unit = "DynamicsSource";
            if ~isa(model, 'NonGravForces.Model') || ~isscalar(model)
                NonGravForces.fail(unit, "inputModel", "model must be a scalar NonGravForces.Model -> build one with NonGravForces.Model(macro, options, ...)");
            end
            o = model.options;
            if o.thermal.model == "transient"
                % The walls carry their temperatures from epoch to epoch, and the integrator calls accel
                % non-monotonically and repeatedly at one t: only a table integrated once along a reference
                % orbit can serve that, and only the native path reads one.
                if ~nv.Native
                    NonGravForces.fail(unit, "inputModel", "thermal.model transient is served by the native evaluator's temperature table only (Model.accel would restart the walls from the steady state at every call) -> construct with 'Native', true and 'ThermalOrbit'");
                end
                if ~(isscalar(nv.ThermalOrbit) && all(isfield(nv.ThermalOrbit, {'epoch', 'r', 'v'})))
                    NonGravForces.fail(unit, "inputModel", "thermal.model transient needs 'ThermalOrbit' (.epoch, .r, .v in GCRS, e.g. io.readGNV1B), along which prepare integrates the wall temperatures once -> pass the arc's reference orbit, starting well before it, or use thermal.model steadyState");
                end
            end
            % REFUSED AT CONSTRUCTION, not when the tables are packed: an option this evaluator
            % cannot serve should cost the caller nothing to discover. It was first written into
            % the packing path, where it fired only after prepare() and a contract test that
            % expected it at construction reported "nothing was thrown" -- the guard existed and
            % the test was still right to fail (HANDOVER O-44).
            %
            % THE REASON THIS ROW USED TO GIVE WAS FALSE, and it is worth not repeating: it said a
            % 64800-cell grid "is not table material". The native path ALREADY carries three EBAF
            % grids of exactly that size, plus their axes and month nodes, for erp.backend ceres --
            % 3.1 MB against the 2.1 MB this backend would need. Memory was never the obstacle.
            % The real reason is that it would buy no physics: EBAF is also 180x360 at one degree
            % and also a measured exitance, so over an integration span the two differ in their
            % SOURCE and their reference level, not in what they can express. What ceresClimatology
            % uniquely offers is an epoch range of 1900-2101 and provenance shared with GROOPS, and
            % both are post-processing concerns. An integrator uses ceres where the months are in
            % hand and knocke where they are not (measured 2026-09-18).
            if o.erp.backend == "ceresClimatology"
                NonGravForces.fail(unit, "unsupportedOption", "erp.backend ceresClimatology has no native evaluator, by design: erp.backend ceres is the same 1-degree measured field and the native path serves it -> use ceres for an integrated arc, knocke if you do not hold the months, or ceresClimatology through Model.accel or series");
            end
            % DTM2020 (2026-10-27) has no native evaluator YET: mex/ngfRhs.c links the density models it
            % evaluates, and CNES's subroutines were not added to it. Refused here, for the same reason as
            % ceresClimatology above: find() below would return empty and pack 0, "exponential".
            if o.drag.density == "dtm2020"
                NonGravForces.fail(unit, "unsupportedOption", "drag.density dtm2020 has no native evaluator (HANDOVER O-52) -> use msis00 or msis21 for an integrated arc, or dtm2020 through Model.accel or series");
            end
            if any(~ismember(nv.Parameters, ["cdScale", "crScale"]))
                NonGravForces.fail(unit, "inputModel", "'Parameters' may name cdScale and crScale -> got %s", strjoin(nv.Parameters, ", "));
            end
            if ~isempty(nv.Gravity) && ~isa(nv.Gravity, 'gjToolbox.DynamicsSource')
                NonGravForces.fail(unit, "inputGravity", "'Gravity' must be a gjToolbox.DynamicsSource -> the gravity field belongs to its own toolbox (P9)");
            end
            if isnat(nv.Epoch) || ~isempty(nv.Epoch.TimeZone)
                NonGravForces.fail(unit, "inputEpoch", "'Epoch' must be an unzoned TAI datetime -> the same one gjToolbox.Options carries; it is never assumed");
            end
            if ~(nv.Margin_s >= 0) || ~(nv.TableStep_s > 0)
                NonGravForces.fail(unit, "inputModel", "'Margin_s' must be >= 0 and 'TableStep_s' > 0 -> got %g and %g", nv.Margin_s, nv.TableStep_s);
            end
            obj.model = model;
            obj.epoch = nv.Epoch;
            obj.attitudeSource = nv.Attitude;
            obj.alongVelocity = nv.AlongVelocity;
            obj.parameters = string(nv.Parameters(:))';
            obj.native = nv.Native;
            obj.gravity = nv.Gravity;
            obj.margin_s = nv.Margin_s;
            obj.tableStep_s = nv.TableStep_s;
            obj.thermalOrbit = nv.ThermalOrbit;
            obj.ParameterRule = strjoin(["fixed", obj.parameters], " ");
        end

        function obj = prepare(obj, context)
            %PREPARE Bind to one run: the tables, the native entry, the breakpoints and the spec.
            arguments
                obj (1,1) NonGravForces.app.DynamicsSource
                context struct = struct()
            end
            unit = "DynamicsSource";
            if isfield(context, 'Epoch') && isa(context.Epoch, 'datetime') && ~isnat(context.Epoch) && context.Epoch ~= obj.epoch
                NonGravForces.fail(unit, "inputEpoch", "the run's epoch %s is not this source's %s -> construct the source with the run's epoch", string(context.Epoch), string(obj.epoch));
            end
            span = [0 0];
            if isfield(context, 'StencilSpan') && numel(context.StencilSpan) == 2
                span = double(context.StencilSpan);
            elseif isfield(context, 'Span') && numel(context.Span) == 2
                span = double(context.Span);
            end
            span = [min(span) max(span)] + [-obj.margin_s, obj.margin_s];
            maxCols = 1;
            if isfield(context, 'MaxColumns') && ~isempty(context.MaxColumns)
                maxCols = max(1, double(context.MaxColumns));
            end
            grid = [];
            if isfield(context, 'Grid') && ~isempty(context.Grid)
                grid = context.Grid;
            end
            if ~isempty(obj.gravity)
                obj.gravity = obj.gravity.prepare(context);
            end
            t = obj.buildTables(span, maxCols, grid);
            obj.tables = t;
            obj.Breakpoints = t.breakpoints;
            obj.BreakpointKinds = repmat("restart", 1, numel(t.breakpoints));
            p = numel(obj.parameters);
            obj.ParameterSpec = struct('Count', p, 'ActiveFrom', -inf(1, p), 'ActiveTo', inf(1, p), ...
                'Id', "NonGravForces:" + strjoin([obj.parameters, string(obj.model.options.preset)], ","));
            obj.NativeEntry = t.entry;
            obj.IsPrepared = true;
        end

        function A = accel(obj, t, R, V)
            %ACCEL The GCRS acceleration of an n-by-m block at t, from the entry or from Model.accel.
            obj.requirePrepared();
            if ~(isscalar(t) && isfinite(t)) || size(R, 1) ~= 3 || ~isequal(size(R), size(V))
                NonGravForces.fail("DynamicsSource", "inputModel", "accel takes a scalar t and 3xm R and V of the same size -> got %s and %s", mat2str(size(R)), mat2str(size(V)));
            end
            if obj.native
                A = rhsGateway('accel', obj.tables.handle, t, R, V);
                return
            end
            m = size(R, 2);
            ep = obj.epochAt(t, m);
            q = [];
            if obj.attitudeSource == "nominal"
                q = NonGravForces.geom.AttitudeProvider.nominal(R, V, obj.alongVelocity);
            end
            A = obj.model.accel(ep, R, V, q, 'Frame', "gcrs");
            if ~isempty(obj.gravity)
                A = A + obj.gravity.accel(t, R, V);
            end
        end

        function tf = hasGradient(obj)
            %HASGRADIENT Only with estimated parameters: without them the integrator's own fallback is honest.
            tf = ~isempty(obj.parameters);
        end

        function [Gr, Gv, Gp] = gradient(obj, t, r, v)
            %GRADIENT Central differences of THIS source for Gr and Gv; Model's analytic parameter partials for Gp.
            obj.requirePrepared();
            Gr = obj.stateJacobian(t, r, v, 1);
            Gv = obj.stateJacobian(t, r, v, 2);
            p = numel(obj.parameters);
            Gp = zeros(3, p);
            if p > 0
                ep = obj.epochAt(t, 1);
                q = [];
                if obj.attitudeSource == "nominal"
                    q = NonGravForces.geom.AttitudeProvider.nominal(r, v, obj.alongVelocity);
                end
                J = obj.model.partials(ep, r, v, q, 'WithRespectTo', obj.parameters, 'Frame', "gcrs");
                Gp = reshape(J, 3, p);
            end
        end

        function obj = release(obj)
            %RELEASE Free the native context; the entry must not outlive it.
            if isfield(obj.tables, 'handle') && ~isempty(obj.tables.handle)
                try
                    rhsGateway('destroy', obj.tables.handle);
                catch
                    % release runs under onCleanup and must not throw
                end
            end
            obj.tables = struct('prepared', false);
            obj.NativeEntry = [];
            if ~isempty(obj.gravity)
                obj.gravity = obj.gravity.release();
            end
            obj.IsPrepared = false;
        end

        function d = diagnostics(obj, t, r, v)
            %DIAGNOSTICS What the native evaluator computed on the way to one acceleration.
            %   For diagnosing a parity difference against Model.accel's info, and for the tests that
            %   pin WHERE two implementations differ rather than only that they do.
            obj.requirePrepared();
            if ~obj.native
                NonGravForces.fail("DynamicsSource", "inputModel", "diagnostics reads the native evaluator -> construct with 'Native', true");
            end
            x = rhsGateway('diagnose', obj.tables.handle, t, r(:), v(:));
            d = struct('drag', x(1:3), 'radiation', x(4:6), 'thermal', x(7:9), 'rho', x(10), ...
                'species', x(11:18), 'Tgas_K', x(19), 'vrel_body', x(20:22), 'lla', x(23:25), ...
                'nu', x(26), 'mass_kg', x(27), 'rotationBodyToGcrs', reshape(x(28:36), 3, 3), 'densityScale', x(37));
        end

        function ep = epochAt(obj, t, m)
            %EPOCHAT The referenceSystem epoch (TAI) of the offsets t, repeated m times.
            e0 = referenceSystem.datetimeToEpoch(obj.epoch, "TAI");
            ep = referenceSystem.daySecondToEpoch([repmat(double(e0.day), m, 1), double(e0.sec) + t(:) .* ones(m, 1)], "TAI");
        end
    end

    methods (Access = private)
        function G = stateJacobian(obj, t, r, v, which)
            % Central differences of this source's own acceleration, which is what makes the
            % gradient honest whichever path accel() takes.
            base = (which == 1) * norm(r) + (which == 2) * norm(v);
            h = eps^(1/3) * max(base, 1);
            G = zeros(3, 3);
            for i = 1:3
                d = zeros(3, 1); d(i) = h;
                if which == 1
                    G(:, i) = (obj.accel(t, r + d, v) - obj.accel(t, r - d, v)) / (2 * h);
                else
                    G(:, i) = (obj.accel(t, r, v + d) - obj.accel(t, r, v - d)) / (2 * h);
                end
            end
        end

        function t = buildTables(obj, span, maxCols, grid)
            %BUILDTABLES Every time-only quantity of the run, from the same units Model.accel calls.
            unit = "DynamicsSource";
            M = obj.model;
            o = M.options;
            if isempty(M.eop)
                NonGravForces.fail(unit, "inputModel", "the Model needs 'Eop' -> the rotation and the co-rotation come from it");
            end
            dt = obj.tableStep_s;
            t0 = floor(span(1) / dt) * dt - dt;
            nT = max(4, ceil((span(2) - t0) / dt) + 1);   % one node past the span is what a cubic needs
            tTab = t0 + (0:nT-1) * dt;
            e0 = referenceSystem.datetimeToEpoch(obj.epoch, "TAI");
            epTab = referenceSystem.daySecondToEpoch([repmat(double(e0.day), nT, 1), double(e0.sec) + tTab(:)], "TAI");
            [Rm, Rd] = referenceSystem.rotationItrsToGcrs(epTab, 'Eop', M.eop);
            t = struct('prepared', true, 't0', t0, 'dt', dt, 'nT', nT);
            t.R = reshape(Rm, 9, nT);
            t.Rdot = reshape(Rd, 9, nT);
            if o.sun.source == "rs.de440"
                if isempty(M.ephemeris)
                    NonGravForces.fail(unit, "depMissing", "sun.source rs.de440 needs the Model's 'Ephemeris' -> pass referenceSystem.readEphemeris(files) at construction");
                end
                [rs, vs] = referenceSystem.geocentricBodyState("sun", epTab, 'Source', "mice", 'Ephemeris', M.ephemeris);
                rm = [];
                if o.srp.lunarEclipse
                    rm = referenceSystem.geocentricBodyState("moon", epTab, 'Source', "mice", 'Ephemeris', M.ephemeris);
                end
            else
                [rs, vs] = referenceSystem.geocentricBodyState("sun", epTab);
                rm = [];
                if o.srp.lunarEclipse
                    rm = referenceSystem.geocentricBodyState("moon", epTab);
                end
            end
            t.sun = rs;
            t.moon = rm;
            t.eclipticNormal = [];
            if o.srp.eclipse == "solaarsCF"
                n = cross(rs, vs, 1);
                t.eclipticNormal = n ./ vecnorm(n);
            end
            t.mass = M.mass.at(epTab);
            t.tPanel = [];
            if o.thermal.model == "transient"
                t.tPanel = obj.thermalTable(e0, tTab);
            end
            % attitude
            t.attitude = obj.attitudeSource;
            t.tq = []; t.q = [];
            if obj.attitudeSource == "provider"
                if isempty(M.attitude)
                    NonGravForces.fail(unit, "inputModel", "attitude 'provider' needs the Model's AttitudeProvider -> construct the Model with 'Attitude', or use 'Attitude', nominal");
                end
                if M.attitude.method ~= "slerp"
                    NonGravForces.fail(unit, "inputModel", "the native evaluator interpolates attitude by slerp; this provider uses %s -> build the provider with 'Method', slerp", M.attitude.method);
                end
                [t.tq, t.q] = obj.attitudeNodes(e0, tTab([1 end]));
            end
            % 3-hour UTC driver bins
            t = obj.driverBins(t, e0, tTab([1 end]), grid);
            % CERES
            t.ceres = struct('nLat', 0, 'nLon', 0, 'nMonth', 0, 'latC', [], 'lonC', [], 'tMonth', [], 'rsut', [], 'rlut', [], 'rsdt', [], 'height', 0);
            if o.erp.backend == "ceres"
                if isempty(M.ceres)
                    NonGravForces.fail(unit, "inputModel", "erp.backend ceres needs the Model's 'Ceres' -> pass io.importCERES(files) at construction");
                end
                c = M.ceres;
                t.ceres = struct('nLat', numel(c.lat_deg), 'nLon', numel(c.lon_deg), 'nMonth', numel(c.mjd), ...
                    'latC', c.lat_deg(:)', 'lonC', c.lon_deg(:)', 'tMonth', c.mjd(:)', ...
                    'rsut', double(c.rsut), 'rlut', double(c.rlut), 'rsdt', double(c.rsdt), 'height', c.referenceHeight_m);
            end
            t = obj.scaleDays(t);
            t.maxColumns = maxCols;
            t.handle = [];
            t.entry = [];
            if obj.native
                t = obj.createNative(t);
            end
        end

        function T = thermalTable(obj, e0, tTab)
            %THERMALTABLE The wall temperatures of thermal.model transient at the table nodes: Model.accel along
            % 'ThermalOrbit', whose transient integration runs once from the steady state at the orbit's first
            % epoch, pchip'd onto the nodes (exact where a node is an orbit epoch). The C reads them by time.
            unit = "DynamicsSource";
            orb = obj.thermalOrbit;
            ep = orb.epoch;
            if string(ep.scale) ~= "TAI", ep = referenceSystem.convertTimeScale(ep, "TAI"); end
            ts = referenceSystem.secondsBetween(e0, ep);
            ts = ts(:)';
            if isempty(ts), ts = NaN; end
            if isnan(ts(1)) || ts(1) > tTab(1) || ts(end) < tTab(end) || any(diff(ts) <= 0) || ~isequal(size(orb.r), [3 numel(ts)]) || ~isequal(size(orb.v), size(orb.r))
                NonGravForces.fail(unit, "dataRange", "'ThermalOrbit' runs from %.1f to %.1f s and the tables need %.1f to %.1f s (span plus 'Margin_s') -> supply the reference orbit (3xK r and v, increasing epochs) over the whole stencil", ...
                    ts(1), ts(end), tTab(1), tTab(end));
            end
            q = [];
            if obj.attitudeSource == "nominal"
                q = NonGravForces.geom.AttitudeProvider.nominal(orb.r, orb.v, obj.alongVelocity);
            end
            [~, info] = obj.model.accel(orb.epoch, orb.r, orb.v, q);
            T = pchip(ts, info.temperature_K, tTab);
        end

        function [tq, q] = attitudeNodes(obj, e0, window)
            %ATTITUDENODES The provider's nodes over the window, sign-continued exactly as
            % referenceSystem.interpolateQuaternions does before it interpolates.
            unit = "DynamicsSource";
            ap = obj.model.attitude;
            nodes = referenceSystem.convertTimeScale(ap.epochNodes, "TAI");
            tn = referenceSystem.secondsBetween(e0, nodes);
            tn = tn(:)';
            i0 = find(tn <= window(1), 1, 'last');
            i1 = find(tn >= window(2), 1, 'first');
            if isempty(i0) || isempty(i1)
                NonGravForces.fail(unit, "dataRange", "the attitude nodes run from %.1f to %.1f s and the run needs %.1f to %.1f s (span plus 'Margin_s') -> supply attitude over the whole stencil, or shrink the margin", tn(1), tn(end), window(1), window(2));
            end
            idx = i0:i1;
            qn = ap.qNodes(:, idx);
            normalized = referenceSystem.normalizeQuaternion(qn);
            reversed = sum(normalized .* qn, 1) < 0;
            normalized(:, reversed) = -normalized(:, reversed);
            pairDot = sum(normalized(:, 2:end) .* normalized(:, 1:end-1), 1);
            stepSign = sign(pairDot);
            stepSign(stepSign == 0) = 1;
            q = normalized .* cumprod([1, stepSign]);
            tq = tn(idx);
        end

        function t = driverBins(obj, t, e0, window, grid)
            %DRIVERBINS The 3-hour UTC bins the space weather resolves: the calendar of each bin and,
            % for the drivers, edges aligned to the run's grid so that the declared breakpoints and the
            % model's switch are the same instant.
            unit = "DynamicsSource";
            M = obj.model;
            needDrivers = M.options.drag.density ~= "exponential" || M.options.drag.wind == "hwm14";
            utc0 = referenceSystem.convertTimeScale(referenceSystem.daySecondToEpoch([double(e0.day) double(e0.sec) + window(1)], "TAI"), "UTC");
            secBin = floor(double(utc0.sec) / 10800) * 10800;
            startUtc = referenceSystem.daySecondToEpoch([double(utc0.day) secBin], "UTC");
            tStart = referenceSystem.secondsBetween(e0, referenceSystem.convertTimeScale(startUtc, "TAI"));
            nb = max(1, ceil((window(2) - tStart) / 10800) + 1);
            t.tb0 = tStart;
            t.nb = nb;
            edges = tStart + (0:nb) * 10800;
            if ~isempty(grid) && isa(grid, 'gjToolbox.Grid')
                % Node-aligned, so a declared breakpoint IS where the drivers change. The grid owns
                % its axis: its indexOf/offset pair rounds to the node, and nothing here rebuilds it.
                for k = 2:nb
                    if grid.contains(edges(k))
                        edges(k) = grid.offset(round(grid.indexOf(edges(k))));
                    end
                end
            end
            t.alignedEdges = edges;
            t.breakpoints = edges(2:end-1);
            t.breakpoints = t.breakpoints(t.breakpoints > window(1) + obj.margin_s & t.breakpoints < window(2) - obj.margin_s);
            [t.f107, t.f107a, t.ap7, t.ap3] = deal(zeros(1, nb), zeros(1, nb), zeros(7, nb), zeros(1, nb));
            [t.doy, t.yy, t.utcSecAtBin, t.utcDay] = deal(zeros(1, nb), zeros(1, nb), zeros(1, nb), zeros(1, nb));
            mid = referenceSystem.daySecondToEpoch([repmat(double(e0.day), nb, 1), double(e0.sec) + (t.tb0 + (0:nb-1)' * 10800 + 1)], "TAI");
            utcMid = referenceSystem.convertTimeScale(mid, "UTC");
            start = referenceSystem.daySecondToEpoch([repmat(double(e0.day), nb, 1), double(e0.sec) + (t.tb0 + (0:nb-1)' * 10800)], "TAI");
            utcStart = referenceSystem.convertTimeScale(start, "UTC");
            dtStart = referenceSystem.epochTo(utcStart, "datetime");
            t.doy = double(day(dtStart(:)', 'dayofyear'));
            t.yy = mod(double(year(dtStart(:)')), 100);
            t.utcSecAtBin = double(utcStart.sec(:))';
            t.utcDay = double(utcStart.day(:))';
            if needDrivers
                if isempty(M.spaceWeather)
                    NonGravForces.fail(unit, "inputModel", "%s needs the Model's 'SpaceWeather' -> pass io.spaceWeather(file) at construction", M.options.drag.density);
                end
                q = M.spaceWeather.query(utcMid);
                t.f107 = q.f107(:)';
                t.f107a = q.f107a(:)';
                t.ap7 = q.apHistory;
                t.ap3 = q.apHistory(2, :);
            end
        end

        function t = scaleDays(obj, t)
            %SCALEDAYS drag.densityScale as one row per UTC day of the driver bins (2026-10-29, O-19): Swarm's
            % bracket where the scale has one (the C interpolates it as env.densityScale does), else ONE factor
            % for the day, taken from the scale's own .at -- the driver scale's factor depends on the day alone,
            % so the table holds exactly what Model.accel evaluates. A day with neither stays "missing" and is
            % refused by the evaluator when a step reaches it, as Model.accel refuses it.
            t.nScaleDay = 0; t.scaleDay0 = 0; t.scale = [];
            M = obj.model;
            if M.options.drag.densityScale == "none"
                return
            end
            S = M.densityScale;
            if isempty(S)
                NonGravForces.fail("DynamicsSource", "inputModel", "drag.densityScale %s needs the Model's 'DensityScale' -> pass env.densityScale or env.densityDriverScale at construction", M.options.drag.densityScale);
            end
            days = min(t.utcDay):max(t.utcDay);
            tab = zeros(6, numel(days));                    % kinds of mex/ngfRhs.h: 0 missing, 1 Swarm's bracket, 2 one factor
            for k = 1:numel(days)
                noon = referenceSystem.daySecondToEpoch([days(k) 43200], "UTC");
                if any(S.kind == ["swarm" "combined"])
                    j = find(S.day == days(k), 1);
                    if ~isempty(j) && ~isnan(S.low(j)) && ~isnan(S.high(j))
                        tab(:, k) = [1; 0; S.low(j); S.high(j); S.hLow_m(j); S.hHigh_m(j)];
                        continue
                    end
                    if S.kind == "swarm", continue, end   % kind 0, missing
                    tab(:, k) = [2; S.fallback.at(noon, 0); 0; 0; 0; 0];
                else
                    tab(:, k) = [2; S.at(noon, 0); 0; 0; 0; 0];
                end
            end
            t.nScaleDay = numel(days); t.scaleDay0 = days(1); t.scale = tab;
        end

        function t = createNative(obj, t)
            %CREATENATIVE Hand the tables to the compiled evaluator and publish its entry.
            unit = "DynamicsSource";
            M = obj.model;
            o = M.options;
            here = fullfile(fileparts(mfilename('fullpath')), 'private');
            if ~isfile(fullfile(here, "rhsGateway." + mexext))
                NonGravForces.fail(unit, "depMissing", "the native evaluator %s is not built -> run NonGravForces.buildMex (target rhs, which needs gjToolbox for csrc/gj_native.h), or construct with 'Native', false", fullfile(here, "rhsGateway." + mexext));
            end
            if o.geom.selfShadow == "raycast"
                NonGravForces.fail(unit, "inputModel", "the native evaluator serves geom.selfShadow none and lut, not raycast (a ray cast per direction is not an integrator's inner loop) -> use lut, or construct with 'Native', false");
            end
            c = NonGravForces.const();
            ell = referenceSystem.ellipsoidParameters("WGS84");
            envDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), '+env', 'private');
            s = struct();
            s.cLight = referenceSystem.iersConstants().speedOfLight; s.sigmaSB = c.sigmaSB; s.kB = c.kB;
            s.au = c.au; s.Rsun = c.Rsun; s.Rmoon = c.Rmoon;
            s.aEarth = ell.a; s.bEarth = ell.b; s.fEarth = ell.f; s.e2Earth = ell.e2; s.tsi = o.srp.tsi; s.amu = c.amu; s.eV = c.eV; s.torr = c.torr;
            s.speciesMass = c.speciesMass(:)' / c.NA;
            s.n = M.macro.n;
            s.area = M.macro.A(:)'; s.normal = M.macro.N; s.centroid = M.macro.C;
            s.optVis = M.macro.optVis; s.optIR = M.macro.optIR;
            s.internalShare = [];
            if o.thermal.internalPower_W > 0
                s.internalShare = M.macro.internalShare(:)';
            end
            s.rcm = M.macro.centreOfMass(:)';
            s.gsi = double(find(["sentman","maxwell","constantCd"] == o.drag.gsi) - 1);
            s.cd = o.drag.cd; s.specular = localOrZero(o.drag.specularFraction); s.twK = o.drag.Tw_K;
            s.cdScale = o.drag.cdScale; s.crScale = o.srp.crScale;
            s.accommodation = double(find(["fixed","pilinski2010","pilinski2011","walker2014","sesam"] == o.drag.accommodation) - 1);
            s.alphaFixed = localOrZero(o.drag.alpha); s.surfaceMassAmu = localOrZero(o.drag.surfaceMass_amu);
            s.density = double(find(["exponential","msis00","msis21"] == o.drag.density) - 1);   % dtm2020 is refused at construction
            s.expRho0 = 0; s.expH0 = 0; s.expH = 1;
            if o.drag.density == "exponential"
                if isempty(M.exponential)
                    NonGravForces.fail(unit, "inputModel", "drag.density exponential needs the Model's 'Exponential' -> pass it at construction");
                end
                s.expRho0 = M.exponential.rho0_kg_m3; s.expH0 = M.exponential.h0_m; s.expH = M.exponential.H_m;
            end
            s.wind = double(find(["none","corotation","hwm14"] == o.drag.wind) - 1);
            s.eclipse = double(find(["cylindrical","conical","solaarsCF"] == o.srp.eclipse) - 1);
            s.lunar = double(o.srp.lunarEclipse);
            % ceresClimatology is refused at CONSTRUCTION (see the constructor), so find() below
            % cannot be handed a name this enum lacks -- which would return empty and pack 0, i.e.
            % "none": a Model that computes Earth radiation silently becoming one that does not.
            s.erp = double(find(["none","pointSource","knocke","ceres"] == o.erp.backend) - 1);
            s.nRings = 0; s.nSectors = 0;
            if any(o.erp.backend == ["knocke", "ceres"])
                s.nRings = o.erp.grid(1); s.nSectors = o.erp.grid(2);
            end
            s.thermal = double(find(["none","steadyState","transient"] == o.thermal.model) - 1);   % transient is the table (mex/ngfRhs.h)
            s.internalPower = o.thermal.internalPower_W;
            s.tPanel = t.tPanel;
            s.shadow = double(o.geom.selfShadow == "lut");
            s.lutResolution = o.geom.lutResolution;
            s.lut = [];
            [s.nChi, s.nArc] = deal(0, 0);
            [s.xChi, s.wChi, s.xArc, s.wArc] = deal([], [], [], []);
            if s.shadow
                s.lut = M.lut.table.f;
                [s.xChi, s.wChi] = localGauss(12);
                [s.xArc, s.wArc] = localGauss(16);
                s.nChi = 12; s.nArc = 16;
            end
            % geom.flowShadow "table": the same table Model.accel reads (geom.flowShadow keeps it for the
            % session), handed to C, which interpolates it in the same order and falls back exactly outside it.
            [s.flowTable, s.ftN, s.ftNs, s.ftNShaded] = deal(0);
            [s.ftHalfWidthTan, s.ftT0, s.ftH, s.ftS0, s.ftHs] = deal(0);
            [s.ftShaded, s.ftPlus, s.ftMinus] = deal([]);
            if s.shadow && o.geom.flowShadow == "table" && o.drag.gsi ~= "constantCd"
                tb = NonGravForces.geom.flowShadow(M.macro, [], [], 'Method', "table", 'Lut', M.lut.table);
                s.flowTable = 1; s.ftN = tb.n; s.ftNs = tb.ns; s.ftNShaded = numel(tb.shaded);
                s.ftHalfWidthTan = tb.halfWidthTan; s.ftT0 = tb.t0; s.ftH = tb.h; s.ftS0 = tb.s0; s.ftHs = tb.hs;
                s.ftShaded = double(tb.shaded); s.ftPlus = tb.plus; s.ftMinus = tb.minus;
            end
            s.t0 = t.t0; s.dt = t.dt; s.nT = t.nT;
            s.R = t.R; s.Rdot = t.Rdot; s.sun = t.sun; s.moon = t.moon; s.eclipticNormal = t.eclipticNormal; s.mass = t.mass;
            s.attitude = double(obj.attitudeSource == "nominal");
            s.alongVelocity = double(obj.alongVelocity);
            s.senseTranspose = double(referenceSystem.attitudeConvention().sense == "vector");
            s.nq = size(t.q, 2); s.tq = t.tq; s.q = t.q;
            s.nb = t.nb; s.tb0 = t.tb0; s.alignedEdges = t.alignedEdges;
            s.f107 = t.f107; s.f107a = t.f107a; s.ap7 = t.ap7; s.ap3 = t.ap3;
            s.doy = t.doy; s.yy = t.yy; s.utcSecAtBin = t.utcSecAtBin; s.utcDay = t.utcDay;
            s.nLat = t.ceres.nLat; s.nLon = t.ceres.nLon; s.nMonth = t.ceres.nMonth;
            s.latC = t.ceres.latC; s.lonC = t.ceres.lonC; s.tMonth = t.ceres.tMonth;
            s.rsut = t.ceres.rsut; s.rlut = t.ceres.rlut; s.rsdt = t.ceres.rsdt; s.ceresHeight = t.ceres.height;
            s.nScaleDay = t.nScaleDay; s.scaleDay0 = t.scaleDay0; s.scale = t.scale;
            s.msis21Dir = char(string(envDir) + filesep);
            s.hwm14Dir = char(envDir);
            s.maxColumns = t.maxColumns;
            [handle, entryPtr] = rhsGateway('create', s);
            t.handle = handle;
            t.entry = struct('AbiVersion', rhsGateway('abi'), 'EntryPtr', entryPtr, 'StateSize', 3, 'MaxColumns', t.maxColumns);
            gjToolbox.DynamicsSource.checkNativeEntry(t.entry, 3);
            if ~isempty(obj.gravity) && ~isempty(obj.gravity.NativeEntry)
                rhsGateway('base', handle, obj.gravity.NativeEntry.EntryPtr);
            elseif ~isempty(obj.gravity)
                NonGravForces.fail(unit, "inputGravity", "'Gravity' publishes no NativeEntry, so the native evaluator cannot add it -> construct with 'Native', false, or give a gravity source with a native entry");
            end
        end
    end
end

function v = localOrZero(x)
%LOCALORZERO An option that is deliberately empty reaches C as a number the model never reads.
if isempty(x)
    v = 0;
else
    v = double(x);
end
end

function [x, w] = localGauss(n)
%LOCALGAUSS Gauss-Legendre nodes and weights, the same Golub-Welsch as geom.flowShadow's.
b = (1:n-1) ./ sqrt(4 * (1:n-1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
x = x';
w = 2 * V(1, i).^2;
end

