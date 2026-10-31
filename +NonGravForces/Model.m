classdef Model
%MODEL Single point of access: non-gravitational forces, torques and partials.
%
%   M = NonGravForces.Model(macro, options, 'Eop', eop, 'SpaceWeather', sw)
%   M = NonGravForces.Model(macro, options, 'Attitude', att, 'Mass', 601.2, 'Eop', eop, 'SpaceWeather', sw)
%
%   [a, info]        = M.accel(epoch, r_gcrs, v_gcrs, q_body2gcrs)   % 3xK m/s^2, body frame
%   [tau, info]      = M.torque(epoch, r, v, q)                      % 3xK N m, body frame
%   [a, tau, info]   = M.evaluate(epoch, r, v, q)                    % both from ONE environment pass
%   T                = M.series(epoch, r, v, q)                      % timetable over chunks
%   [J, names, info] = M.partials(epoch, r, v, q, 'WithRespectTo', ["r","v","cdScale"])
%   B                = M.budget(epoch, r, v, q)                      % per-effect table
%
% METHOD NAME-VALUE ARGUMENTS
%   Every method takes 'Frame', "body" (default) or "gcrs". These are the rest, and none of them
%   appeared in this help until 2026-09-18: a caller could pass them, the code validated them, and
%   the class that is this toolbox's single point of access described them nowhere.
%
%   Method      Name           Size   Type    Units   Description
%   ------      ----           ----   ----    -----   -----------
%   torque      B_body         3xK    double  T       Magnetic field at the satellite, BODY frame.
%   evaluate                                          Required by the magnetic torque and by nothing
%                                                     else; default [] and the term is omitted. Get
%                                                     it from env.magneticField, rotated to the body.
%   torque      Dipole_Am2     3x1    double  A m^2   Residual magnetic dipole of the spacecraft,
%   evaluate                                          body frame; default []. The torque is
%                                                     m x B, so BOTH of these are needed or neither
%                                                     is: one alone leaves the term out silently.
%                                                     torques.magnetic must also be true.
%   accel,      T0_K           Nx1    double  K       Initial panel temperatures for a TRANSIENT
%                                                     thermal model, one per panel; default [] and
%                                                     the steady state is used as the start. Ignored
%                                                     by thermal.model "none" and "steadyState".
%   series      Quantities     1x2    string  -       "accel" (default), "torque", or both. Selects
%                                                     which columns the timetable carries; a torque
%                                                     column costs the torque reduction, which accel
%                                                     alone does not pay for.
%   partials    Method         1x1    string  -       "" (default: the Options' partials.method),
%                                                     "analytic" or "finiteDiff". analytic uses
%                                                     kernel.partials' closed forms; finiteDiff is a
%                                                     central difference in the parameter, for
%                                                     cross-checking the analytic form. A complex
%                                                     step is NOT offered and is refused by name:
%                                                     MATLAB's erf, erfc and every norm in the chain
%                                                     discard the imaginary part (measured
%                                                     2026-09-08).
%   partials    WithRespectTo  1xP    string  -       Any of "r", "v" (default ["r" "v"]),
%                                                     "cdScale", "crScale", "alpha". alpha is
%                                                     refused for a two-state surface, which
%                                                     computes its own coverage.
%
%   The Model owns the packed macro model, the attitude provider, the mass, the
%   Earth-orientation, space-weather, ephemeris and CERES sources, the resolved
%   Options and the self-shadow table. It computes forces and torques and nothing
%   else: time, frames and attitude kinematics are referenceSystem, reading is
%   gracefoReading, integration is gjToolbox (through app.DynamicsSource).
%
%   ONE KERNEL (principle P2). accel and torque are two reductions of the same
%   per-panel force array F_k: F = sum_k F_k and tau = sum_k (c_k - r_cm) x F_k,
%   computed once per epoch. budget is the same array kept per effect instead of
%   summed, and it adds back to accel bitwise. evaluate returns both reductions
%   from one pass; series uses it when both are asked for.
%
%   WHAT IT COMPOSES, in order, per epoch:
%     1. the relative velocity            env.wind (co-rotation from referenceSystem)
%     2. the geodetic height              referenceSystem.cartesianToGeodetic on the ITRS position
%     3. density and composition          env.density with the space-weather drivers, times the observed
%                                         daily scale when drag.densityScale is "swarm" (env.densityScale,
%                                         observed) or "driver" (env.densityDriverScale, predicted)
%     4. the radiation source list        env.radiationSources (Sun, Earth tiles, eclipse)
%     5. occlusion                        geom.selfShadow per radiation source; geom.flowShadow for the gas
%                                         (Method "table" when geom.flowShadow is "table", the standard preset)
%     6. panel temperatures               geom.thermalState with the per-source occlusion
%     7. accommodation                    kernel.accommodation (fixed, pilinski2010/2011, walker2014)
%     8. per-panel forces                 kernel.dragPanelForce, radiationPanelForce, thermalRecoil
%     9. reductions                       kernel.assemble, and kernel.bodyTorques for the torque
%   Every step is the unit that owns that fact; the Model adds no physics. Until
%   2026-09-15 steps 5 and 6 fed a radiation-beam fraction to the gas (deleting the
%   edge-on shear) and the Sun's fraction to every source's heating (deleting the
%   Earth heating of panels facing away from the Sun): HANDOVER LP-C, F1 and F2.
%
%   FRAME. Results are BODY frame by default (principle P4: the body frame is
%   what an accelerometer reads and what the sensor layer maps from). Pass
%   'Frame', "gcrs" for the inertial frame, which is what an integrator wants.
%
% CONSTRUCTOR INPUTS
%   Name     Size  Type                            Units  Description
%   ----     ----  ----                            -----  -----------
%   macro    1x1   NonGravForces.geom.MacroModel   -      Packed panel geometry
%   options  1x1   NonGravForces.Options           -      Resolved options (default "standard")
%
% CONSTRUCTOR NAME-VALUE ARGUMENTS
%   Name          Size  Type            Units  Description
%   ----          ----  ----            -----  -----------
%   Attitude      1x1   AttitudeProvider -     BODY->GCRS quaternions, used when a method is
%                                              called without q. Default [].
%   Mass          1x1   double|struct   kg     A constant mass, or an io.massModel result.
%                                              Default: the macro model's mass_kg (geom.graceFO
%                                              carries the handbook value); a macro model without
%                                              one is REFUSED, never given GRACE-C's (LP-C F10).
%   Eop           1x1   struct          -      referenceSystem.readEop result or a file path.
%                                              REQUIRED to evaluate: the geodetic height needs
%                                              the ITRS rotation and the co-rotation its rate.
%   SpaceWeather  1x1   struct          -      io.spaceWeather result. Required by msis00, msis21 and dtm2020.
%   Ephemeris     1x1   struct          -      referenceSystem.readEphemeris result, required
%                                              when sun.source is "rs.de440".
%   Exponential   1x1   struct          -      .rho0_kg_m3 .h0_m .H_m, required when
%                                              drag.density is "exponential".
%   Ceres         1x1   struct          -      io.importCERES result, required when erp.backend is "ceres".
%   Climatology   1x1   struct          -      .albedo and .longwave io.importGroopsGrid results,
%                                              required when erp.backend is "ceresClimatology".
%   DensityScale  1x1   struct          -      env.densityScale result (drag.densityScale "swarm"; with
%                                              'Fallback' "combined") or env.densityDriverScale result
%                                              ("driver"), required by those
%                                              and refused otherwise; its kind must match the option and it
%                                              must have been formed for this Model's drag.density. Density and every species
%                                              density are multiplied by S.at(epoch, geodetic height), so
%                                              the gas force scales exactly (it is linear in them).
%
% OUTPUT (accel / torque / evaluate)
%   a, tau   3xK   double         acceleration [m/s^2] or torque [N m]
%   info     1x1   struct         .frame, .drag .radiation .thermal (3xK each, the
%                                 same frame), .bodyTorque (torque only), .rho (1xK),
%                                 .nu (1xK shadow factor), .mass_kg (1xK), .vrel_body,
%                                 .sources (M), .temperature_K (NxK), .flowShadow,
%                                 .accommodation, .lla, .densityScale (1xK, ones with
%                                 drag.densityScale "none"), .backends, .provenance
%
% ACCURACY
%   The Model adds no approximation: every number comes from the unit that owns
%   it, and each of those states its own error. What the Model does own is the
%   ORDER, which is fixed above, and the resolved options, which travel in
%   info.provenance (rule R2). The dominant error is the thermospheric density,
%   10-30 % (HANDOVER Part C.4), not anything here. Cost per epoch is recorded by
%   tools/benchmark.m.
%
% ERRORS
%   NonGravForces:Model:inputMacro      macro is not a MacroModel
%   NonGravForces:Model:inputOptions    options is not an Options object, or its choices contradict each other
%                                       (alpha beside a non-fixed accommodation, maxwell without a specular
%                                       fraction, walker2014 without a surface mass, internal power without shares,
%                                       geom.flowShadow table beside selfShadow raycast, a 'DensityScale' with
%                                       drag.densityScale "none", of the other kind, or formed for another drag.density)
%   NonGravForces:Model:inputMass       no 'Mass' and no mass on the macro model, or 'Mass' not a positive mass / massModel
%   NonGravForces:Model:inputState      epoch, r, v or q disagree in size, or are not finite
%   NonGravForces:Model:inputAttitude   no q was given and the Model has no attitude provider
%   NonGravForces:Model:dependEop       evaluated without 'Eop'
%   NonGravForces:Model:dependWeather   evaluated without 'SpaceWeather' or 'Exponential' for the density backend, or without 'SpaceWeather' for drag.wind hwm14
%   NonGravForces:Model:dependDensityScale  evaluated with drag.densityScale "swarm", "driver" or "combined" and no 'DensityScale'
%   NonGravForces:Model:dependCeres     evaluated with erp.backend ceres and no 'Ceres', or
%                                       ceresClimatology and no 'Climatology'
%   NonGravForces:Model:unsupported     a method option this Model cannot honour (partials of alpha for a
%                                       computed accommodation, transient thermal under parfor, unknown quantities)
%   NonGravForces:<unit>:*              propagated unwrapped from the unit that refused
%
% REFERENCES
%   HANDOVER Part D.3 (workflow), Part G LP-6 and LP-C; theory sections 2, 3, 5 and 7;
%   Siemes C. et al. (2023) Adv. Space Res.
%
% EXAMPLE
%   eop = referenceSystem.readEop(eopFile);
%   sw  = NonGravForces.io.spaceWeather(swFile);
%   [ep, r, v] = NonGravForces.io.readGNV1B(gnv, 'Eop', eop);
%   [~, q]     = NonGravForces.io.readSCA1B(sca);
%   S = NonGravForces.env.densityScale(swarmFiles, sw, 'Fallback', NonGravForces.env.densityDriverScale(sw, cls));   % the standard preset's density scale
%   M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("standard"), ...
%                           'Eop', eop, 'SpaceWeather', sw, 'DensityScale', S);
%   [a, info] = M.accel(ep, r, v, q);
%   vecnorm(info.drag)          % ~1e-8 m/s^2 at solar minimum
%
% See also NonGravForces.Options, NonGravForces.geom.MacroModel, NonGravForces.app.DynamicsSource, NonGravForces.kernel.assemble.
%
% KEYWORDS model accel torque series partials budget single point of access
% Created 2026-08-27 by Claude Fable 5 (Anthropic); flow shadowing, per-source heating, honest options 2026-09-15 by Claude Opus 5 (Anthropic).
    properties (SetAccess = private)
        macro
        options (1,1) NonGravForces.Options
        attitude
        mass struct
        eop
        spaceWeather
        ephemeris
        exponential
        ceres
        climatology
        densityScale
        lut struct
        provenance struct
    end
    methods
        function obj = Model(macro, options, nv)
            arguments
                macro
                options = NonGravForces.Options("standard")
                nv.Attitude = []
                nv.Mass = []
                nv.Eop = []
                nv.SpaceWeather = []
                nv.Ephemeris = []
                nv.Exponential = []
                nv.Ceres = []
                nv.Climatology = []
                nv.DensityScale = []
            end
            if ~isa(options, 'NonGravForces.Options') || ~isscalar(options)
                NonGravForces.fail("Model", "inputOptions", ...
                    "options must be a scalar NonGravForces.Options -> build one with NonGravForces.Options(preset)");
            end
            if ~isa(macro, 'NonGravForces.geom.MacroModel')
                NonGravForces.fail("Model", "inputMacro", ...
                    "macro must be a NonGravForces.geom.MacroModel -> build one with geom.graceFO or geom.boxWing");
            end
            localCheckOptions(options, macro);
            obj.macro = macro;
            obj.options = options;
            obj.attitude = nv.Attitude;
            obj.mass = localMass(nv.Mass, macro);
            obj.eop = nv.Eop;
            obj.spaceWeather = nv.SpaceWeather;
            obj.ephemeris = nv.Ephemeris;
            obj.exponential = nv.Exponential;
            obj.ceres = nv.Ceres;
            obj.climatology = nv.Climatology;
            localCheckDensityScale(options, nv.DensityScale);
            obj.densityScale = nv.DensityScale;
            % The self-shadow table costs a ray cast per cube-map node, so it is built once here.
            if options.geom.selfShadow == "lut"
                obj.lut = struct('built', true, 'table', ...
                    NonGravForces.geom.selfShadow(macro, [], 'Method', "lut", 'LutResolution', options.geom.lutResolution));
            else
                obj.lut = struct('built', false, 'table', []);
            end
            obj.provenance = NonGravForces.provenance("Model", 'Options', options.toStruct());
        end

        function [a, info] = accel(obj, epoch, r, v, q, nv)
            %ACCEL Non-gravitational acceleration, body frame by default.
            arguments
                obj
                epoch (1,1) struct
                r double
                v double
                q double = []
                nv.Frame (1,1) string {mustBeMember(nv.Frame, ["body","gcrs"])} = "body"
                nv.T0_K double = []
            end
            st = obj.environment(epoch, r, v, q, nv.T0_K);
            [a, info] = obj.reduceAccel(st, nv.Frame);
            if nargout > 1, info.provenance = localStamp(obj.provenance, "Model.accel"); end
        end

        function [tau, info] = torque(obj, epoch, r, v, q, nv)
            %TORQUE Non-gravitational surface torque about the centre of mass, plus the enabled body torques.
            arguments
                obj
                epoch (1,1) struct
                r double
                v double
                q double = []
                nv.Frame (1,1) string {mustBeMember(nv.Frame, ["body","gcrs"])} = "body"
                nv.B_body double = []
                nv.Dipole_Am2 double = []
            end
            st = obj.environment(epoch, r, v, q, []);
            [tau, info] = obj.reduceTorque(st, nv.Frame, nv.B_body, nv.Dipole_Am2);
            if nargout > 1, info.provenance = localStamp(obj.provenance, "Model.torque"); end
        end

        function [a, tau, info] = evaluate(obj, epoch, r, v, q, nv)
            %EVALUATE Acceleration AND torque from one environment pass (principle P2 at the method level).
            arguments
                obj
                epoch (1,1) struct
                r double
                v double
                q double = []
                nv.Frame (1,1) string {mustBeMember(nv.Frame, ["body","gcrs"])} = "body"
                nv.B_body double = []
                nv.Dipole_Am2 double = []
                nv.T0_K double = []
            end
            st = obj.environment(epoch, r, v, q, nv.T0_K);
            [a, info] = obj.reduceAccel(st, nv.Frame);
            [tau, ti] = obj.reduceTorque(st, nv.Frame, nv.B_body, nv.Dipole_Am2);
            info.bodyTorque = ti.bodyTorque;
            info.torque = struct('drag', ti.drag, 'radiation', ti.radiation, 'thermal', ti.thermal);
            info.provenance = localStamp(obj.provenance, "Model.evaluate");
        end

        function B = budget(obj, epoch, r, v, q, nv)
            %BUDGET Per-effect, per-axis accelerations as a table; the effects add back to accel bitwise.
            arguments
                obj
                epoch (1,1) struct
                r double
                v double
                q double = []
                nv.Frame (1,1) string {mustBeMember(nv.Frame, ["body","gcrs"])} = "body"
            end
            [a, info] = obj.accel(epoch, r, v, q, 'Frame', nv.Frame);
            eff = ["drag"; "radiation"; "thermal"; "total"];
            vals = {info.drag, info.radiation, info.thermal, a};
            rms = zeros(numel(eff), 3); mx = zeros(numel(eff), 3);
            for i = 1:numel(eff)
                rms(i, :) = sqrt(mean(vals{i}.^2, 2))';
                mx(i, :) = max(abs(vals{i}), [], 2)';
            end
            B = table(eff, rms(:, 1), rms(:, 2), rms(:, 3), mx(:, 1), mx(:, 2), mx(:, 3), vecnorm(rms, 2, 2), ...
                'VariableNames', {'effect', 'rms_x', 'rms_y', 'rms_z', 'max_x', 'max_y', 'max_z', 'rms_norm'});
            B.Properties.VariableUnits = ["", repmat("m/s^2", 1, 7)];
            B.Properties.Description = sprintf('NonGravForces.Model.budget, %s frame, %d epochs, preset %s', ...
                nv.Frame, info.K, obj.options.preset);
        end

        function T = series(obj, epoch, r, v, q, nv)
            %SERIES accel (and optionally torque) over many epochs, in chunks, optionally under parfor.
            arguments
                obj
                epoch (1,1) struct
                r double
                v double
                q double = []
                nv.Quantities (1,:) string = "accel"
                nv.Frame (1,1) string {mustBeMember(nv.Frame, ["body","gcrs"])} = "body"
            end
            if ~all(ismember(nv.Quantities, ["accel", "torque"]))
                NonGravForces.fail("Model", "unsupported", "'Quantities' must be accel and/or torque -> got %s", strjoin(nv.Quantities, ", "));
            end
            K = numel(epoch.sec);
            obj.checkState(epoch, r, v, q, K);
            chunk = max(1, round(obj.options.series.chunk));
            first = 1:chunk:K;
            last = min(first + chunk - 1, K);
            nC = numel(first);
            aC = cell(1, nC); tC = cell(1, nC);
            wantA = any(nv.Quantities == "accel"); wantT = any(nv.Quantities == "torque");
            transient = obj.options.thermal.model == "transient";
            if obj.options.series.useParallel && transient
                NonGravForces.fail("Model", "unsupported", "thermal.model transient integrates the panel temperatures from one chunk into the next -> set series.useParallel false, or use steadyState");
            end
            if obj.options.series.useParallel
                parfor c = 1:nC
                    [aC{c}, tC{c}] = obj.chunkEval(epoch, r, v, q, first(c), last(c), wantA, wantT, nv.Frame, []); %#ok<PFBNS>
                end
            else
                T0 = [];
                for c = 1:nC
                    [aC{c}, tC{c}, T0] = obj.chunkEval(epoch, r, v, q, first(c), last(c), wantA, wantT, nv.Frame, T0);
                end
            end
            t = referenceSystem.epochTo(epoch, "datetime");
            T = timetable(t(:));
            if wantA
                a = [aC{:}];
                T.ax = a(1, :)'; T.ay = a(2, :)'; T.az = a(3, :)';
            end
            if wantT
                tq = [tC{:}];
                T.tx = tq(1, :)'; T.ty = tq(2, :)'; T.tz = tq(3, :)';
            end
            T.Properties.Description = sprintf('NonGravForces.Model.series, %s frame, %d epochs in %d chunk(s), parallel %d, preset %s', ...
                nv.Frame, K, nC, obj.options.series.useParallel, obj.options.preset);
        end

        function [J, names, info] = partials(obj, epoch, r, v, q, nv)
            %PARTIALS Derivatives of the acceleration with respect to the state and the estimated parameters.
            arguments
                obj
                epoch (1,1) struct
                r double
                v double
                q double = []
                nv.WithRespectTo (1,:) string = ["r", "v"]
                nv.Method (1,1) string = ""
                nv.Frame (1,1) string {mustBeMember(nv.Frame, ["body","gcrs"])} = "body"
            end
            known = ["r", "v", "cdScale", "crScale", "alpha"];
            axisName = ["x", "y", "z"];
            bad = setdiff(nv.WithRespectTo, known);
            if ~isempty(bad)
                NonGravForces.fail("Model", "unsupported", "'WithRespectTo' entries %s are not model parameters -> use %s", ...
                    strjoin(bad, ", "), strjoin(known, ", "));
            end
            method = nv.Method;
            if method == "", method = obj.options.partials.method; end
            if ~any(method == ["analytic", "finiteDiff"])
                NonGravForces.fail("Model", "unsupported", ...
                    "'Method' %s is not offered -> use analytic or finiteDiff (a complex step cannot pass MATLAB's erf/erfc or any norm, measured 2026-09-08)", method);
            end
            if any(nv.WithRespectTo == "alpha") && obj.options.drag.accommodation ~= "fixed"
                NonGravForces.fail("Model", "unsupported", ...
                    "alpha is a parameter only of accommodation fixed; %s computes it from the environment -> differentiate cdScale instead", obj.options.drag.accommodation);
            end
            K = numel(epoch.sec);
            obj.checkState(epoch, r, v, q, K);
            names = string.empty(1, 0);
            cols = {};
            scalars = intersect(nv.WithRespectTo, ["cdScale", "crScale", "alpha"], 'stable');
            if ~isempty(scalars)
                if method == "analytic"
                    st = obj.environment(epoch, r, v, q, []);          % one pass for every scalar parameter
                end
            end
            for p = nv.WithRespectTo
                switch p
                    case {"r", "v"}
                        base = r; if p == "v", base = v; end
                        scale = median(vecnorm(base));
                        for i = 1:3
                            h = eps^(1/3) * max(scale, 1);      % central differences: truncation ~ h^2, rounding ~ eps/h
                            dp = zeros(3, K); dp(i, :) = h;
                            if p == "r"
                                aP = obj.accel(epoch, r + dp, v, q, 'Frame', nv.Frame);
                                aM = obj.accel(epoch, r - dp, v, q, 'Frame', nv.Frame);
                            else
                                aP = obj.accel(epoch, r, v + dp, q, 'Frame', nv.Frame);
                                aM = obj.accel(epoch, r, v - dp, q, 'Frame', nv.Frame);
                            end
                            cols{end+1} = (aP - aM) / (2 * h); %#ok<AGROW>
                            names(end+1) = p + axisName(i); %#ok<AGROW>
                        end
                    otherwise
                        if method == "finiteDiff"
                            cols{end+1} = obj.parameterPartialNumeric(p, epoch, r, v, q, nv.Frame); %#ok<AGROW>
                        else
                            cols{end+1} = obj.parameterPartial(p, st, nv.Frame); %#ok<AGROW>
                        end
                        names(end+1) = p; %#ok<AGROW>
                end
            end
            J = zeros(3, numel(cols), K);
            for i = 1:numel(cols)
                J(:, i, :) = reshape(cols{i}, 3, 1, K);
            end
            info.method = method;
            info.frame = nv.Frame;
            info.stateStep = struct('kind', "central difference", 'note', "h = eps^(1/3) * median(|r| or |v|); truncation h^2, rounding eps/h");
            info.analytic = intersect(nv.WithRespectTo, ["cdScale", "crScale", "alpha"]);
            info.K = K;
            info.provenance = NonGravForces.provenance("Model.partials", 'Options', struct('method', method, 'wrt', strjoin(nv.WithRespectTo, ",")));
        end
    end

    methods (Access = private)
        function st = environment(obj, epoch, r, v, q, T0)
            %ENVIRONMENT The one pass every method reduces from: environment, occlusion, per-panel forces.
            K = numel(epoch.sec);
            q = obj.resolveAttitude(epoch, q, K);
            obj.checkState(epoch, r, v, q, K);
            o = obj.options;
            if isempty(obj.eop)
                NonGravForces.fail("Model", "dependEop", "evaluating needs 'Eop' at construction -> the geodetic height and the co-rotation both come from the ITRS rotation, and neither may be guessed");
            end
            % 1. relative velocity, GCRS then body
            windArgs = {};
            if o.drag.wind == "hwm14"
                if isempty(obj.spaceWeather)
                    NonGravForces.fail("Model", "dependWeather", "drag.wind hwm14 needs 'SpaceWeather' at construction for the 3-hour ap -> pass io.spaceWeather(file)");
                end
                windArgs = {'SpaceWeather', obj.spaceWeather};
            end
            [vrel, windInfo] = NonGravForces.env.wind(o.drag.wind, epoch, r, v, 'Eop', obj.eop, windArgs{:});
            R = referenceSystem.rotationFromQuaternion(q);                 % body -> GCRS
            vrelBody = reshape(pagemtimes(pagetranspose(R), reshape(vrel, 3, 1, K)), 3, K);
            % 2. geodetic height from the ITRS position
            [Ritrs, ~] = referenceSystem.rotationItrsToGcrs(epoch, 'Eop', obj.eop);
            rItrs = reshape(pagemtimes(pagetranspose(Ritrs), reshape(r, 3, 1, K)), 3, K);
            lla = referenceSystem.cartesianToGeodetic(rItrs', 'Ellipsoid', "WGS84")';   % the dependency is Nx3
            % 3. density and composition
            [rho, species, Tgas, fScale] = obj.callDensity(epoch, lla);
            % 4. radiation sources
            srcArgs = {};
            if ~isempty(obj.ephemeris), srcArgs = [srcArgs, {'Ephemeris', obj.ephemeris}]; end
            if o.erp.backend == "ceres"
                if isempty(obj.ceres)
                    NonGravForces.fail("Model", "dependCeres", "erp.backend ceres needs 'Ceres' at construction -> pass io.importCERES(files); the Earth's fluxes are never assumed");
                end
                srcArgs = [srcArgs, {'Ceres', obj.ceres, 'RotationItrsToGcrs', Ritrs}];
            elseif o.erp.backend == "ceresClimatology"
                if isempty(obj.climatology)
                    NonGravForces.fail("Model", "dependCeres", "erp.backend ceresClimatology needs 'Climatology' at construction -> pass struct('albedo', io.importGroopsGrid(albedoFile), 'longwave', io.importGroopsGrid(longwaveFile)); the Earth's fluxes are never assumed");
                end
                srcArgs = [srcArgs, {'Climatology', obj.climatology, 'RotationItrsToGcrs', Ritrs}];
            end
            [src, srcInfo] = NonGravForces.env.radiationSources(epoch, r, q, o, srcArgs{:});
            % 5. occlusion: per radiation source, and flux-weighted for the gas
            kinetic = o.drag.gsi ~= "constantCd";
            [flow, fRad] = obj.occlusion(vrelBody, Tgas, src, K, kinetic);
            % 6. panel temperatures, heated source by source through that source's own occlusion
            switch o.thermal.model
                case "none"
                    Tpanel = zeros(obj.macro.n, K);
                otherwise
                    thArgs = {'Exposed', fRad};
                    if o.thermal.internalPower_W > 0
                        thArgs = [thArgs, {'InternalPower_W', o.thermal.internalPower_W, 'InternalShare', obj.macro.internalShare}];
                    end
                    if o.thermal.model == "transient"
                        if isempty(T0)
                            T0 = NonGravForces.geom.thermalState(obj.macro, localFirstEpoch(src), thArgs{1}, fRad(:, 1, :), thArgs{3:end});
                        end
                        thArgs = [thArgs, {'Model', "transient", 'Epoch', epoch, 'T0_K', T0(:)'}];
                    end
                    Tpanel = NonGravForces.geom.thermalState(obj.macro, src, thArgs{:});
            end
            % 7-8. accommodation and per-panel forces
            [Fdrag, acc] = obj.dragForces(rho, species, Tgas, vrelBody, flow);
            st.Fdrag = o.drag.cdScale * Fdrag;
            st.Frad = o.srp.crScale * NonGravForces.kernel.radiationPanelForce(obj.macro, src, 'Shadow', fRad);
            if o.thermal.model == "none"
                st.Fthermal = zeros(3, obj.macro.n, K);
            else
                st.Fthermal = NonGravForces.kernel.thermalRecoil(obj.macro, Tpanel);
            end
            st.K = K; st.R = R; st.r = r; st.rho = rho; st.species = species; st.Tgas = Tgas;
            st.vrelBody = vrelBody; st.lla = lla; st.src = src; st.nu = src.nu; st.densityScale = fScale;
            st.flow = flow; st.shadowRad = fRad; st.Tpanel = Tpanel; st.acc = acc;
            st.mass_kg = obj.mass.at(epoch);
            st.backends = struct('density', o.drag.density, 'wind', windInfo.model, 'gsi', o.drag.gsi, ...
                'accommodation', o.drag.accommodation, 'eclipse', o.srp.eclipse, 'erp', o.erp.backend, ...
                'sun', srcInfo.sunSource, 'thermal', o.thermal.model, 'selfShadow', o.geom.selfShadow);
        end

        function [a, info] = reduceAccel(obj, st, frame)
            m = reshape(st.mass_kg, 1, []);
            aD = NonGravForces.kernel.assemble(st.Fdrag, obj.macro) ./ m;
            aR = NonGravForces.kernel.assemble(st.Frad, obj.macro) ./ m;
            aT = NonGravForces.kernel.assemble(st.Fthermal, obj.macro) ./ m;
            a = aD + aR + aT;
            if frame == "gcrs"
                [a, aD, aR, aT] = obj.toGcrs(st.R, a, aD, aR, aT);
            end
            info = obj.stateInfo(st, frame);
            info.drag = aD; info.radiation = aR; info.thermal = aT;
        end

        function [tau, info] = reduceTorque(obj, st, frame, B_body, dipole)
            [~, tD] = NonGravForces.kernel.assemble(st.Fdrag, obj.macro);
            [~, tR] = NonGravForces.kernel.assemble(st.Frad, obj.macro);
            [~, tT] = NonGravForces.kernel.assemble(st.Fthermal, obj.macro);
            o = obj.options.torques;
            if o.gravityGradient || o.magnetic
                rBody = reshape(pagemtimes(pagetranspose(st.R), reshape(st.r, 3, 1, [])), 3, []);
                tB = NonGravForces.kernel.bodyTorques(obj.macro, rBody, ...
                    'GravityGradient', o.gravityGradient, 'Magnetic', o.magnetic, ...
                    'B_body', B_body, 'Dipole_Am2', dipole);
            else
                tB = zeros(3, st.K);
            end
            tau = tD + tR + tT + tB;
            if frame == "gcrs"
                [tau, tD, tR, tT, tB] = obj.toGcrs(st.R, tau, tD, tR, tT, tB);
            end
            info = obj.stateInfo(st, frame);
            info.drag = tD; info.radiation = tR; info.thermal = tT; info.bodyTorque = tB;
        end

        function [rho, species, Tgas, f] = callDensity(obj, epoch, lla)
            %CALLDENSITY The density backend with the drivers the Model was given, never a default set,
            %times the daily scale when drag.densityScale is "swarm" (observed), "driver" (predicted) or "combined".
            o = obj.options;
            if o.drag.density == "exponential"
                if isempty(obj.exponential)
                    NonGravForces.fail("Model", "dependWeather", "drag.density exponential needs 'Exponential' at construction (.rho0_kg_m3 .h0_m .H_m) -> no scale height is assumed");
                end
                [rho, species, T] = NonGravForces.env.density("exponential", epoch, lla, [], 'Exponential', obj.exponential);
            else
                if isempty(obj.spaceWeather)
                    NonGravForces.fail("Model", "dependWeather", "drag.density %s needs 'SpaceWeather' at construction -> pass io.spaceWeather(file); the drivers are never assumed", o.drag.density);
                end
                [rho, species, T] = NonGravForces.env.density(o.drag.density, epoch, lla, obj.spaceWeather);
            end
            Tgas = T(2, :);
            f = ones(1, numel(rho));
            if o.drag.densityScale ~= "none"
                if isempty(obj.densityScale)
                    NonGravForces.fail("Model", "dependDensityScale", "drag.densityScale %s needs 'DensityScale' at construction -> pass env.densityScale(swarmFiles, sw) for swarm, env.densityDriverScale(sw, cls) for driver, env.densityScale(swarmFiles, sw, 'Fallback', D) for combined", o.drag.densityScale);
                end
                f = reshape(obj.densityScale.at(epoch, lla(3, :)), 1, []);
                rho = rho .* f;
                species = species .* f;               % every species alike: composition and T are the model's
            end
        end

        function [F, acc] = dragForces(obj, rho, species, Tgas, vrelBody, flow)
            %DRAGFORCES The gas-surface model and the accommodation the options name.
            o = obj.options;
            shadowArgs = {};
            if ~isempty(flow), shadowArgs = {'Shadow', flow}; end
            if o.drag.gsi == "constantCd"
                acc = [];
                F = NonGravForces.kernel.dragPanelForce(obj.macro, rho, [], [], vrelBody, ...
                    'Model', "constantCd", 'Cd', o.drag.cd, shadowArgs{:});
                return
            end
            switch o.drag.accommodation
                case "fixed"
                    acc = NonGravForces.kernel.accommodation("fixed", species, Tgas, 'Alpha', o.drag.alpha);
                case {"pilinski2010", "pilinski2011"}
                    acc = NonGravForces.kernel.accommodation(o.drag.accommodation, species, Tgas);
                case "walker2014"
                    acc = NonGravForces.kernel.accommodation("walker2014", species, Tgas, 'SurfaceMass_amu', o.drag.surfaceMass_amu);
                case "sesam"
                    % SESAM's pressure is a RAM pressure, so it needs the flow. vrelBody is already
                    % the gas-relative velocity (env.wind), which is exactly the V_r of its eq. (15).
                    acc = NonGravForces.kernel.accommodation("sesam", species, Tgas, 'RelativeSpeed_ms', vecnorm(vrelBody, 2, 1));
            end
            extra = {};
            if o.drag.gsi == "maxwell", extra = {'SpecularFraction', o.drag.specularFraction}; end
            F = NonGravForces.kernel.dragPanelForce(obj.macro, rho, species, Tgas, vrelBody, ...
                'Model', o.drag.gsi, 'Alpha', acc, 'Tw_K', o.drag.Tw_K, extra{:}, shadowArgs{:});
        end

        function [flow, fRad] = occlusion(obj, vrelBody, Tgas, src, K, kinetic)
            %OCCLUSION geom.flowShadow for the gas stream, geom.selfShadow for every radiation source.
            N = obj.macro.n; M = size(src.dir, 2);
            method = obj.options.geom.selfShadow;
            if method == "none"
                flow = []; fRad = ones(N, K, M);
                return
            end
            args = {'Method', method};
            if method == "lut", args = [args, {'Lut', obj.lut.table}]; end   % built once at construction
            T = [];
            if kinetic, T = Tgas; end                                        % constantCd reads the parallel beam only
            flowArgs = args;
            if obj.options.geom.flowShadow == "table" && kinetic
                flowArgs{2} = "table";                                       % the lut ratios, tabulated (geom.flowShadow)
            end
            flow = NonGravForces.geom.flowShadow(obj.macro, vrelBody, T, flowArgs{:});
            dirAll = reshape(src.dir, 3, M * K);                            % every source direction of every epoch
            fRad = permute(reshape(NonGravForces.geom.selfShadow(obj.macro, dirAll, args{:}), N, M, K), [1 3 2]);
        end

        function q = resolveAttitude(obj, epoch, q, K)
            if ~isempty(q), return; end
            if isempty(obj.attitude)
                NonGravForces.fail("Model", "inputAttitude", "no q was given and this Model has no 'Attitude' provider -> pass q, or construct with 'Attitude', geom.AttitudeProvider(...)");
            end
            q = obj.attitude.at(epoch);
            if size(q, 2) ~= K
                NonGravForces.fail("Model", "inputState", "the attitude provider returned %d quaternions for %d epochs", size(q, 2), K);
            end
        end

        function checkState(obj, epoch, r, v, q, K) %#ok<INUSD>
            if nargin < 6, K = numel(epoch.sec); end
            if ~(ismatrix(r) && isequal(size(r), [3 K]) && isequal(size(v), [3 K]))
                NonGravForces.fail("Model", "inputState", "r and v must be 3xK with K = %d epochs -> got %s and %s", K, mat2str(size(r)), mat2str(size(v)));
            end
            if ~isempty(q) && ~isequal(size(q), [4 K])
                NonGravForces.fail("Model", "inputState", "q must be 4xK with K = %d epochs -> got %s", K, mat2str(size(q)));
            end
            if any(~isfinite(r(:))) || any(~isfinite(v(:))) || (~isempty(q) && any(~isfinite(q(:))))
                NonGravForces.fail("Model", "inputState", "the state contains NaN/Inf -> screen it with qc before evaluating");
            end
        end

        function [a, b, c, d, e] = toGcrs(~, R, a, b, c, d, e)
            %TOGCRS Rotate up to five body-frame 3xK arrays into GCRS with the same attitude.
            rot = @(x) reshape(pagemtimes(R, reshape(x, 3, 1, [])), 3, []);
            a = rot(a);
            if nargin > 3, b = rot(b); end
            if nargin > 4, c = rot(c); end
            if nargin > 5, d = rot(d); end
            if nargin > 6, e = rot(e); end
        end

        function info = stateInfo(~, st, frame)
            info = struct('frame', frame, 'K', st.K, 'rho', st.rho, 'nu', st.nu, ...
                'mass_kg', st.mass_kg, 'vrel_body', st.vrelBody, 'sources', size(st.src.dir, 2), ...
                'temperature_K', st.Tpanel, 'flowShadow', st.flow, 'accommodation', st.acc, 'lla', st.lla, ...
                'densityScale', st.densityScale, 'backends', st.backends);
        end

        function [a, tq, Tlast] = chunkEval(obj, epoch, r, v, q, i0, i1, wantA, wantT, frame, T0)
            %CHUNKEVAL One chunk of series from ONE environment pass; the transient temperatures carry over.
            idx = i0:i1;
            ep = referenceSystem.daySecondToEpoch([epoch.day(idx) epoch.sec(idx)], epoch.scale);
            qc = []; if ~isempty(q), qc = q(:, idx); end
            st = obj.environment(ep, r(:, idx), v(:, idx), qc, T0);
            a = []; tq = [];
            if wantA, a = obj.reduceAccel(st, frame); end
            if wantT, tq = obj.reduceTorque(st, frame, [], []); end
            Tlast = st.Tpanel(:, end);
        end

        function d = parameterPartial(obj, p, st, frame)
            %PARAMETERPARTIAL Analytic derivatives of the acceleration with respect to a scalar parameter.
            m = reshape(st.mass_kg, 1, []);
            o = obj.options;
            switch p
                case "cdScale"
                    dF = st.Fdrag / o.drag.cdScale;                 % the drag force is linear in it
                case "crScale"
                    dF = st.Frad / o.srp.crScale;                   % as is the radiation force
                case "alpha"
                    if o.drag.gsi == "constantCd"
                        dF = zeros(size(st.Fdrag));                 % constantCd has no accommodation coefficient
                    else
                        shadowArgs = {};
                        if ~isempty(st.flow), shadowArgs = {'Shadow', st.flow}; end
                        extra = {};
                        if o.drag.gsi == "maxwell", extra = {'SpecularFraction', o.drag.specularFraction}; end
                        dF = o.drag.cdScale * NonGravForces.kernel.partials("drag", "alpha", obj.macro, ...
                            'Rho', st.rho, 'Species', st.species, 'T_K', st.Tgas, 'Vrel_body', st.vrelBody, ...
                            'Model', o.drag.gsi, 'Alpha', st.acc, 'Tw_K', o.drag.Tw_K, extra{:}, shadowArgs{:});
                    end
            end
            d = NonGravForces.kernel.assemble(dF, obj.macro) ./ m;
            if frame == "gcrs", d = obj.toGcrs(st.R, d); end
        end

        function d = parameterPartialNumeric(obj, p, epoch, r, v, q, frame)
            %PARAMETERPARTIALNUMERIC Central difference in a scalar option, for cross-checking the analytic form.
            base = obj.options;
            field = struct('cdScale', "drag.cdScale", 'crScale', "srp.crScale", 'alpha', "drag.alpha");
            name = field.(p);
            parts = split(name, "."); x0 = base.(parts(1)).(parts(2));
            h = max(1e-6 * abs(x0), 1e-8);
            if p == "alpha"
                h = min(h, min(x0, 1 - x0) / 2);                  % stay inside [0, 1]: the option refuses the edge
                if h <= 0
                    NonGravForces.fail("Model", "unsupported", "alpha = %g sits on the edge of [0, 1], where no central difference exists -> use the analytic method", x0);
                end
            end
            up = obj; up.options = base.set(name, x0 + h);
            dn = obj; dn.options = base.set(name, x0 - h);
            aP = up.accel(epoch, r, v, q, 'Frame', frame);
            aM = dn.accel(epoch, r, v, q, 'Frame', frame);
            d = (aP - aM) / (2 * h);
        end
    end
end

function mass = localMass(m, macro)
%LOCALMASS A constant mass, a massModel, or the macro model's own; never an assumed spacecraft.
if isempty(m)
    if isnan(macro.mass_kg)
        NonGravForces.fail("Model", "inputMass", ...
            "the macro model '%s' carries no mass and no 'Mass' was given -> pass 'Mass', <kg> or io.massModel(...), or build the macro model with 'Mass_kg'", macro.name);
    end
    mass = localConstant(macro.mass_kg, "macro model " + macro.name);
elseif isnumeric(m) && isscalar(m) && isfinite(m) && m > 0
    mass = localConstant(m, "caller");
elseif isstruct(m) && isscalar(m) && isfield(m, 'at') && isa(m.at, 'function_handle')
    mass = m;
else
    NonGravForces.fail("Model", "inputMass", "'Mass' must be a positive mass in kg or an io.massModel result -> got a %s", class(m));
end
end

function mass = localConstant(m, source)
mass = struct('kind', "constant", 'value_kg', m, 'source', source, 'at', @(epoch) repmat(m, 1, numel(epoch.sec)));
end

function localCheckOptions(o, macro)
%LOCALCHECKOPTIONS Choices that contradict each other are refused at construction, not ignored at evaluation.
d = o.drag;
if d.gsi ~= "constantCd"
    if d.accommodation == "fixed" && isempty(d.alpha)
        NonGravForces.fail("Model", "inputOptions", "drag.accommodation fixed needs drag.alpha -> set('drag.alpha', 0.85) or choose pilinski2010 / pilinski2011 / walker2014");
    end
    if d.accommodation ~= "fixed" && ~isempty(d.alpha)
        NonGravForces.fail("Model", "inputOptions", "drag.alpha is read only by accommodation fixed, and %s computes its own -> set('drag.alpha', []) so the value is not silently ignored", d.accommodation);
    end
    if d.accommodation == "walker2014" && isempty(d.surfaceMass_amu)
        NonGravForces.fail("Model", "inputOptions", "drag.accommodation walker2014 needs drag.surfaceMass_amu -> the sources disagree (263.3223 amu Mehta et al. 2022, 65 amu Pilinski 2013); state which one you use");
    end
    if d.accommodation ~= "walker2014" && ~isempty(d.surfaceMass_amu)
        NonGravForces.fail("Model", "inputOptions", "drag.surfaceMass_amu is read only by accommodation walker2014, and %s does not (sesam fixes 65 amu) -> set('drag.surfaceMass_amu', []) so the value is not silently ignored", d.accommodation);
    end
    if d.gsi == "maxwell" && isempty(d.specularFraction)
        NonGravForces.fail("Model", "inputOptions", "drag.gsi maxwell needs drag.specularFraction -> set('drag.specularFraction', 0.1), say");
    end
end
if o.geom.flowShadow == "table" && o.geom.selfShadow == "raycast"
    NonGravForces.fail("Model", "inputOptions", "geom.flowShadow table tabulates the lut occlusion, and geom.selfShadow asks for the ray cast -> set('geom.flowShadow', ""exact"") to ray-cast the stream too, or use selfShadow lut");
end
if o.thermal.internalPower_W > 0 && any(isnan(macro.internalShare))
    NonGravForces.fail("Model", "inputOptions", "thermal.internalPower_W > 0 needs every panel's internalShare on the macro model '%s' -> give the panels internalShare, or set the power to 0", macro.name);
end
end

function localCheckDensityScale(o, S)
%LOCALCHECKDENSITYSCALE A density scale is used only when asked for, and only on the density it was formed for.
if isempty(S)
    return
end
if o.drag.densityScale == "none"
    NonGravForces.fail("Model", "inputOptions", "'DensityScale' was given and drag.densityScale is ""none"" -> set('drag.densityScale', ""swarm"", ""driver"" or ""combined"") or leave 'DensityScale' out, so it is not silently ignored");
end
if ~isstruct(S) || ~isfield(S, 'at') || ~isfield(S, 'model') || ~isfield(S, 'kind')
    NonGravForces.fail("Model", "inputOptions", "'DensityScale' must be an env.densityScale or env.densityDriverScale result -> build it with one of them");
end
if S.kind ~= o.drag.densityScale
    NonGravForces.fail("Model", "inputOptions", "the density scale is of kind %s and drag.densityScale is %s -> pass env.densityScale for swarm, env.densityDriverScale for driver", S.kind, o.drag.densityScale);
end
if S.model ~= o.drag.density
    NonGravForces.fail("Model", "inputOptions", "the density scale was formed for %s and drag.density is %s -> form it with 'Model', ""%s""", S.model, o.drag.density, o.drag.density);
end
end

function p = localStamp(p, unit)
%LOCALSTAMP The Model's own provenance record, renamed and re-timed: the option hash was formed once.
p.unit = unit;
p.created = datetime('now', 'TimeZone', 'UTC');
end

function s = localFirstEpoch(src)
%LOCALFIRSTEPOCH The source list of the first epoch only, for the steady state that starts a transient.
s = src;
s.dir = src.dir(:, :, 1); s.P = src.P(:, 1); s.nu = src.nu(1); s.K = 1;
end
