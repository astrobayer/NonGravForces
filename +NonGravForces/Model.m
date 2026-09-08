classdef Model
%MODEL Single point of access: non-gravitational forces, torques and partials.
%
%   M = NonGravForces.Model(macro, options, 'Eop', eop, 'SpaceWeather', sw)
%   M = NonGravForces.Model(macro, options, 'Attitude', att, 'Mass', massModel, 'Eop', eop, 'SpaceWeather', sw)
%
%   [a, info]      = M.accel(epoch, r_gcrs, v_gcrs, q_body2gcrs)   % 3xK m/s^2, body frame
%   [tau, info]    = M.torque(epoch, r, v, q)                      % 3xK N m, body frame
%   T              = M.series(epoch, r, v, q)                      % timetable over chunks
%   [J, names, info] = M.partials(epoch, r, v, q, 'WithRespectTo', ["r","v","cdScale"])
%   f              = M.rhs()                                       % handle for an integrator
%   B              = M.budget(epoch, r, v, q)                      % per-effect table
%
%   The Model owns the packed macro model, the attitude provider, the mass model,
%   the Earth-orientation and space-weather sources, the resolved Options and the
%   self-shadow look-up table. It computes forces and torques and nothing else:
%   time, frames and attitude kinematics are referenceSystem, reading is
%   gracefoReading, integration is gjToolbox, estimation is lsqvce.
%
%   ONE KERNEL (principle P2). accel and torque are two reductions of the same
%   per-panel force array F_k: F = sum_k F_k and tau = sum_k (c_k - r_cm) x F_k,
%   computed once per epoch by one private evaluation. budget is the same array
%   kept per effect instead of summed, and it adds back to accel bitwise.
%
%   WHAT IT COMPOSES, in order, per epoch:
%     1. the relative velocity            env.wind (co-rotation from referenceSystem)
%     2. the geodetic height              referenceSystem.cartesianToGeodetic on the ITRS position
%     3. density and composition          env.density with the space-weather drivers
%     4. the radiation source list        env.radiationSources (Sun, Earth tiles, eclipse)
%     5. self-shadowing                   geom.selfShadow for the gas stream and every source
%     6. panel temperatures               geom.thermalState
%     7. per-panel forces                 kernel.dragPanelForce, radiationPanelForce, thermalRecoil
%     8. reductions                       kernel.assemble, and kernel.bodyTorques for the torque
%   Every step is the unit that owns that fact; the Model adds no physics.
%
%   FRAME. Results are BODY frame by default (principle P4: the body frame is
%   what an accelerometer reads and what the sensor layer maps from). Pass
%   'Frame', "gcrs" for the inertial frame, which is what an integrator wants.
%
% CONSTRUCTOR INPUTS
%   Name     Size  Type                            Units  Description
%   ----     ----  ----                            -----  -----------
%   macro    1x1   NonGravForces.geom.MacroModel   -      Packed panel geometry
%   options  1x1   NonGravForces.Options           -      Resolved options
%
% CONSTRUCTOR NAME-VALUE ARGUMENTS
%   Name          Size  Type    Units  Description
%   ----          ----  ----    -----  -----------
%   Attitude      1x1   geom.AttitudeProvider  -  BODY->GCRS quaternions, used when a
%                                          method is called without q. Default [].
%   Mass          1x1   struct  kg     io.massModel result. Default: the handbook
%                                      constant for GRACE-C.
%   Eop           1x1   struct  -      referenceSystem.readEop result or a file path.
%                                      REQUIRED to evaluate: the geodetic height needs
%                                      the ITRS rotation and the co-rotation needs its
%                                      rate. No default, because a wrong Earth
%                                      orientation is invisible in the answer (R5).
%   SpaceWeather  1x1   struct  -      io.spaceWeather result. Required by every
%                                      density backend but "exponential".
%   Ephemeris     1x1   struct  -      referenceSystem.readEphemeris result, required
%                                      when sun.source is "rs.de440" (O-7).
%   Exponential   1x1   struct  -      .rho0_kg_m3 .h0_m .H_m, required when
%                                      drag.density is "exponential".
%
% OUTPUT (accel / torque)
%   a, tau   3xK   double         acceleration [m/s^2] or torque [N m]
%   info     1x1   struct         .frame, .drag .radiation .thermal (3xK each, the
%                                 same frame), .bodyTorque (torque only), .rho (1xK),
%                                 .nu (1xK shadow factor), .mass_kg (1xK), .vrel_body,
%                                 .sources (M), .temperature_K (NxK), .shadowDrag,
%                                 .backends (struct of the units actually used),
%                                 .provenance
%
% ACCURACY
%   The Model adds no approximation: every number comes from the unit that owns
%   it, and each of those states its own error. What the Model does own is the
%   ORDER, which is fixed above, and the resolved options, which travel in
%   info.provenance (rule R2). Measured end to end on the real 2020-01-01
%   GRACE-C day, standard preset, 523 km: drag 8.5e-9, radiation 2.6e-8 and
%   thermal 6.1e-9 m/s^2, the last against 6.4 nm/s^2 cross-track RMS published
%   for GRACE (Siemes et al. 2023). The dominant error is the thermospheric
%   density, 10-30 % (HANDOVER Part E), not anything here.
%
% ERRORS
%   NonGravForces:Model:inputMacro      macro is not a MacroModel
%   NonGravForces:Model:inputOptions    options is not an Options object
%   NonGravForces:Model:inputState      epoch, r, v or q disagree in size, or are not finite
%   NonGravForces:Model:inputAttitude   no q was given and the Model has no attitude provider
%   NonGravForces:Model:dependEop       evaluated without 'Eop'
%   NonGravForces:Model:dependWeather   evaluated without 'SpaceWeather' for a driven backend
%   NonGravForces:Model:unsupported     an option combination this level does not deliver
%   NonGravForces:<unit>:*              propagated unwrapped from the unit that refused
%
% REFERENCES
%   HANDOVER Part D.3 (workflow) and Part G, WP-14 to WP-17; theory sections 2, 3,
%   5 and 7; Siemes C. et al. (2023) Adv. Space Res.
%
% EXAMPLE
%   eop = referenceSystem.readEop(eopFile);
%   sw  = NonGravForces.io.spaceWeather(swFile);
%   [ep, r, v] = NonGravForces.io.readGNV1B(gnv, 'Eop', eop);
%   [~, q]     = NonGravForces.io.readSCA1B(sca);
%   M = NonGravForces.Model(NonGravForces.geom.graceFO("C"), NonGravForces.Options("standard"), ...
%                           'Eop', eop, 'SpaceWeather', sw);
%   [a, info] = M.accel(ep, r, v, q);
%   vecnorm(info.drag)          % ~1e-8 m/s^2 at solar minimum
%
% See also NonGravForces.Options, NonGravForces.geom.MacroModel, NonGravForces.env.wind, NonGravForces.kernel.assemble.
%
% KEYWORDS model accel torque series partials budget single point of access
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    properties (SetAccess = private)
        macro
        options (1,1) NonGravForces.Options
        attitude
        mass struct
        eop
        spaceWeather
        ephemeris
        exponential
        lut struct
        provenance struct
    end
    methods
        function obj = Model(macro, options, nv)
            arguments
                macro
                options = NonGravForces.Options("standard")
                nv.Attitude = []
                nv.Mass struct = struct('kind', "constant", 'value_kg', NaN)
                nv.Eop = []
                nv.SpaceWeather = []
                nv.Ephemeris = []
                nv.Exponential = []
            end
            if ~isa(options, 'NonGravForces.Options') || ~isscalar(options)
                NonGravForces.fail("Model", "inputOptions", ...
                    "options must be a scalar NonGravForces.Options -> build one with NonGravForces.Options(preset)");
            end
            if ~isa(macro, 'NonGravForces.geom.MacroModel')
                NonGravForces.fail("Model", "inputMacro", ...
                    "macro must be a NonGravForces.geom.MacroModel -> build one with geom.graceFO or geom.boxWing");
            end
            obj.macro = macro;
            obj.options = options;
            obj.attitude = nv.Attitude;
            if isnan(nv.Mass.value_kg) && nv.Mass.kind == "constant"
                obj.mass = NonGravForces.io.massModel("C");     % the handbook constant, stated not assumed
            else
                obj.mass = nv.Mass;
            end
            obj.eop = nv.Eop;
            obj.spaceWeather = nv.SpaceWeather;
            obj.ephemeris = nv.Ephemeris;
            obj.exponential = nv.Exponential;
            % The self-shadow look-up table costs a ray cast per icosphere vertex, so it is built
            % once here rather than on every evaluation (measured 2026-09-08: rebuilding it made a
            % 5-epoch accel take 1.7 s).
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
            end
            st = obj.evaluate(epoch, r, v, q);
            m = reshape(st.mass_kg, 1, []);
            aD = NonGravForces.kernel.assemble(st.Fdrag, obj.macro) ./ m;
            aR = NonGravForces.kernel.assemble(st.Frad, obj.macro) ./ m;
            aT = NonGravForces.kernel.assemble(st.Fthermal, obj.macro) ./ m;
            a = aD + aR + aT;
            if nv.Frame == "gcrs"
                [a, aD, aR, aT] = obj.toGcrs(st.R, a, aD, aR, aT);
            end
            info = obj.stateInfo(st, nv.Frame);
            info.drag = aD; info.radiation = aR; info.thermal = aT;
            info.provenance = NonGravForces.provenance("Model.accel", 'Options', obj.options.toStruct());
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
            st = obj.evaluate(epoch, r, v, q);
            [~, tD] = NonGravForces.kernel.assemble(st.Fdrag, obj.macro);
            [~, tR] = NonGravForces.kernel.assemble(st.Frad, obj.macro);
            [~, tT] = NonGravForces.kernel.assemble(st.Fthermal, obj.macro);
            o = obj.options.torques;
            if o.gravityGradient || o.magnetic
                rBody = reshape(pagemtimes(pagetranspose(st.R), reshape(st.r, 3, 1, [])), 3, []);
                tB = NonGravForces.kernel.bodyTorques(obj.macro, rBody, ...
                    'GravityGradient', o.gravityGradient, 'Magnetic', o.magnetic, ...
                    'B_body', nv.B_body, 'Dipole_Am2', nv.Dipole_Am2);
            else
                tB = zeros(3, st.K);
            end
            tau = tD + tR + tT + tB;
            if nv.Frame == "gcrs"
                [tau, tD, tR, tT, tB] = obj.toGcrs(st.R, tau, tD, tR, tT, tB);
            end
            info = obj.stateInfo(st, nv.Frame);
            info.drag = tD; info.radiation = tR; info.thermal = tT; info.bodyTorque = tB;
            info.provenance = NonGravForces.provenance("Model.torque", 'Options', obj.options.toStruct());
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
            if obj.options.series.useParallel
                parfor c = 1:nC
                    [aC{c}, tC{c}] = obj.chunkEval(epoch, r, v, q, first(c), last(c), wantA, wantT, nv.Frame); %#ok<PFBNS>
                end
            else
                for c = 1:nC
                    [aC{c}, tC{c}] = obj.chunkEval(epoch, r, v, q, first(c), last(c), wantA, wantT, nv.Frame);
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
            if method == "complexStep"
                NonGravForces.fail("Model", "unsupported", ...
                    "complexStep cannot differentiate this chain -> MATLAB's erf and erfc raise notFullReal on a complex argument and abs/vecnorm discard the imaginary part, so no complex step survives the gas-surface coefficients or the speed; use analytic or finiteDiff (HANDOVER, 2026-09-08)");
            end
            K = numel(epoch.sec);
            obj.checkState(epoch, r, v, q, K);
            names = string.empty(1, 0);
            cols = {};
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
                        cols{end+1} = obj.parameterPartial(p, epoch, r, v, q, nv.Frame, method); %#ok<AGROW>
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

        function f = rhs(obj)
            %RHS Function handle of the non-gravitational acceleration in GCRS, for an integrator.
            %   f(epoch, r, v, q) -> 3xK m/s^2 in GCRS. The Model supplies non-gravitational
            %   forces only: the caller's integrator adds gravity and everything else.
            %   The gjToolbox.Dynamics contract could not be confirmed on this host,
            %   which has no gjToolbox (HANDOVER O-4), so the shape of this handle is
            %   what Part B.3 records and nothing here proves the consumer accepts it.
            f = @(epoch, r, v, q) obj.accel(epoch, r, v, q, 'Frame', "gcrs");
        end
    end

    methods (Access = private)
        function st = evaluate(obj, epoch, r, v, q)
            %EVALUATE The one pass every method reduces from: environment, geometry, per-panel forces.
            K = numel(epoch.sec);
            q = obj.resolveAttitude(epoch, q, K);
            obj.checkState(epoch, r, v, q, K);
            o = obj.options;
            if isempty(obj.eop)
                NonGravForces.fail("Model", "dependEop", "evaluating needs 'Eop' at construction -> the geodetic height and the co-rotation both come from the ITRS rotation, and neither may be guessed");
            end
            % 1. relative velocity, GCRS then body
            [vrel, windInfo] = NonGravForces.env.wind(o.drag.wind, epoch, r, v, 'Eop', obj.eop);
            R = referenceSystem.rotationFromQuaternion(q);                 % body -> GCRS
            vrelBody = reshape(pagemtimes(pagetranspose(R), reshape(vrel, 3, 1, K)), 3, K);
            % 2. geodetic height from the ITRS position
            [Ritrs, ~] = referenceSystem.rotationItrsToGcrs(epoch, 'Eop', obj.eop);
            rItrs = reshape(pagemtimes(pagetranspose(Ritrs), reshape(r, 3, 1, K)), 3, K);
            lla = referenceSystem.cartesianToGeodetic(rItrs', 'Ellipsoid', "WGS84")';   % the dependency is Nx3
            % 3. density and composition
            [rho, species, Tgas] = obj.callDensity(epoch, lla);
            % 4. radiation sources
            srcArgs = {};
            if ~isempty(obj.ephemeris), srcArgs = {'Ephemeris', obj.ephemeris}; end
            [src, srcInfo] = NonGravForces.env.radiationSources(epoch, r, q, o, srcArgs{:});
            % 5. self-shadowing
            [fDrag, fRad] = obj.shadowFractions(vrelBody, src, K);
            % 6. panel temperatures
            if o.thermal.model == "none"
                Tpanel = zeros(obj.macro.n, K);
            else
                % geom.thermalState takes ONE exposure per panel and epoch, while self-shadowing is
                % per source. The Sun carries the great majority of the absorbed flux, so its column
                % is the one passed; the Earth tiles' self-shadowing does not reach the panel
                % temperature (HANDOVER O-27).
                Tpanel = NonGravForces.geom.thermalState(obj.macro, src, 'Model', o.thermal.model, ...
                    'InternalPower_W', o.thermal.internalPower_W, 'Exposed', fRad(:, :, 1));
            end
            % 7. per-panel forces
            st.Fdrag = o.drag.cdScale * obj.dragForces(rho, species, Tgas, vrelBody, fDrag);
            st.Frad = o.srp.crScale * NonGravForces.kernel.radiationPanelForce(obj.macro, src, 'Shadow', fRad);
            if o.thermal.model == "none"
                st.Fthermal = zeros(3, obj.macro.n, K);
            else
                st.Fthermal = NonGravForces.kernel.thermalRecoil(obj.macro, Tpanel);
            end
            st.K = K; st.R = R; st.r = r; st.rho = rho; st.species = species; st.Tgas = Tgas;
            st.vrelBody = vrelBody; st.lla = lla; st.src = src; st.nu = src.nu;
            st.shadowDrag = fDrag; st.shadowRad = fRad; st.Tpanel = Tpanel;
            st.mass_kg = obj.mass.at(epoch);
            st.backends = struct('density', o.drag.density, 'wind', windInfo.model, 'gsi', o.drag.gsi, ...
                'eclipse', o.srp.eclipse, 'erp', o.erp.backend, 'sun', srcInfo.sunSource, ...
                'thermal', o.thermal.model, 'selfShadow', o.geom.selfShadow);
        end

        function [rho, species, Tgas] = callDensity(obj, epoch, lla)
            %CALLDENSITY The density backend with the drivers the Model was given, never a default set.
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
        end

        function F = dragForces(obj, rho, species, Tgas, vrelBody, fDrag)
            %DRAGFORCES The gas-surface model the options name, with the arguments that model needs.
            o = obj.options;
            switch o.drag.gsi
                case "constantCd"
                    F = NonGravForces.kernel.dragPanelForce(obj.macro, rho, [], [], vrelBody, ...
                        'Model', "constantCd", 'Cd', o.drag.cd, 'Shadow', fDrag);
                case "maxwell"
                    NonGravForces.fail("Model", "unsupported", "drag.gsi maxwell needs a specular fraction, which Options does not carry -> use sentman or constantCd; adding the fraction to Options would be a second authority for a number no preset sets");
                otherwise
                    % sentman here, and cll / schamberg refuse inside gsiCoefficients (O-17)
                    alpha = NonGravForces.kernel.accommodation(o.drag.accommodation, species, Tgas, 'Alpha', o.drag.alpha);
                    F = NonGravForces.kernel.dragPanelForce(obj.macro, rho, species, Tgas, vrelBody, ...
                        'Model', o.drag.gsi, 'Alpha', alpha, 'Tw_K', o.drag.Tw_K, 'Shadow', fDrag);
            end
        end

        function [fDrag, fRad] = shadowFractions(obj, vrelBody, src, K)
            %SHADOWFRACTIONS Self-shadowing for the gas stream and for every radiation source.
            N = obj.macro.n; M = size(src.dir, 2);
            if obj.options.geom.selfShadow == "none"
                fDrag = ones(N, K); fRad = ones(N, K, M);
                return
            end
            method = obj.options.geom.selfShadow;
            args = {'Method', method};
            if method == "lut", args = [args, {'Lut', obj.lut.table}]; end   % built once at construction
            dirGas = vrelBody ./ vecnorm(vrelBody);                 % the gas arrives from +v_rel
            fDrag = NonGravForces.geom.selfShadow(obj.macro, dirGas, args{:});
            dirAll = reshape(permute(src.dir, [1 3 2]), 3, K * M);  % every source direction of every epoch
            fFlat = NonGravForces.geom.selfShadow(obj.macro, dirAll, args{:});
            fRad = reshape(fFlat, N, K, M);
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

        function info = stateInfo(obj, st, frame) %#ok<INUSL>
            info = struct('frame', frame, 'K', st.K, 'rho', st.rho, 'nu', st.nu, ...
                'mass_kg', st.mass_kg, 'vrel_body', st.vrelBody, 'sources', size(st.src.dir, 2), ...
                'temperature_K', st.Tpanel, 'shadowDrag', st.shadowDrag, 'lla', st.lla, ...
                'backends', st.backends);
        end

        function [a, tq] = chunkEval(obj, epoch, r, v, q, i0, i1, wantA, wantT, frame)
            %CHUNKEVAL One chunk of series; the slice must give the same numbers as the whole.
            idx = i0:i1;
            ep = referenceSystem.daySecondToEpoch([epoch.day(idx) epoch.sec(idx)], epoch.scale);
            qc = []; if ~isempty(q), qc = q(:, idx); end
            a = []; tq = [];
            if wantA, a = obj.accel(ep, r(:, idx), v(:, idx), qc, 'Frame', frame); end
            if wantT, tq = obj.torque(ep, r(:, idx), v(:, idx), qc, 'Frame', frame); end
        end

        function d = parameterPartial(obj, p, epoch, r, v, q, frame, method)
            %PARAMETERPARTIAL Analytic derivatives of the acceleration with respect to a scalar parameter.
            if method == "finiteDiff"
                d = obj.parameterPartialNumeric(p, epoch, r, v, q, frame);
                return
            end
            st = obj.evaluate(epoch, r, v, q);
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
                        alpha = NonGravForces.kernel.accommodation(o.drag.accommodation, st.species, st.Tgas, 'Alpha', o.drag.alpha);
                        dF = o.drag.cdScale * NonGravForces.kernel.partials("drag", "alpha", obj.macro, ...
                            'Rho', st.rho, 'Species', st.species, 'T_K', st.Tgas, 'Vrel_body', st.vrelBody, ...
                            'Model', o.drag.gsi, 'Alpha', alpha, 'Tw_K', o.drag.Tw_K, 'Shadow', st.shadowDrag);
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
            up = obj; up.options = base.set(name, x0 + h);
            dn = obj; dn.options = base.set(name, x0 - h);
            aP = up.accel(epoch, r, v, q, 'Frame', frame);
            aM = dn.accel(epoch, r, v, q, 'Frame', frame);
            d = (aP - aM) / (2 * h);
        end
    end
end
