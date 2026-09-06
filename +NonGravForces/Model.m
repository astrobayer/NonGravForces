classdef Model
%MODEL Single point of access: non-gravitational forces, torques and partials.
%
%   M = NonGravForces.Model(macro, options)
%   M = NonGravForces.Model(macro, options, 'Attitude', att, 'Mass', massModel)
%
%   [a, info]   = M.accel(epoch, r_gcrs, v_gcrs, q_body2gcrs)   % 3xK m/s^2
%   [tau, info] = M.torque(epoch, r, v, q)                      % 3xK N m body
%   T           = M.series(epoch, r, v, q)                      % timetable
%   [J, info]   = M.partials(epoch, r, v, q, 'WithRespectTo', ["r","v","cd","cr"])
%   f           = M.rhs()        % a(t,r,v) handle, non-grav only, for gjToolbox
%   B           = M.budget(...)  % per-effect, per-axis RMS table
%
%   The Model owns: the packed macro model, the attitude provider, the mass
%   model, the resolved Options, and every cache (self-shadow LUT, density
%   tables, radiance grids). It computes forces and torques and nothing else:
%   time, frames and attitude kinematics are referenceSystem, reading is
%   gracefoReading, integration is gjToolbox, estimation is lsqvce.
%
%   ONE KERNEL (principle P2). accel and torque are two reductions of the same
%   per-panel force array F_k: F = sum_k F_k, tau = sum_k (r_k - r_cm) x F_k.
%   series is accel+torque over chunks, optionally under parfor, with results
%   required to be bit-identical regardless of chunking or parallelism.
%
% CONSTRUCTOR INPUTS
%   Name     Size  Type                     Units  Description
%   ----     ----  ----                     -----  -----------
%   macro    1x1   NonGravForces.geom.MacroModel  -  Packed panel geometry
%   options  1x1   NonGravForces.Options    -      Resolved options
%
% NAME-VALUE
%   Attitude 1x1   NonGravForces.geom.AttitudeProvider  BODY->GCRS quaternions
%   Mass     1x1   struct (io.massModel)    kg     Mass model, default handbook
%
% ERRORS
%   NonGravForces:Model:inputMacro      macro is not a MacroModel
%   NonGravForces:Model:inputOptions    options is not an Options object
%   NonGravForces:Model:notImplemented  evaluation methods, until LP-6 lands
%
% ACCURACY
%   See the ACCURACY block of each kernel and environment unit; the Model
%   adds no approximation of its own. Provenance (info.provenance) records
%   the resolved options and backends actually used.
%
% EXAMPLE
%   macro = NonGravForces.geom.graceFO("C");
%   M = NonGravForces.Model(macro, NonGravForces.Options("standard"));
%
% See also NonGravForces.Options, NonGravForces.geom.MacroModel.
%
% KEYWORDS model accel torque series partials single point of access
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    properties (SetAccess = private)
        macro
        options (1,1) NonGravForces.Options
        attitude
        mass struct
        provenance struct
    end
    methods
        function obj = Model(macro, options, nv)
            arguments
                macro
                options (1,1) NonGravForces.Options = NonGravForces.Options("standard")
                nv.Attitude = []
                nv.Mass struct = struct('kind', "constant", 'value_kg', NaN)
            end
            if ~isa(macro, 'NonGravForces.geom.MacroModel')
                NonGravForces.fail("Model", "inputMacro", ...
                    "macro must be a NonGravForces.geom.MacroModel -> build one with geom.graceFO or geom.boxWing");
            end
            obj.macro = macro;
            obj.options = options;
            obj.attitude = nv.Attitude;
            obj.mass = nv.Mass;
            obj.provenance = NonGravForces.provenance("Model", 'Options', options.toStruct());
        end
        function varargout = accel(~, varargin)
            varargout = cell(1, nargout);
            NonGravForces.fail("Model", "notImplemented", "accel is specified in LP-6 (WP-14) -> not delivered in 0.1.0");
        end
        function varargout = torque(~, varargin)
            varargout = cell(1, nargout);
            NonGravForces.fail("Model", "notImplemented", "torque is specified in LP-6 (WP-14) -> not delivered in 0.1.0");
        end
        function varargout = series(~, varargin)
            varargout = cell(1, nargout);
            NonGravForces.fail("Model", "notImplemented", "series is specified in LP-6 (WP-15) -> not delivered in 0.1.0");
        end
        function varargout = partials(~, varargin)
            varargout = cell(1, nargout);
            NonGravForces.fail("Model", "notImplemented", "partials is specified in LP-6 (WP-16) -> not delivered in 0.1.0");
        end
        function varargout = rhs(~, varargin)
            varargout = cell(1, nargout);
            NonGravForces.fail("Model", "notImplemented", "rhs is specified in LP-6 (WP-17) -> not delivered in 0.1.0");
        end
        function varargout = budget(~, varargin)
            varargout = cell(1, nargout);
            NonGravForces.fail("Model", "notImplemented", "budget is specified in LP-8 (WP-30) -> not delivered in 0.1.0");
        end
    end
end
