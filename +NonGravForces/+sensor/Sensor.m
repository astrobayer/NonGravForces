classdef (Abstract) Sensor
%SENSOR Abstract measurement model: maps true non-gravitational dynamics to an observable.
%
%   The physics core (NonGravForces.Model) produces the TRUE quantities:
%   a_ng (3xK, body frame), tau (3xK), and - from the attitude provider via
%   referenceSystem - omega and omegaDot (3xK). A Sensor maps these to what a
%   particular instrument reads, and provides the partials of that reading
%   with respect to its own calibration parameters. Calibration (app.calibrate)
%   and quality control (qc.screen) are written against THIS interface, so a
%   new instrument (optical, MEMS, cold-atom, drag-free residual) is one new
%   subclass and zero changes elsewhere.
%
% ABSTRACT METHODS (every subclass implements)
%   y          = observe(obj, dyn)           dyn: struct .a .tau .omega .omegaDot
%                                            .epoch; y: struct .value (MxK)
%                                            .epoch .frame .mask (logical MxK)
%   [J, names] = partials(obj, dyn)          d(y)/d(theta), M*K x P, theta named
%   S          = noiseModel(obj, f_Hz)       one-sided ASD [unit/sqrt(Hz)] at f
%   n          = numParameters(obj)
%
% PROPERTIES
%   name        string   instrument label
%   axes        string   e.g. ["x","y","z"] in the sensor frame
%   sampleRate  double   Hz
%
% ERRORS
%   Subclasses use NonGravForces.fail("<Subclass>", ...).
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-29; Zahzam et al. 2022 (hybrid CAI/electrostatic).
%
% See also NonGravForces.sensor.Electrostatic, NonGravForces.sensor.AtomInterferometer.
%
% KEYWORDS sensor measurement model interface quantum accelerometer
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    properties (SetAccess = protected)
        name (1,1) string = ""
        axes (1,:) string = ["x","y","z"]
        sampleRate (1,1) double = 1
    end
    methods (Abstract)
        y = observe(obj, dyn)
        [J, names] = partials(obj, dyn)
        S = noiseModel(obj, f_Hz)
        n = numParameters(obj)
    end
end
