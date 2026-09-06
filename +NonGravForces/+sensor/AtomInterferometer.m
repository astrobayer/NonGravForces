classdef AtomInterferometer < NonGravForces.sensor.Sensor
%ATOMINTERFEROMETER Cold-atom interferometer: phase = k_eff a T^2 with dead time and rotation. NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-29 (level L6). The constructor raises
%   NonGravForces:AtomInterferometer:notImplemented so that a caller written against the
%   future API fails loudly (rule R4).
%
% ERRORS
%   NonGravForces:AtomInterferometer:notImplemented   Always, from the constructor.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-29.
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-29
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    methods
        function obj = AtomInterferometer(varargin)
            NonGravForces.fail("AtomInterferometer", "notImplemented", "AtomInterferometer is specified in WP-29 -> not delivered in 0.1.0, see HANDOVER.md");
        end
        function y = observe(~, ~), y = []; end
        function [J, names] = partials(~, ~), J = []; names = string.empty; end
        function S = noiseModel(~, ~), S = []; end
        function n = numParameters(~), n = 0; end
    end
end
