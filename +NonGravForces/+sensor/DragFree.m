classdef DragFree < NonGravForces.sensor.Sensor
%DRAGFREE Drag-free residual acceleration model. NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-29 (level L6). The constructor raises
%   NonGravForces:DragFree:notImplemented so that a caller written against the
%   future API fails loudly (rule R4).
%
% ERRORS
%   NonGravForces:DragFree:notImplemented   Always, from the constructor.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-29.
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-29
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    methods
        function obj = DragFree(varargin) %#ok<INUSD>
            NonGravForces.fail("DragFree", "notImplemented", "DragFree is specified in WP-29 -> not delivered in 2.0.0, see HANDOVER.md");
        end
        function y = observe(obj, dyn), y = []; end %#ok<INUSD>
        function [J, names] = partials(obj, dyn), J = []; names = string.empty; end %#ok<INUSD>
        function S = noiseModel(obj, f_Hz), S = []; end %#ok<INUSD>
        function n = numParameters(obj), n = 0; end %#ok<INUSD>
    end
end
