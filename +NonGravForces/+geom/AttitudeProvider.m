classdef AttitudeProvider < handle
%ATTITUDEPROVIDER BODY->GCRS quaternions at query epochs (SCA1B or nominal) via referenceSystem.interpolateQuaternions. NOT YET IMPLEMENTED (stub).
%
%   Specified in HANDOVER.md, WP-11 (level L3). The constructor raises
%   NonGravForces:AttitudeProvider:notImplemented so that a caller written against the
%   future API fails loudly (rule R4).
%
% ERRORS
%   NonGravForces:AttitudeProvider:notImplemented   Always, from the constructor.
%
% REFERENCES
%   HANDOVER.md Rev 2.0, WP-11.
%
% See also NonGravForces.Model.
%
% KEYWORDS stub WP-11
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    methods
        function obj = AttitudeProvider(varargin)
            NonGravForces.fail("AttitudeProvider", "notImplemented", "AttitudeProvider is specified in WP-11 -> not delivered in 0.1.0, see HANDOVER.md");
        end
    end
end
