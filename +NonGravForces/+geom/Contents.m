% GEOM  NonGravForces.geom - L3 layer.
% Version 0.15.0  15-Oct-2026
%
%   MacroModel           - Panel geometry frozen and validated once
%   graceFO              - GRACE-FO macro model from JPL D-56935
%   grace                - GRACE-A/B macro model from the GRACE Product Specification (Bettadpur 2012)
%   boxWing              - Generic box-wing macro model factory
%   selfShadow           - Per-panel unobstructed fraction toward a direction: exact ray cast or cube-map table
%   flowShadow           - Self-shadowing of a free-molecular gas stream: flux-weighted occlusion ratios per species
%   thermalState         - Per-panel temperatures: steady state or transient (Wöske)
%   AttitudeProvider     - BODY->GCRS quaternions at query epochs (SCA1B or nominal) via referenceSystem.interpolateQuaternions
%   checkAttitudeSense   - Physics test of quaternion sense (nadir and line-of-sight)
