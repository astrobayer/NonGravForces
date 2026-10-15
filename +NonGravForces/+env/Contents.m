% ENV  NonGravForces.env - L2 layer.
% Version 0.15.0  15-Oct-2026
%
%   density              - Thermosphere density and composition (exponential msis00 msis21)
%   wind                 - Relative-velocity builder: co-rotation from referenceSystem plus HWM14 or a given wind field
%   eclipseFactor        - Shadow function nu from 0 to 1: cylindrical conical solaarsCF plus lunar
%   earthRadiance        - Earth albedo and longwave emission on cap tiles (pointSource knocke ceres)
%   magneticField        - IGRF main geomagnetic field at Earth-fixed positions (tesla; any frame)
%   radiationSources     - Source list {dir_body P band} from Sun and Earth radiance for one epoch batch
