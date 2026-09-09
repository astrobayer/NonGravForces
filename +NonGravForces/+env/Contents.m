% ENV  NonGravForces.env - L2 layer.
% Version 0.9.1  08-Sep-2026
%
%   density              - Thermosphere density and composition dispatcher (exponential msis00 msis21 dtm2020 jb2008)
%   wind                 - Relative-velocity builder: co-rotation plus optional HWM14 wind [stub]
%   eclipseFactor        - Shadow function nu from 0 to 1: cylindrical conical solaarsCF plus lunar
%   earthRadiance        - Earth albedo and IR radiance grids (pointSource knocke ceres) with ADM anisotropy
%   magneticField        - IGRF-14 geomagnetic field in ITRS and body frame [stub]
%   radiationSources     - Source list {dir_body P band} from Sun and Earth radiance for one epoch batch
