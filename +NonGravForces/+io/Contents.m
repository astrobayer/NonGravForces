% IO  NonGravForces.io - L1 layer.
% Version 0.20.1  09-Nov-2026
%
%   readACC1B            - ACC1B/ACT1B/ACH1B reader over gracefoReading (3xK SRF accelerations + aux)
%   readSCA1B            - SCA1B reader: ICRF->SRF quaternions conjugated once to BODY->GCRS
%   readGNV1B            - GNV1B reader: ITRS state rotated to GCRS via referenceSystem
%   readTHR1B            - THR1B thruster event reader (impulse priors for qc)
%   readMAS1B            - MAS1B tank-mass reader (NOT spacecraft mass)
%   readMAG1B            - MAG1B reader: measured field (SI) and magnetorquer branch and net currents
%   massModel            - Spacecraft mass model: handbook constant or MAS1B-driven variable mass
%   spaceWeather         - Space-weather provider with query(epoch): F10.7 F30 S/M/Y10.7 Ap Kp Dst
%   importCERES          - CERES EBAF TOA monthly fluxes (obs4MIPs netCDF) as lat x lon x month grids
%   importGroopsGrid     - GROOPS binary gridded data: grid, solid angles and values
%   igrfCoefficients     - IGRF-14 coefficient loader with the model's own interpolation and SV rule
%   readSwarmDensity     - Swarm POD-only and TOLEOS accelerometer thermosphere densities from VirES HAPI CSV
%   readIMU1B            - IMU1B reader: filtered accumulated gyro angles (rad), one record per gyro sample
