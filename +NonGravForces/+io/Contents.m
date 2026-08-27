% IO  NonGravForces.io - L1 layer.
% Version 2.0.0  27-Aug-2026
%
%   readACC1B            - ACC1B/ACT1B/ACH1B reader over gracefoReading (3xK SRF accelerations + aux) [stub]
%   readSCA1B            - SCA1B reader: ICRF->SRF quaternions conjugated once to BODY->GCRS [stub]
%   readGNV1B            - GNV1B reader: ITRS state rotated to GCRS via referenceSystem [stub]
%   readTHR1B            - THR1B thruster event reader (impulse priors for qc) [stub]
%   readMAS1B            - MAS1B tank-mass reader (NOT spacecraft mass) [stub]
%   readMAG1B            - MAG1B magnetometer/magnetorquer reader (magnetic torque priors) [stub]
%   massModel            - Spacecraft mass model: handbook constant or MAS1B-driven variable mass [stub]
%   spaceWeather         - Space-weather provider with query(epoch): F10.7 F30 S/M/Y10.7 Ap Kp Dst [stub]
%   importCERES          - CERES SYN1deg/EBAF netCDF import into radiance grids [stub]
%   igrfCoefficients     - IGRF-14 coefficient loader [stub]
