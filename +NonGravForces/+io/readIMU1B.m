function [epoch, angle_rad, aux] = readIMU1B(files, options)
%READIMU1B IMU1B reader: the filtered, accumulated angle of each GRACE-FO gyro, in radians, one record per sample.
%
%   [epoch, angle_rad, aux] = NonGravForces.io.readIMU1B(file)
%   [epoch, angle_rad, aux] = NonGravForces.io.readIMU1B(files, 'TimeRange', tr, 'QualityBitMask', bits)
%
%   Reads GRACE-FO Level-1B IMU data through gracefoReading. The product is NOT a
%   rate: each record carries one gyro's FiltAng, the angle that gyro has
%   ACCUMULATED about its own sensitive axis, passed through the instrument's
%   2nd-order low-pass filter, in degrees (JPL D-56935; the RL04 header's own
%   description). The three gyros of 2020-01-01 GRACE-C (gyro_id 1..3) are
%   sampled in turn, about 1 ms apart, at 8 Hz each, so a record's epoch belongs
%   to its gyro alone and the reader does NOT merge them into 3xK rows: that
%   would invent common epochs the instrument does not have. Rates are the
%   caller's (geom.attitudeQuality differences the angles per gyro).
%
%   The gyro axes are NOT the SRF axes (on 2020-01-01 C gyros 1 and 3 both
%   accumulate the orbital pitch rate, with opposite signs, gyro 2 almost
%   nothing), and no alignment table is applied here: the handbook's is not
%   held by this repository, and geom.attitudeQuality fits each gyro against the
%   star-camera rates instead of assuming an axis (on an Earth-pointing arc that
%   fit determines a residual, not the axis; see its help).
%
%   Angles are returned as published, converted to radians: a reset or a wrap
%   of the accumulator (the header's valid range is +-5758 deg) shows as a jump
%   and is left to the caller, which must not difference across it.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One or more IMU1B files of ONE satellite, chronological; .gz accepted.
%
% NAME-VALUE ARGUMENTS
%   Name            Size  Type              Units  Description
%   ----            ----  ----              -----  -----------
%   TimeRange       1x2   double|datetime   -      Forwarded to readL1b: decimal MJD or GPS datetime. Default [].
%   QualityBitMask  1xQ   double            -      qualflg bits removing a record; forwarded. Default [].
%
% OUTPUT ARGUMENTS
%   Name       Size  Type    Units  Description
%   ----       ----  ----    -----  -----------
%   epoch      1x1   struct  -      referenceSystem epoch struct, K entries, "GPS" (each record's own time)
%   angle_rad  1xK   double  rad    FiltAng of the record's gyro, deg2rad of the published value
%   aux        1x1   struct  -      .gyroId (1xK), .gyros (sorted ids present), .timeRef (1xK string,
%                                   "G" = GPS), .angleRange_rad ([-1 1]*deg2rad(5758), the header's
%                                   valid range), .product, .satellite, .files, .numRecords,
%                                   .sampling_s (median spacing of ALL records, i.e. the gyro
%                                   interleave, not one gyro's 0.125 s), .obsNames, .units,
%                                   .qualflgBits, .provenance
%
% ACCURACY
%   The conversion is deg2rad of the published decimal value (15 significant digits, so
%   1e-14 deg on a 5000 deg accumulator). Nothing is interpolated, merged or filtered here.
%
% ERRORS
%   NonGravForces:readIMU1B:dataProduct     a file is not an IMU1B product
%   NonGravForces:readIMU1B:dataSatellite   the files mix GRACE-C and GRACE-D
%   NonGravForces:readIMU1B:dataEmpty       no record survives the filters
%   NonGravForces:readIMU1B:dataColumn      FiltAng or gyro_id is missing
%   NonGravForces:readIMU1B:dataTimeRef     a record is not tagged in GPS time (time_ref ~= "G")
%   gracefoReading:*, referenceSystem:*     propagated unwrapped
%
% REFERENCES
%   JPL D-56935 (GRACE-FO Level-1 Data Product User Handbook), IMU1B; the RL04 header of
%   IMU1B_2020-01-01_C_04 (variables rcvtime_intg, rcvtime_frac, time_ref, GRACEFO_id,
%   gyro_id, FiltAng, qualflg; "8-Hz IMU data").
%
% EXAMPLE
%   root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
%   [ep, a, aux] = NonGravForces.io.readIMU1B(fullfile(root, 'data', 'samples', 'multiday', 'IMU1B_2020-01-01_C_04.txt'));
%   g1 = aux.gyroId == 1;                       % one gyro's accumulated angle
%
% See also NonGravForces.geom.attitudeQuality, NonGravForces.io.readSCA1B, gracefoReading.readL1b.
%
% KEYWORDS IMU1B gyro angle rate attitude inertial measurement unit
% Created 2026-10-30 by Claude (Anthropic) for Matthias Weigelt.
arguments
    files (1,:) string {mustBeNonempty}
    options.TimeRange = []
    options.QualityBitMask (1,:) double = []
end
unit = "readIMU1B";
[epoch, obs, raw, aux] = l1bLoad(unit, files, "IMU1B", options);
iA = l1bColumns(unit, raw, "FiltAng");
if ~isfield(raw, 'gyro_id')
    NonGravForces.fail(unit, "dataColumn", "gyro_id is missing from %s -> the file is not an RL04 IMU1B product", strjoin(files, ", "));
end
angle_rad = deg2rad(obs(:, iA)');
aux.gyroId = double(raw.gyro_id(:)');
aux.gyros = unique(aux.gyroId);
aux.timeRef = strings(1, 0);
if isfield(raw, 'time_ref')
    tr = string(raw.time_ref);
    aux.timeRef = reshape(tr, 1, []);
    if any(aux.timeRef ~= "G")
        NonGravForces.fail(unit, "dataTimeRef", "%d record(s) are not tagged in GPS time (time_ref %s) -> their epochs cannot be placed beside SCA1B's", ...
            nnz(aux.timeRef ~= "G"), strjoin(unique(aux.timeRef(aux.timeRef ~= "G")), ","));
    end
end
aux.angleRange_rad = [-1 1] * deg2rad(5758);
aux.provenance = NonGravForces.provenance("io.readIMU1B", 'Files', aux.files, 'Options', options);
end
