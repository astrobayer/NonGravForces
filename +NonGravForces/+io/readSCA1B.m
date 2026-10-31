function [epoch, q_body2gcrs, aux] = readSCA1B(files, options)
%READSCA1B SCA1B reader: ICRF->SRF quaternions conjugated ONCE to BODY->GCRS, 4xK scalar first.
%
%   [epoch, q, aux] = NonGravForces.io.readSCA1B(file)
%   [epoch, q, aux] = NonGravForces.io.readSCA1B(files, 'TimeRange', tr, 'QualityBitMask', bits)
%
%   Reads GRACE-FO star-camera attitude (SCA1B) through gracefoReading and
%   returns the attitude in the form the rest of the toolbox stores
%   (HANDOVER D.4): unit quaternions BODY->GCRS, scalar first, one column per
%   epoch. The product stores the rotation ICRF->SRF (quatangle = scalar,
%   quati/j/kcoeff = vector part); the conjugation happens HERE, exactly once,
%   through referenceSystem.invertQuaternion, so no later unit conjugates
%   again. GCRS and ICRF share their orientation by definition (the GCRS is
%   the ICRS-oriented geocentric system), so no frame bias is applied.
%
%   The SRF is the toolbox body frame (geom.graceFO panels are given in it).
%   The sense is proven on real data every run: the body +Z axis rotated into
%   the GCRS lies within 1 degree of nadir on the 2020-01-01 sample, and 96
%   degrees off when the conjugation is omitted (tests/tNgfIo.m).
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One or more SCA1B files of ONE satellite,
%                                 chronological; .gz accepted.
%
% NAME-VALUE ARGUMENTS
%   Name            Size  Type              Units  Description
%   ----            ----  ----              -----  -----------
%   TimeRange       1x2   double|datetime   -      Forwarded to readL1b: decimal
%                                                  MJD or GPS datetime. Default [].
%   QualityBitMask  1xQ   double            -      qualflg bits removing a record;
%                                                  forwarded. Default [].
%
% OUTPUT ARGUMENTS
%   Name          Size  Type    Units  Description
%   ----          ----  ----    -----  -----------
%   epoch         1x1   struct  -      referenceSystem epoch struct, K entries, "GPS".
%   q_body2gcrs   4xK   double  -      Unit quaternions BODY->GCRS, [q0;q1;q2;q3],
%                                      canonical sign q0 >= 0 (as
%                                      referenceSystem.invertQuaternion returns).
%   aux           1x1   struct  -      .frame "BODY->GCRS", .convention
%                                      (referenceSystem.attitudeConvention()),
%                                      .qualRss (1xK, qual_rss column),
%                                      .scaId (1xK), .normDeviationMax (max
%                                      |norm-1| of the file quaternions),
%                                      .qualflgBits, .product, .satellite,
%                                      .files, .numRecords, .sampling_s,
%                                      .obsNames, .units, .provenance.
%
% ACCURACY
%   The conjugate of a unit quaternion is exact; invertQuaternion normalises
%   first, so the returned norm is 1 to 2.2e-16 while the file's own norm
%   deviates by up to 1e-9 (ASCII rounding; 9e-11 measured on the sample). A
%   file quaternion whose norm deviates by more than 1e-6 is refused, not
%   repaired (gracefoReading.checkQuaternionNorm tolerance). The sign
%   canonicalisation can flip q -> -q between neighbouring epochs where q0
%   crosses zero; a rotation matrix is unaffected and
%   referenceSystem.interpolateQuaternions restores continuity itself.
%
% ERRORS
%   NonGravForces:readSCA1B:dataProduct      a file is not an SCA1B product
%   NonGravForces:readSCA1B:dataSatellite    the files mix GRACE-C and GRACE-D
%   NonGravForces:readSCA1B:dataEmpty        no record survives the filters
%   NonGravForces:readSCA1B:dataColumn       a quaternion column is missing
%   NonGravForces:readSCA1B:dataQuaternion   a file quaternion is not of unit norm
%   gracefoReading:*, referenceSystem:*      propagated unwrapped
%
% REFERENCES
%   JPL D-56935, SCA1B format ("rotation from the inertial frame to the SRF").
%   referenceSystem.attitudeConvention: scalar-first Hamilton, passive.
%
% EXAMPLE
%   f = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'samples', 'SCA1B_2020-01-01_C_04.txt');
%   [ep, q, aux] = NonGravForces.io.readSCA1B(f);
%   R = referenceSystem.rotationFromQuaternion(q);   % 3x3xK BODY->GCRS
%
% See also NonGravForces.io.readGNV1B, NonGravForces.geom.AttitudeProvider, referenceSystem.invertQuaternion.
%
% KEYWORDS SCA1B attitude quaternion star camera conjugate
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
arguments
    files (1,:) string {mustBeNonempty}
    options.TimeRange = []
    options.QualityBitMask (1,:) double = []
end
unit = "readSCA1B";
[epoch, obs, raw, aux] = l1bLoad(unit, files, "SCA1B", options);
iQ = l1bColumns(unit, raw, ["quatangle","quaticoeff","quatjcoeff","quatkcoeff"]);
iR = l1bColumns(unit, raw, "qual_rss", false);
[ok, stats] = gracefoReading.checkQuaternionNorm(obs, raw);
if ~all(ok)
    NonGravForces.fail(unit, "dataQuaternion", "%d of %d quaternions deviate from unit norm by up to %.2e (tolerance %.1e) -> the file is corrupt or not SCA1B; re-download it", ...
        stats.numViolations, stats.n, stats.maxDeviation, stats.tolerance);
end
q_icrf2srf = obs(:, iQ)';
[q_body2gcrs, convention] = referenceSystem.invertQuaternion(q_icrf2srf);
aux.frame = "BODY->GCRS";
aux.convention = convention;
if iR > 0, aux.qualRss = obs(:, iR)'; else, aux.qualRss = []; end
if isfield(raw, 'sca_id'), aux.scaId = double(raw.sca_id(:)'); else, aux.scaId = []; end
aux.normDeviationMax = stats.maxDeviation;
aux.provenance = NonGravForces.provenance("io.readSCA1B", 'Files', aux.files, 'Options', options);
end
