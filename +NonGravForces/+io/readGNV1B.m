function [epoch, r, v, aux] = readGNV1B(files, options)
%READGNV1B GNV1B reader: Earth-fixed state rotated to the GCRS through referenceSystem, 3xK.
%
%   [epoch, r, v, aux] = NonGravForces.io.readGNV1B(file, 'Eop', eop)
%   [epoch, r, v, aux] = NonGravForces.io.readGNV1B(files, 'Eop', eopFile, 'TimeRange', tr)
%   [epoch, r, v, aux] = NonGravForces.io.readGNV1B(file, 'Frame', "ITRS")
%
%   Reads GRACE-FO GNV1B orbit files through gracefoReading. GNV1B IS
%   EARTH-FIXED (coord_ref = 'E', ITRS), not inertial - the pitfall recorded
%   in HANDOVER C.1. With the default 'Frame' "GCRS" the state is rotated by
%   referenceSystem.rotationItrsToGcrs and applied with
%   referenceSystem.applyRotation, which carries the frame-rate term
%   (Rdot r, about 470 m/s) into the velocity; this function never forms
%   omega x r by hand. The rotation needs Earth orientation parameters, and
%   there is no default series: pass 'Eop' (rule R5). With 'Frame' "ITRS" the
%   state is returned as printed and no EOP is needed.
%
%   EOP SERIES NOTE (re-measured 2026-09-07 against referenceSystem 2.1.216,
%   HANDOVER O-10 closed): a series carrying rows before 1972-01-01 is now accepted.
%   referenceSystem trims those rows and warns
%   (referenceSystem:interpolateEop:*), because UTC before that date is a rate-based
%   fractional offset it does not model; until 2.1.216 it refused outright. The
%   shipped eopc04_20.1962-now therefore works read whole, with 3652 of 23562 rows
%   ignored. finals2000A.all and the trimmed C04 agree to 9 mm on the
%   2020-01-01 arc.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One or more GNV1B files of ONE satellite,
%                                 chronological; .gz accepted.
%
% NAME-VALUE ARGUMENTS
%   Name            Size  Type              Units  Description
%   ----            ----  ----              -----  -----------
%   Eop             1x1   struct|string     -      EOP series from
%                                                  referenceSystem.readEop, or
%                                                  the path of a series file
%                                                  (forwarded as 'EopFile').
%                                                  REQUIRED for Frame "GCRS";
%                                                  ignored for "ITRS". Default
%                                                  [] (none).
%   Frame           1x1   string            -      "GCRS" (default) or "ITRS".
%   TimeRange       1x2   double|datetime   -      Forwarded to readL1b. Default [].
%   QualityBitMask  1xQ   double            -      Forwarded to readL1b. Default [].
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   epoch    1x1   struct  -      referenceSystem epoch struct, K entries, "GPS".
%   r        3xK   double  m      Position in aux.frame.
%   v        3xK   double  m/s    Velocity in aux.frame (GCRS: includes Rdot r).
%   aux      1x1   struct  -      .frame ("GCRS"|"ITRS"), .sigma_r_m (3xK,
%                                 xpos_err.., ITRS as printed), .sigma_v_m_s
%                                 (3xK), .rotation (info struct of
%                                 rotationItrsToGcrs; [] for ITRS),
%                                 .eopSeries (eop.meta.series or the file
%                                 path; "" for ITRS), .qualflgBits, .product,
%                                 .satellite, .files, .numRecords,
%                                 .sampling_s, .obsNames, .units, .provenance.
%
% ACCURACY
%   The rotation is referenceSystem's (IAU 2006/2000A, CIO based, tidal EOP
%   corrections on by default; Rdot by central differences, step 200 s, see
%   its ACCURACY block). Round trip GCRS->ITRS->GCRS reproduces the file to
%   5e-9 m on the sample (double rounding at 7e6 m); the velocity agrees with
%   a five-point finite difference of the GCRS positions to 4.4e-5 m/s
%   (dominated by the 1 Hz sampling and the GNV1B position noise), against
%   470 m/s if the Rdot term were dropped. The formal errors are NOT rotated:
%   rotating a per-axis sigma without its covariance is meaningless, so they
%   stay in the ITRS and are labelled so.
%
% ERRORS
%   NonGravForces:readGNV1B:inputEop        Frame "GCRS" without a usable Eop
%   NonGravForces:readGNV1B:dataProduct     a file is not a GNV1B product
%   NonGravForces:readGNV1B:dataSatellite   the files mix satellites
%   NonGravForces:readGNV1B:dataEmpty       no record survives the filters
%   NonGravForces:readGNV1B:dataColumn      a position/velocity column is missing
%   NonGravForces:readGNV1B:dataFrame       coord_ref is not 'E' for every record
%   referenceSystem:*, gracefoReading:*     propagated unwrapped (a pre-1972 EOP
%                                           series is trimmed with a warning,
%                                           not refused; see the note above)
%
% REFERENCES
%   JPL D-56935, GNV1B format; IERS Conventions (2010) Ch. 5 via referenceSystem.
%
% EXAMPLE
%   [~, loc] = referenceSystem.dataFolder();                  % loc.EOP = <pool>/EOP
%   eop = referenceSystem.readEop(fullfile(loc.EOP, 'finals2000A.all'));
%   f = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'samples', 'GNV1B_2020-01-01_C_04.txt');
%   [ep, r, v, aux] = NonGravForces.io.readGNV1B(f, 'Eop', eop);   % GCRS
%
% See also NonGravForces.io.readSCA1B, referenceSystem.rotationItrsToGcrs, referenceSystem.applyRotation.
%
% KEYWORDS GNV1B orbit ITRS GCRS rotation EOP
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
arguments
    files (1,:) string {mustBeNonempty}
    options.Eop = []
    options.Frame (1,1) string {mustBeMember(options.Frame, ["GCRS","ITRS"])} = "GCRS"
    options.TimeRange = []
    options.QualityBitMask (1,:) double = []
end
unit = "readGNV1B";
% Eop is validated BEFORE reading so a wrong call fails in milliseconds.
eopArg = {};
eopSeries = "";
if options.Frame == "GCRS"
    e = options.Eop;
    if isstruct(e) && all(isfield(e, ["mjd","xp","yp","dut1","lod","dX","dY","meta"]))
        eopArg = {'Eop', e};
        eopSeries = string(e.meta.series);
    elseif (isstring(e) || ischar(e)) && isscalar(string(e)) && isfile(string(e))
        eopArg = {'EopFile', char(string(e))};
        eopSeries = string(e);
    else
        NonGravForces.fail(unit, "inputEop", ...
            "Frame GCRS needs Earth orientation parameters -> pass 'Eop', referenceSystem.readEop(file) or the file path (finals2000A.all, or a C04 series from 1972 on), or request Frame ITRS");
    end
end
[epoch, obs, raw, aux] = l1bLoad(unit, files, "GNV1B", options);
iR = l1bColumns(unit, raw, ["xpos","ypos","zpos"]);
iV = l1bColumns(unit, raw, ["xvel","yvel","zvel"]);
iSr = l1bColumns(unit, raw, ["xpos_err","ypos_err","zpos_err"], false);
iSv = l1bColumns(unit, raw, ["xvel_err","yvel_err","zvel_err"], false);
ref = unique(string(raw.coord_ref));
if ~(isscalar(ref) && ref == "E")
    NonGravForces.fail(unit, "dataFrame", "coord_ref is %s, this reader knows only 'E' (Earth-fixed) -> check the product; an inertial GNV1B variant needs a new code path", ...
        strjoin(ref, ","));
end
rI = obs(:, iR)';
vI = obs(:, iV)';
if options.Frame == "GCRS"
    [R, Rdot, ~, info] = referenceSystem.rotationItrsToGcrs(epoch, eopArg{:});
    [r, v] = referenceSystem.applyRotation(R, Rdot, rI, vI);
    aux.rotation = info;
else
    r = rI; v = vI;
    aux.rotation = [];
end
aux.frame = options.Frame;
aux.eopSeries = eopSeries;
if all(iSr > 0), aux.sigma_r_m = obs(:, iSr)'; else, aux.sigma_r_m = []; end
if all(iSv > 0), aux.sigma_v_m_s = obs(:, iSv)'; else, aux.sigma_v_m_s = []; end
optRec = options; optRec.Eop = eopSeries;   % the hash must not swallow a 23k-row table
aux.provenance = NonGravForces.provenance("io.readGNV1B", 'Files', aux.files, 'Options', optRec);
end
