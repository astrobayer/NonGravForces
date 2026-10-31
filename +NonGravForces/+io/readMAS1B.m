function [epoch, tank, aux] = readMAS1B(files, options)
%READMAS1B MAS1B reader: cold-gas TANK masses (kg) with GPS epochs. NOT the spacecraft mass.
%
%   [epoch, tank, aux] = NonGravForces.io.readMAS1B(file)
%   [epoch, tank, aux] = NonGravForces.io.readMAS1B(files, 'TimeRange', tr)
%
%   Reads GRACE-FO MAS1B files through gracefoReading. MAS1B REPORTS THE
%   PROPELLANT (cold gas) MASS: about 15 kg per tank, 30 kg in total, hourly.
%   It is not, and must never be used as, the spacecraft mass (601.214 kg for
%   GRACE-C, 601.209 kg for GRACE-D, NonGravForces.const): the predecessor
%   toolbox made that mistake and produced accelerations 20-40 times too large
%   (HANDOVER C.1). This reader therefore refuses any total above 100 kg, and
%   NonGravForces.io.massModel is the only unit that turns these numbers into
%   a spacecraft mass.
%
%   MAS1B is a prod_flag product: each record carries only the columns its
%   flag marks present; absent values are NaN in the outputs and the flags are
%   in aux.prodFlagBits. On the real 2020 files only gas_mass_tnk1 and
%   gas_mass_tnk2 are present.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One or more MAS1B files of ONE satellite; .gz accepted.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type              Units  Description
%   ----        ----  ----              -----  -----------
%   TimeRange   1x2   double|datetime   -      Forwarded to readL1b. Default [].
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   epoch    1x1   struct  -      referenceSystem epoch struct, K entries, "GPS".
%   tank     1x1   struct  kg     .total_kg (1xK: tank1 + tank2 where both are
%                                 present, else mass_tnk), .tank_kg (2xK:
%                                 gas_mass_tnk1, gas_mass_tnk2), .thruster_kg
%                                 (2xK: gas_mass_thr1, _thr2), .tankTotal_kg
%                                 (1xK mass_tnk), .thrusterTotal_kg (1xK
%                                 mass_thr), .sigma_kg (1xK mass_tnk_err),
%                                 .present (struct of logicals per column: any
%                                 record has it). NaN where the flag says absent.
%   aux      1x1   struct  -      .prodFlagBits, .qualflgBits, .product,
%                                 .satellite, .files, .numRecords, .sampling_s,
%                                 .obsNames, .units, .provenance.
%
% ACCURACY
%   Exact (file values). No interpolation between the hourly samples: that is
%   massModel's job, and it documents its own rule.
%
% ERRORS
%   NonGravForces:readMAS1B:dataProduct     a file is not a MAS1B product
%   NonGravForces:readMAS1B:dataSatellite   the files mix satellites
%   NonGravForces:readMAS1B:dataEmpty       no record survives TimeRange
%   NonGravForces:readMAS1B:dataColumn      gas_mass_tnk1 or gas_mass_tnk2 missing
%                                           from the header
%   NonGravForces:readMAS1B:dataMissing     no record has a usable total (both
%                                           tank columns and mass_tnk absent)
%   NonGravForces:readMAS1B:dataRange       a total is >= 100 kg or <= 0 kg
%   NonGravForces:readMAS1B:dataTimeRef     time_ref is not 'G' for every record
%   gracefoReading:*                        propagated unwrapped
%
% REFERENCES
%   JPL D-56935, MAS1B format; HANDOVER C.1 (v1 defect: MAS1B as mass).
%
% EXAMPLE
%   f = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'samples', 'MAS1B_2020-01-01_C_04.txt');
%   [ep, tank] = NonGravForces.io.readMAS1B(f); tank.total_kg(1)   % ~30.5
%
% See also NonGravForces.io.massModel, NonGravForces.const.
%
% KEYWORDS MAS1B tank mass cold gas propellant
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
arguments
    files (1,:) string {mustBeNonempty}
    options.TimeRange = []
end
unit = "readMAS1B";
loadOpt = struct('TimeRange', options.TimeRange, 'QualityBitMask', []);
[epoch, obs, raw, aux] = l1bLoad(unit, files, "MAS1B", loadOpt);
if isfield(raw, 'time_ref')
    tr = unique(string(raw.time_ref));
    if ~(isscalar(tr) && tr == "G")
        NonGravForces.fail(unit, "dataTimeRef", "time_ref is %s, expected 'G' (GPS) -> a non-GPS time reference needs a new code path", strjoin(tr, ","));
    end
end
cols = ["gas_mass_tnk1","gas_mass_tnk2","gas_mass_thr1","gas_mass_thr2","mass_tnk","mass_thr","mass_tnk_err"];
idx = l1bColumns(unit, raw, cols(1:2));
idxOpt = l1bColumns(unit, raw, cols(3:end), false);
idx = [idx idxOpt];
K = size(obs, 1);
val = nan(numel(cols), K);
val(idx > 0, :) = obs(:, idx(idx > 0))';
tank.tank_kg          = val(1:2, :);
tank.thruster_kg      = val(3:4, :);
tank.tankTotal_kg     = val(5, :);
tank.thrusterTotal_kg = val(6, :);
tank.sigma_kg         = val(7, :);
bothTanks = all(isfinite(tank.tank_kg), 1);
tank.total_kg = nan(1, K);
tank.total_kg(bothTanks) = sum(tank.tank_kg(:, bothTanks), 1);
useAgg = ~bothTanks & isfinite(tank.tankTotal_kg);
tank.total_kg(useAgg) = tank.tankTotal_kg(useAgg);
tank.present = struct();
for k = 1:numel(cols)
    tank.present.(cols(k)) = any(isfinite(val(k, :)));
end
if ~any(isfinite(tank.total_kg))
    NonGravForces.fail(unit, "dataMissing", "no record carries gas_mass_tnk1+2 or mass_tnk (prod_flag) -> this file has no usable tank total; check aux.prodFlagBits");
end
finiteTotal = tank.total_kg(isfinite(tank.total_kg));
if any(finiteTotal >= 100) || any(finiteTotal <= 0)
    NonGravForces.fail(unit, "dataRange", "tank total %.3f kg is outside (0, 100) kg: MAS1B is propellant mass, about 30 kg -> this is not a GRACE-FO MAS1B file, or the columns are not gas masses", ...
        finiteTotal(find(finiteTotal >= 100 | finiteTotal <= 0, 1)));
end
if isfield(raw, 'prodFlagBits'), aux.prodFlagBits = logical(raw.prodFlagBits); else, aux.prodFlagBits = false(K, 0); end
aux.provenance = NonGravForces.provenance("io.readMAS1B", 'Files', aux.files, 'Options', options);
end
