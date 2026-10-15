function T = eventTable(varargin)
%EVENTTABLE [DEPRECATED, forwards to timeSeriesPreparation.qc.eventTable] Canonical event table (schema authority) for quality-control results.
%
%   T = NonGravForces.qc.eventTable()                       % empty, schema only
%   T = NonGravForces.qc.eventTable(tStart, tEnd, class, axis, amplitude, sigma, source, confidence)
%
%   DEPRECATED (N-01, 2026-10-14). The qc stage moved to timeSeriesPreparation
%   (CONCEPT 10.1 N-01, X-05), which ported this schema UNCHANGED. This unit now
%   only forwards: it warns NonGravForces:eventTable:deprecated and calls
%   timeSeriesPreparation.qc.eventTable with the same arguments. It is deleted
%   in the release after the one that carries the forwarding. io.readTHR1B
%   builds its table with the target directly. The target's errors propagate
%   unwrapped (timeSeriesPreparation:eventTable:*).
%
% INPUT ARGUMENTS
%   Name      Size  Type  Units Description
%   ----      ----  ----  ----- -----------
%   varargin  1x0 or 1x8 cell - the eight columns of timeSeriesPreparation.qc.eventTable, or none
%
% NAME-VALUE ARGUMENTS
%   None.
%
% OUTPUT ARGUMENTS
%   T         Nx8   table  -    the target's table, unchanged
%
% ACCURACY
%   None of its own: the target's table, unchanged.
%
% ERRORS
%   NonGravForces:eventTable:depMissing   timeSeriesPreparation is not on the path
%   timeSeriesPreparation:eventTable:*    the target's errors, unwrapped
%
% WARNINGS
%   NonGravForces:eventTable:deprecated   on every call; call timeSeriesPreparation.qc.eventTable
%
% REFERENCES
%   timeSeriesPreparation docs/CONCEPT.md section 10.1 rows N-01, X-05.
%
% EXAMPLE
%   % Call the target directly; this name disappears in the release after next.
%   T = timeSeriesPreparation.qc.eventTable();
%
% See also timeSeriesPreparation.qc.eventTable.
%
% KEYWORDS deprecated forwarder event table schema
% Created 2026-09-07 by Claude Fable 5 (Anthropic); reduced to a forwarder (N-01) 2026-10-14 by Claude (Anthropic).
arguments (Repeating)
    varargin
end
if isempty(which('timeSeriesPreparation.qc.eventTable'))
    NonGravForces.fail("eventTable", "depMissing", "qc.eventTable forwards to timeSeriesPreparation.qc.eventTable and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation");
end
warning("NonGravForces:eventTable:deprecated", "NonGravForces.qc.eventTable is deprecated and forwards to timeSeriesPreparation.qc.eventTable -> call that directly; this name is removed in the release after the one that carries the forwarding");
T = timeSeriesPreparation.qc.eventTable(varargin{:});
end
