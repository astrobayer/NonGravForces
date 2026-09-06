classdef NgfTestCase < matlab.unittest.TestCase
%NGFTESTCASE Base class for NonGravForces tests: shared checks and the record store.
%
%   classdef tNgfSomething < NgfTestCase
%
%   Seven test categories are mandatory per public unit, expressed as method
%   tags and counted by auditTestCategories:
%     contract      shapes, rejection, error identifiers, option validation
%     reference     agreement with an authority outside this project
%     precision     a numerical claim with a tolerance justified in a comment
%     vectorization K epochs at once == the same epochs one at a time
%     robustness    NaN, gaps, empty input, flags: defined behaviour, no crash
%     metamorphic   a relation that must hold without a known answer
%     speed         a ratio or a budget, measured in this run on this host
%
%   Every bound that passes leaves its number in the record store
%   (verifyBoundAndRecord), so the run report carries what was measured.
%
% KEYWORDS test base class categories record
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    methods
        function verifyErrorId(tc, fcn, identifier)
            %VERIFYERRORID fcn() must throw exactly IDENTIFIER (and must throw).
            try
                fcn();
                tc.verifyFail(sprintf('expected error %s, nothing was thrown', identifier));
            catch ME
                where = "";
                for k = 1:min(3, numel(ME.stack))
                    where = where + sprintf(' <- %s:%d', ME.stack(k).name, ME.stack(k).line);
                end
                tc.verifyEqual(ME.identifier, char(identifier), ...
                    sprintf('wrong identifier: got %s, message: %s%s', ME.identifier, ME.message, where));
            end
        end
        function rec = verifyBoundAndRecord(tc, name, measured, bound, units, note)
            %VERIFYBOUNDANDRECORD measured <= bound, and the number is recorded.
            arguments
                tc; name (1,1) string; measured (1,1) double; bound (1,1) double
                units (1,1) string = ""; note (1,1) string = ""
            end
            rec = struct('name', name, 'measured', measured, 'bound', bound, 'units', units, 'note', note, 'class', string(class(tc)));
            NgfTestCase.recordStore('add', rec);
            tc.verifyLessThanOrEqual(measured, bound, sprintf('%s: %.3e exceeds %.3e %s (%s)', name, measured, bound, units, note));
        end
        function verifyRelClose(tc, actual, expected, relTol, name)
            %VERIFYRELCLOSE relative closeness with an absolute floor at eps-level.
            arguments
                tc; actual double; expected double; relTol (1,1) double; name (1,1) string = "value"
            end
            scale = max(abs(expected(:)));
            if scale == 0, scale = 1; end
            err = max(abs(actual(:) - expected(:))) / scale;
            tc.verifyBoundAndRecord(name, err, relTol, "rel", "max |a-e| / max |e|");
        end
    end
    methods (Static)
        function out = recordStore(action, rec)
            persistent store
            if isempty(store), store = struct('name',{},'measured',{},'bound',{},'units',{},'note',{},'class',{}); end
            switch action
                case 'add',   store(end+1) = rec; out = [];
                case 'get',   out = store;
                case 'reset', store = struct('name',{},'measured',{},'bound',{},'units',{},'note',{},'class',{}); out = [];
            end
        end
        function T = recordTable()
            s = NgfTestCase.recordStore('get');
            if isempty(s), T = table(); else, T = struct2table(s, 'AsArray', true); end
        end
    end
end
