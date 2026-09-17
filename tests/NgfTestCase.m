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
    methods (TestMethodSetup)
        function hushExtrapolation(tc)
            %HUSHEXTRAPOLATION Silence O-42's extrapolation warning for every test by default.
            %   Most of this suite evaluates adsorption models at GRACE-FO densities on purpose,
            %   three orders below where those isotherms were fitted, so after the warn ruling of
            %   2026-09-17 they emitted about twenty-five copies of the very warning one test
            %   asserts. A run log full of an EXPECTED warning is exactly where an unexpected one
            %   hides, so it is off by default and `loud` turns it on where it is the subject.
            st = warning('off', 'NonGravForces:accommodation:outsideFit');
            tc.addTeardown(@() warning(st));
        end
    end

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
        function varargout = loud(~, fcn)
            %LOUD Evaluate fcn with the accommodation extrapolation warning ENABLED for this call.
            %   The suppression above is per test method, so the one test that asserts the warning
            %   turns it back on here, which also makes the intent visible at the assertion.
            %   RESTORED WHEN THIS FUNCTION RETURNS, not at teardown: with addTeardown the warning
            %   stayed on for the REST of the calling test, and the three later calls in that same
            %   method went on printing it -- the flood this helper exists to prevent, moved rather
            %   than removed (measured 2026-09-17: 25 warnings became 3, not 0).
            st = warning('on', 'NonGravForces:accommodation:outsideFit');
            restore = onCleanup(@() warning(st));                                 %#ok<NASGU>
            [varargout{1:nargout}] = fcn();
        end
        function rec = verifyBoundAndRecord(tc, name, measured, bound, units, note)
            %VERIFYBOUNDANDRECORD measured <= bound, and the number is recorded.
            arguments
                tc; name (1,1) string; measured (1,1) double; bound (1,1) double
                units (1,1) string = ""; note (1,1) string = ""
            end
            rec = struct('name', name, 'measured', measured, 'bound', bound, 'units', units, 'note', note, 'class', string(class(tc)));
            NgfTestCase.recordStore('add', rec);
            % A NaN measurement must never pass. MATLAB's max, min and sum omit NaN by default, so a
            % unit that returns one can reach here as a clean number; on 2026-09-08 exactly that hid a
            % NaN from geom.selfShadow behind a bound that read 0 (rule R4).
            if ~all(isfinite(measured))
                tc.verifyFail(sprintf('%s: the measurement is %s, not a number (%s). A NaN here usually means the unit under test returned one and an omitnan reduction hid it.', name, mat2str(measured), note));
                return
            end
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
