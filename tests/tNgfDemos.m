classdef tNgfDemos < NgfTestCase
%TNGFDEMOS Every demo in demos/ runs to completion under the gate.
%
%   WHY THIS EXISTS. A demo nobody executes rots exactly as a hand-written
%   documentation page does, and this project watched that happen: until
%   2026-09-09 doc/html/Sensors.html credited sensor.Electrostatic with a filter
%   the class deliberately does not build, because nothing regenerated it. The
%   answer there was to generate the page; the answer here is to RUN the demo.
%
%   WHAT IS ASSERTED is completion, not output. A demo's job is to show a
%   working sequence of public calls; the numbers it prints are already bounded
%   by the unit tests of the units it calls, and asserting them twice would put
%   a second authority in a file whose purpose is illustration. What this
%   catches is the failure a demo actually suffers: an API that moved.
%
%   The list is discovered from the folder, so a new demo is covered the moment
%   it is added and nobody has to remember to register it.
%
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
    properties (TestParameter)
        demo = tNgfDemos.demoList()
    end
    methods (Static)
        function s = demoList()
            here = fileparts(mfilename('fullpath'));
            d = dir(fullfile(fileparts(here), 'demos', 'demo*.m'));
            s = cellfun(@(n) n(1:end-2), {d.name}, 'UniformOutput', false);
        end
    end
    methods (Test, TestTags = {'contract'})
        function demoRuns(tc, demo)
            % Demos print by design -- they are top-level scripts, not library code -- and the
            % output is swallowed so the gate's own report stays readable. evalc also means a
            % demo that starts printing megabytes is a memory problem here rather than a wall of
            % text in the run log.
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.assertTrue(isfile(fullfile(root, 'demos', demo + ".m")), ...
                'the demo list is built from the folder, so a missing file here means the folder changed mid-run');
            old = pwd;
            cleanup = onCleanup(@() cd(old));
            cd(fullfile(root, 'demos'));
            out = evalc(demo);
            tc.verifyGreaterThan(strlength(out), 0, ...
                'a demo that prints nothing shows nothing; give it at least one line of result');
            tc.verifyEmpty(regexp(out, 'NaN|Inf', 'once'), ...
                'a demo must not present NaN or Inf as a result');
        end
    end
end
