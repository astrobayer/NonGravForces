classdef tNgfTutorials < NgfTestCase
%TNGFTUTORIALS Every first-time tutorial in tutorials/ runs to completion on the committed data.
%
%   WHY THIS EXISTS. The tutorials are documentation in three places -- the scripts, their published
%   Help-browser pages (tools/publishTutorials.m) and the tutorials chapter of the user guide PDF,
%   built from a RECORDED run (tools/runTutorials.m, tools/makeTutorialsTex.py). A recording proves
%   that the code ran once; it does not prove that it still runs. This class runs every tutorial
%   whose header says "% Requirements: none" under the gate, so an API that moves under a tutorial
%   reddens the gate the day it moves. A tutorial that needs the data pool says
%   "% Requirements: data pool" and is recorded on the CI runner only.
%
%   WHAT IS ASSERTED is completion and a non-empty printout, as tNgfDemos does for demos: the
%   numbers a tutorial prints are bounded by the unit tests of the units it calls, and asserting
%   them twice would put a second authority in a file whose purpose is explanation.
%
%   The list is discovered from the folder, so a new tutorial is covered the moment it is added.
%
% Created 2026-11-12 by Claude (Anthropic) for Matthias Weigelt, after lgdToolbox's
% tDocumentation/tutorialRuns (R-055).
    properties (TestParameter)
        tutorial = tNgfTutorials.tutorialList("none")
    end
    methods (Static)
        function s = tutorialList(requirement)
            % The tutorials whose header line "% Requirements: <x>" names this requirement.
            here = fileparts(mfilename('fullpath'));
            d = dir(fullfile(fileparts(here), 'tutorials', 'ngfTutorial*.m'));
            s = {};
            for k = 1:numel(d)
                req = tNgfTutorials.requirementOf(fullfile(d(k).folder, d(k).name));
                if req == requirement
                    s{end + 1} = d(k).name(1:end-2); %#ok<AGROW>
                end
            end
        end
        function req = requirementOf(file)
            tok = regexp(fileread(file), '(?m)^%\s*Requirements:\s*([a-z ]+?)\s*$', 'tokens', 'once');
            if isempty(tok), req = ""; else, req = string(tok{1}); end
        end
        function out = runTutorial(file)
            % In a function workspace of its own, with figures invisible and closed afterwards.
            old = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            cleanup = onCleanup(@() tNgfTutorials.restore(old));
            out = evalc('tNgfTutorials.execute(file)');
        end
        function execute(ngfTut__file)
            run(ngfTut__file);
        end
        function restore(old)
            close all force
            set(groot, 'DefaultFigureVisible', old);
        end
    end
    methods (Test, TestTags = {'contract'})
        function tutorialRuns(tc, tutorial)
            % Every offline tutorial runs on the current code. What the PDF shows of it
            % (doc/src/tutorials/*/run.json) is checked against the script by
            % tools/makeTutorialsTex.py --check (auditPackage.py check 15).
            root = fileparts(fileparts(mfilename('fullpath')));
            file = fullfile(root, 'tutorials', tutorial + ".m");
            tc.assertTrue(isfile(file), 'the list is built from the folder; a missing file means it changed mid-run');
            out = tNgfTutorials.runTutorial(file);
            tc.verifyGreaterThan(strlength(out), 0, 'a tutorial that prints nothing shows nothing');
        end
        function tutorialRequirementsDeclared(tc)
            % Every tutorial declares what it needs, one of three words; the gate runs the "none" ones
            % and there is at least one of them.
            root = fileparts(fileparts(mfilename('fullpath')));
            d = dir(fullfile(root, 'tutorials', 'ngfTutorial*.m'));
            tc.assertNotEmpty(d, 'tutorials/ holds no ngfTutorial*.m');
            for k = 1:numel(d)
                req = tNgfTutorials.requirementOf(fullfile(d(k).folder, d(k).name));
                tc.verifyTrue(ismember(req, ["none" "data pool" "network"]), ...
                    sprintf('%s: its "%% Requirements:" line must say none, data pool or network', d(k).name));
            end
            tc.verifyNotEmpty(tNgfTutorials.tutorialList("none"));
        end
    end
end
