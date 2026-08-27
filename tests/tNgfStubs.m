classdef tNgfStubs < NgfTestCase
%TNGFSTUBS Every stub in doc/functionlist.csv raises NonGravForces:<name>:notImplemented.
%
%   This is the mechanism behind "fail loudly": a stub that resolved to
%   something else, or a catalogue row without a file, is caught here.
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    properties (TestParameter)
        stub = tNgfStubs.stubList()
    end
    methods (Static)
        function s = stubList()
            here = fileparts(mfilename('fullpath'));
            fl = readtable(fullfile(fileparts(here), 'doc', 'functionlist.csv'), 'TextType', 'string');
            rows = fl(fl.status == "stub", :);
            s = cellfun(@(p, n) char(p + "." + n), num2cell(rows.package), num2cell(rows.name), 'UniformOutput', false);
        end
    end
    methods (Test, TestTags = {'contract'})
        function stubRaises(tc, stub)
            parts = split(string(stub), ".");
            fcn = str2func("NonGravForces." + parts(1) + "." + parts(2));
            tc.verifyErrorId(fcn, "NonGravForces:" + parts(2) + ":notImplemented");
        end
    end
end
