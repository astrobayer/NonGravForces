function ex = testCategoryExemptions()
%TESTCATEGORYEXEMPTIONS One reason per exempted (unit, category). Rule T2.
%
%   An exemption claims a test is IMPOSSIBLE, not inconvenient. contract may
%   never be exempted for a unit with an arguments block.
%
% KEYWORDS test exemptions
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
ex = struct();
ex.toolboxVersion.reference     = "the version is the authority; there is nothing outside the project to compare it with";
ex.toolboxVersion.precision     = "no numerical claim";
ex.toolboxVersion.vectorization = "scalar by definition";
ex.toolboxVersion.metamorphic   = "no input relation exists";
ex.toolboxVersion.speed         = "no numerical work";
ex.dependencyReport.reference   = "measures the host; the host is the reference";
ex.dependencyReport.precision   = "no numerical claim";
ex.dependencyReport.vectorization = "scalar by definition";
ex.dependencyReport.metamorphic = "no input relation exists";
ex.dependencyReport.speed       = "no numerical work";
ex.fail.reference               = "MATLAB MException semantics are the reference; asserted under contract";
ex.fail.precision               = "no numerical claim";
ex.fail.vectorization           = "scalar by definition";
ex.fail.speed                   = "no numerical work";
ex.provenance.reference         = "no external authority for a metadata record";
ex.provenance.precision         = "no numerical claim";
ex.provenance.vectorization     = "scalar by definition";
ex.provenance.speed             = "no numerical work";
ex.Options.reference            = "preset tables are the authority (HANDOVER D.4)";
ex.Options.precision            = "no numerical claim";
ex.Options.vectorization        = "scalar object";
ex.Options.speed                = "no numerical work";
ex.eventTable.reference         = "schema is the authority";
ex.eventTable.precision         = "no numerical claim";
ex.eventTable.vectorization     = "table constructor; N rows tested under contract";
ex.eventTable.speed             = "no numerical work";
ex.graceFO.vectorization        = "factory returns one object";
ex.graceFO.speed                = "no numerical work";
ex.grace.vectorization          = "factory returns one object";
ex.grace.speed                  = "no numerical work";
ex.MacroModel.vectorization     = "constructor of one object";
ex.fitCalibration.vectorization = "one call fits one arc; a campaign of arcs is a different estimator (N-09), and the three axes are asserted independent under metamorphic";
ex.MacroModel.speed             = "no numerical work";
ex.const.vectorization          = "scalar struct";
ex.const.speed                  = "no numerical work";
ex.const.metamorphic            = "constants have no input";
ex.buildMex.reference           = "the compiler and the vendored C are the reference; the density table test is the reference of what it builds";
ex.buildMex.precision           = "no numerical claim of its own";
ex.buildMex.vectorization       = "one build";
ex.buildMex.metamorphic         = "no input relation exists";
ex.buildMex.speed               = "a compile is not timed";
ex.boxWing.vectorization        = "factory returns one object";
ex.boxWing.speed                = "no numerical work";
ex.accommodation.speed          = "two closed forms per epoch; there is no algorithm whose cost a ratio could guard";
ex.assemble.reference           = "two sums and a cross product have no authority outside arithmetic; closed forms are asserted under precision";
end
