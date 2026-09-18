classdef tNgfTools < NgfTestCase
%TNGFTOOLS The development instruments in tools/ run to completion under the gate.
%
%   WHY THIS EXISTS. tools/benchmark.m was written on 2026-09-09 to regenerate the
%   per-epoch cost of HANDOVER O-29, because "a number that cannot be regenerated
%   is a number nobody can check". It was then not run for nine days, and when it
%   was finally run it did not work: its native row asked app.DynamicsSource for
%   attitude from -60 s while the fixture's attitude starts at 0, so prepare()
%   refused and the function errored for EVERY value of 'Epochs', its own default
%   included. Written, committed, cited in Model.m's help, described in HANDOVER,
%   executed never. In those nine days its stale 88 ms was copied into four more
%   files and one shipped PDF chapter.
%
%   An instrument nobody runs is the same disease as a number nobody can check,
%   and it has the same cure. tests/tNgfDemos.m already says it for demos: "The
%   answer there was to generate the page; the answer here is to RUN the demo."
%   This is that answer for tools/.
%
%   WHAT IS ASSERTED is completion and the SHAPE the help promises -- never a
%   duration. The whole argument for keeping benchmark.m out of the gate stands:
%   it measures the machine as much as the code, and a gate that reddens because
%   a laptop was busy teaches people to ignore gates. Nothing here looks at a
%   millisecond. What this catches is the failure the instrument actually
%   suffered: a collaborator's signature moved and nobody noticed.
%
%   'Write' IS FALSE, and that is load-bearing rather than tidy: with it true the
%   gate would overwrite tools/benchmark_reference.csv, which auditPackage.py
%   check 13 treats as the authority for every per-epoch figure in the
%   repository. A test run must not restate the fact it is supposed to protect.
%
% Created 2026-09-18 by Claude Opus 5 (Anthropic).
    methods (Test, TestTags = {'contract'})
        function benchmarkRunsAndKeepsItsShape(tc)
            % The smallest run that still exercises every path: one batch size, one timed
            % repetition. The cost is dominated by reading the arc and building the two Models,
            % not by the timing loop, so asking for more would buy no coverage.
            T = benchmark('Epochs', 1, 'Repeats', 1, 'Write', false);

            tc.verifyClass(T, 'table');
            tc.verifyEqual(string(T.Properties.VariableNames), ...
                ["path" "epochs" "seconds" "perEpoch_ms" "note"], ...
                'the help documents these five columns; a caller selects by name');
            tc.verifyGreaterThanOrEqual(height(T), 5, ...
                'four batched paths plus at least the compiled single-epoch row');

            % Every batched path must have produced a finite measurement. The two POINTWISE rows
            % may legitimately be NaN -- app.DynamicsSource needs gjToolbox and the rhs MEX, and
            % its absence is REPORTED as a row rather than degraded into silence (R4) -- so the
            % assertion is that a NaN row SAYS why, not that no NaN exists.
            for k = 1:height(T)
                if ~isfinite(T.perEpoch_ms(k))
                    tc.verifySubstring(char(T.note(k)), 'NOT MEASURED', ...
                        sprintf('row "%s" has no measurement and no reason -> a path this host cannot measure must say so in its note', T.path(k)));
                else
                    tc.verifyGreaterThan(T.perEpoch_ms(k), 0, ...
                        sprintf('row "%s" measured a non-positive per-epoch cost', T.path(k)));
                end
            end

            % The batch size travels with the measurement (R2): a per-epoch cost without its K is
            % the exact shape of claim that made O-29 survive ten days of quotation.
            tc.verifyTrue(all(T.epochs >= 1), 'every row records the batch size it was timed at');
        end

        function benchmarkRefusesAnArcItDoesNotHave(tc)
            % The committed GNV1B excerpt is an excerpt, so a large enough batch has no states to
            % time. The count is deliberately not written here -- benchmark.m reads it from the
            % file and puts it in the message; 100000 is past any excerpt this fixture will hold.
            % Until 2026-09-18 this reached ep.day(1:n) and raised MATLAB's "Index exceeds the
            % number of array elements" -- a message naming neither the fixture nor the limit, from
            % no unit the caller can look up. Paired with the call that must SUCCEED, because a
            % refusal test whose sibling does not run proves only that everything refuses.
            tc.verifyError(@() benchmark('Epochs', 100000, 'Write', false), ...
                "NonGravForces:benchmark:inputRange");
            tc.verifyError(@() benchmark('Epochs', 0, 'Write', false), ...
                "NonGravForces:benchmark:inputRange");
        end
    end
end
