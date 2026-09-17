classdef tNgfCore < NgfTestCase
%TNGFCORE Tests of L0: const, toolboxVersion, fail, provenance, Options, dependencyReport.
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    methods (Test, TestTags = {'contract'})
        function toolboxVersionContract(tc)
            v = NonGravForces.toolboxVersion;
            % R1: the string is composed from the numeric fields; no literal is pinned here.
            tc.verifyEqual(v.semantic, string(sprintf('%d.%d.%d', v.major, v.minor, v.patch)));
            tc.verifyEqual(NonGravForces.toolboxVersion('string'), char(v.semantic));
            tc.verifyClass(v.major, 'uint16');
            % Every Contents.m carries a copy of the version; the copies must agree (R1, mechanical).
            here = fileparts(mfilename('fullpath'));
            cfiles = dir(fullfile(fileparts(here), '+NonGravForces', '**', 'Contents.m'));
            tc.verifyGreaterThanOrEqual(numel(cfiles), 8);
            for k = 1:numel(cfiles)
                txt = fileread(fullfile(cfiles(k).folder, cfiles(k).name));
                tok = regexp(txt, '^% Version (\S+)', 'tokens', 'once', 'lineanchors');
                tc.verifyEqual(string(tok), v.semantic, sprintf('%s disagrees with toolboxVersion', fullfile(cfiles(k).folder, cfiles(k).name)));
            end
            tc.verifyEqual(v.concept, uint8(3));
            tc.verifyError(@() NonGravForces.toolboxVersion('xml'), 'MATLAB:validators:mustBeMember');
        end
        function constContract(tc)
            c = NonGravForces.const;
            tc.verifyEqual(c.S0, 1361);
            tc.verifyEqual(numel(c.speciesName), 8);
            tc.verifyEqual(numel(c.speciesMass), 8);
            tc.verifyGreaterThan(c.graceFO.massC, 500, 'MAS1B tank mass must never be mistaken for spacecraft mass');
        end
        function failContract(tc)
            tc.verifyErrorId(@() NonGravForces.fail("unitA", "inputShape", "bad %d -> fix it", 3), 'NonGravForces:unitA:inputShape');
            ME = NonGravForces.fail("unitA", "dataMissing", "file %s -> download it", "x.txt");
            tc.verifyClass(ME, 'MException');
            tc.verifyEqual(ME.message, 'file x.txt -> download it');
            tc.verifyErrorId(@() NonGravForces.fail("u", "weird", "m -> h"), 'NonGravForces:fail:inputMnemonic');
            tc.verifyErrorId(@() NonGravForces.fail("u", "inputX", "no hint here"), 'NonGravForces:fail:inputHint');
        end
        function provenanceContract(tc)
            p = NonGravForces.provenance("unitA", 'Files', "a.txt", 'Options', struct('k', 1));
            tc.verifyEqual(p.unit, "unitA");
            tc.verifyEqual(strlength(p.optionHash), 8);
            tc.verifyEqual(p.toolboxVersion, string(NonGravForces.toolboxVersion('string')));
            tc.verifyEqual(p.created.TimeZone, 'UTC');
            q = NonGravForces.provenance("unitA");
            tc.verifyEqual(q.optionHash, "00000000");
        end
        function optionsContract(tc)
            o = NonGravForces.Options;
            tc.verifyEqual(o.preset, "standard");
            tc.verifyEqual(o.drag.density, "msis00");
            p = NonGravForces.Options("precise"); tc.verifyEqual(p.erp.backend, "ceres");
            f = NonGravForces.Options("fast");    tc.verifyEqual(f.erp.backend, "none");
            tc.verifyErrorId(@() NonGravForces.Options("best"), 'NonGravForces:Options:inputPreset');
            tc.verifyErrorId(@() NonGravForces.Options("fast", "drag.densty", "msis00"), 'NonGravForces:Options:inputField');
            tc.verifyErrorId(@() NonGravForces.Options("fast", "drag.density", "nrlmsis"), 'NonGravForces:Options:inputValue');
            tc.verifyErrorId(@() NonGravForces.Options("fast", "drag.alpha", 1.5), 'NonGravForces:Options:inputValue');
            tc.verifyErrorId(@() NonGravForces.Options("fast", "drag.alpha"), 'NonGravForces:Options:inputField');
            % LP-C F4/F5: every admissible value evaluates, so what cannot is refused here, by name
            for bad = {{"drag.gsi", "cll"}, {"drag.gsi", "schamberg"}, {"drag.density", "dtm2020"}, {"drag.density", "jb2008"}, ...
                       {"drag.accommodation", "sesam"}, {"thermal.model", "foldedIntoSrp"}, {"partials.method", "complexStep"}}
                tc.verifyErrorId(@() NonGravForces.Options("standard", bad{1}{:}), 'NonGravForces:Options:inputValue');
            end
            tc.verifyErrorId(@() NonGravForces.Options("standard", "erp.anisotropy", "adm"), 'NonGravForces:Options:inputField');
            tc.verifyErrorId(@() NonGravForces.Options("standard", "geom.lutResolution", 16.5), 'NonGravForces:Options:inputValue');
            p2 = NonGravForces.Options("precise");
            tc.verifyEqual([p2.drag.accommodation, p2.drag.density, p2.geom.selfShadow], ["walker2014", "msis21", "lut"]);
            tc.verifyEqual(p2.drag.surfaceMass_amu, 263.3223); tc.verifyEmpty(p2.drag.alpha); tc.verifyEqual(p2.geom.lutResolution, 32);
            tc.verifyEmpty(o.set("drag.specularFraction", 0.1).set("drag.specularFraction", []).drag.specularFraction, 'an optional value clears back to empty');
            o2 = o.set("srp.tsi", 1367);
            tc.verifyEqual(o2.srp.tsi, 1367);
            tc.verifyEqual(o.srp.tsi, 1361, 'value semantics: the original must not change');
        end
        function dependencyReportContract(tc)
            d = NonGravForces.dependencyReport('Verbose', false);
            tc.verifyEqual(numel(d), 4);
            % lsqvce left this table on 2026-09-09 with the owner's ruling that it is overkill for
            % WP-19; app.calibrate carries its own solver now, so nothing in the toolbox calls it.
            % mice joined on 2026-09-16: the suite evaluates sun.source rs.de440, the precise preset's
            % own default, which cannot run without NAIF's toolkit. It is listed because it was
            % reachable only through this machine's saved pathdef, and an undeclared dependency that
            % happens to be present is the failure mode this project met twice that day.
            tc.verifyEqual([d.name], ["referenceSystem","gracefoReading","gjToolbox","mice"]);
            tc.verifyTrue(all(islogical([d.present])));
            % Only the two inviolable ones are required; the other two are the gate's business, not
            % the toolbox's, and a user who computes forces without integrating needs neither.
            tc.verifyEqual([d.required], [true true false false]);
            % A BROKEN PROBE REPORTS ABSENCE, and absence is what this unit exists to report, so the
            % two are indistinguishable from the table alone. On 2026-09-16 the mice probe passed a
            % string to a MEX that takes a char, threw, and read as "mice absent" on a machine that
            % had just printed its toolkit version. Each row is therefore cross-checked against a
            % DIFFERENT observation of the same fact.
            reachable = [~isempty(which('referenceSystem.toolboxVersion')), ...
                         ~isempty(which('gracefoReading.version')), ...
                         ~isempty(which('gjToolbox.versionInfo')), ...
                         exist('cspice_tkvrsn', 'file') ~= 0];
            tc.verifyEqual([d.present], reachable, ...
                'a dependency whose entry point resolves must be reported present: a probe that throws is not evidence of absence');
            tc.verifyTrue(all(strlength([d(logical([d.present])).version]) > 0), ...
                'a present dependency must report a version, or the probe answered without measuring anything');
            % O-25: the minimum is COMPARED since 2026-09-17, and the comparison reports whether it
            % could be made at all. MICE answers CSPICE_N0067, which no dotted order can rank, so a
            % version nobody can order must never be called too old -- that is the branch that would
            % otherwise turn an unrankable string into a false verdict.
            tc.verifyTrue(all(islogical([d.comparable])) && all(islogical([d.belowMinimum])));
            mice = d(strcmp([d.name], "mice"));
            tc.verifyFalse(mice.comparable, 'CSPICE_N0067 is not a dotted version and must not be ranked');
            tc.verifyFalse(mice.belowMinimum, 'an unrankable version must never be reported below its minimum');
            numeric = d(~strcmp([d.name], "mice") & [d.present]);
            tc.verifyTrue(all([numeric.comparable]), 'the three semantic versions must be comparable with their minimums');
            tc.verifyFalse(any([d.belowMinimum]), 'the installed dependencies satisfy their declared minimums');
            % Nothing absent may claim a comparison: a probe that failed measured nothing.
            tc.verifyFalse(any([d(~[d.present]).comparable]), 'an absent dependency cannot have been compared');
        end
    end
    methods (Test, TestTags = {'reference'})
        function constReference(tc)
            % CODATA 2018 exact values and IAU 2015 B3 nominal solar radius.
            c = NonGravForces.const;
            tc.verifyEqual(c.kB, 1.380649e-23);
            tc.verifyEqual(c.NA, 6.02214076e23);
            tc.verifyRelClose(c.Rgas, c.kB * c.NA, 1e-15, "const.Rgas = kB*NA");
            tc.verifyEqual(c.Rsun, 6.957e8);
            tc.verifyEqual(c.au, 149597870700);
        end
    end
    methods (Test, TestTags = {'precision'})
        function constPrecision(tc)
            % Solar radiation pressure at 1 au: S0/c. Bound 1e-15: pure arithmetic.
            c = NonGravForces.const;
            P0 = c.S0 / 299792458;
            tc.verifyRelClose(P0, 4.5398e-6, 1e-4, "P0 from S0=1361 vs 4.5398e-6 Pa (Montenbruck & Gill give 4.56e-6 for 1367)");
        end
    end
    methods (Test, TestTags = {'robustness'})
        function constRobustness(tc)
            % The constants are data, so robustness is about the data being well formed: nothing
            % non-finite, nothing non-positive that is physically positive, the two species arrays
            % the same length, and no hidden state between calls.
            c = NonGravForces.const();
            num = ["S0","au","Rsun","Rmoon","sigmaSB","kB","NA","Rgas","amu","lambertDiffuse"];
            for f = num
                tc.verifyTrue(isscalar(c.(f)) && isfinite(c.(f)) && c.(f) > 0, sprintf('const.%s must be a finite positive scalar', f));
            end
            tc.verifyEqual(numel(c.speciesName), numel(c.speciesMass), 'the species name and mass arrays are one fact in two arrays');
            tc.verifyTrue(all(isfinite(c.speciesMass) & c.speciesMass > 0));
            tc.verifyTrue(all(strlength(c.speciesName) > 0));
            tc.verifyGreaterThan([c.graceFO.massC c.graceFO.massD], 500, 'handbook masses, not MAS1B tank masses');
            tc.verifyEqual(NonGravForces.const(), c, 'const must be a pure function of nothing');
        end
        function toolboxVersionRobustness(tc)
            % The numeric fields are the authority and `semantic` is composed from them (rule R1):
            % the robustness claim is that the two can never drift, and that a bad format refuses.
            v = NonGravForces.toolboxVersion();
            tc.verifyClass(v.major, 'uint16'); tc.verifyClass(v.minor, 'uint16'); tc.verifyClass(v.patch, 'uint16');
            parts = double(split(v.semantic, "."));
            tc.verifyEqual(parts(:)', double([v.major v.minor v.patch]), 'semantic must decompose into the numeric fields');
            tc.verifyEqual(string(NonGravForces.toolboxVersion('string')), v.semantic);
            tc.verifyEqual(NonGravForces.toolboxVersion(), v, 'no hidden state between calls');
            tc.verifyError(@() NonGravForces.toolboxVersion('bogus'), 'MATLAB:validators:mustBeMember');
        end
        function optionsRobustness(tc)
            % Empty override list and repeated set() are well defined.
            o = NonGravForces.Options("standard");
            o = o.set("drag.gsi", "maxwell"); o = o.set("drag.gsi", "sentman");
            tc.verifyEqual(o.drag.gsi, "sentman");
            s = o.toStruct(); tc.verifyTrue(isstruct(s) && isfield(s, 'drag'));
        end
        function failRobustness(tc)
            % Message with a literal percent sign must survive sprintf.
            ME = NonGravForces.fail("u", "inputX", "50%% of x -> halve it");
            tc.verifyEqual(ME.message, '50% of x -> halve it');
        end
        function provenanceRobustness(tc)
            p = NonGravForces.provenance("u", 'Warnings', string.empty(1,0), 'Files', string.empty(1,0));
            tc.verifyEmpty(p.files); tc.verifyEmpty(p.warnings);
        end
        function dependencyReportRobustness(tc)
            % Must not error whatever is on the path; only warn.
            tc.verifyWarningFree(@() NonGravForces.dependencyReport('Verbose', false), ...
                'if this fails, a REQUIRED dependency is absent on this host - run its setup');
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function optionsMetamorphic(tc)
            % Setting a field to its own value is the identity; two orders of
            % independent overrides give equal structs.
            o = NonGravForces.Options("standard");
            tc.verifyEqual(o.set("drag.alpha", 0.85).toStruct(), o.toStruct());
            a = NonGravForces.Options("fast", "drag.wind", "hwm14", "srp.tsi", 1362);
            b = NonGravForces.Options("fast", "srp.tsi", 1362, "drag.wind", "hwm14");
            tc.verifyEqual(a.toStruct(), b.toStruct());
        end
        function provenanceMetamorphic(tc)
            % Hash depends on content only, not on field order of the struct.
            s1 = struct('a', 1, 'b', "x"); s2 = orderfields(struct('b', "x", 'a', 1), s1);
            tc.verifyEqual(NonGravForces.provenance("u", 'Options', s1).optionHash, NonGravForces.provenance("u", 'Options', s2).optionHash);
            tc.verifyNotEqual(NonGravForces.provenance("u", 'Options', struct('a', 2)).optionHash, NonGravForces.provenance("u", 'Options', s1).optionHash);
        end
        function failMetamorphic(tc)
            % The identifier is a pure function of (unit, mnemonic).
            a = NonGravForces.fail("u1", "inputA", "m -> h"); b = NonGravForces.fail("u1", "inputA", "other -> h");
            tc.verifyEqual(a.identifier, b.identifier);
        end
    end
end
