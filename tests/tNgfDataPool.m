classdef tNgfDataPool < NgfTestCase
%TNGFDATAPOOL Tests of NonGravForces.dataPool and of setupNonGravForces writing every data file into the pool.
%   The contract (toolboxInstaller docs/DATA_POOL.md, owner 2026-11-08): the pool root is the explicit
%   argument, else DATAPOOL, else READING_DATA_POOL, else a refusal; everything downloaded that is data
%   goes below it, laid out by product, and nothing falls back to the toolbox folder. Each test sets
%   both environment variables itself and restores them on teardown, so the host's own pool (the CI
%   runner's READING_DATA_POOL) neither helps nor hinders.
% Created 2026-11-08 by Claude Opus 5.5 (Anthropic).

    methods (Test, TestTags = {'contract'})
        function dataPoolContract(tc)
            tc.environment("E:\poolA", "E:\poolB");
            p = NonGravForces.dataPool();
            tc.verifyClass(p, 'string');
            tc.verifySize(p, [1 1]);
            tc.verifyEqual(p, "E:\poolA", 'DATAPOOL comes before READING_DATA_POOL');
            tc.verifyEqual(NonGravForces.dataPool("F:\mine"), "F:\mine", 'the explicit argument comes first');
            tc.environment("", "E:\poolB");
            tc.verifyEqual(NonGravForces.dataPool(), "E:\poolB", 'READING_DATA_POOL is the fallback');
            tc.environment("", "");
            tc.verifyErrorId(@() NonGravForces.dataPool(), "NonGravForces:dataPool:dataPoolUnset");
            tc.verifyErrorId(@() NonGravForces.dataPool("", 'Option', "Pool"), "NonGravForces:dataPool:dataPoolUnset");
            tc.verifyEqual(NonGravForces.dataPool("", 'Required', false), "", 'Required false names no pool rather than refusing');
            tc.verifyEqual(NonGravForces.dataPool("F:\mine", 'Required', false), "F:\mine");
        end
    end

    methods (Test, TestTags = {'robustness'})
        function dataPoolRobustness(tc)
            % blanks are no name; surrounding blanks are not part of a name; inner ones are
            tc.environment("   ", "E:\poolB");
            tc.verifyEqual(NonGravForces.dataPool(), "E:\poolB");
            tc.verifyEqual(NonGravForces.dataPool("  F:\my pool  "), "F:\my pool");
            tc.verifyEqual(NonGravForces.dataPool('F:\char'), "F:\char", 'a char argument is a name too');
            tc.environment(" ", sprintf('\t'));
            tc.verifyErrorId(@() NonGravForces.dataPool("  "), "NonGravForces:dataPool:dataPoolUnset");
            % a pool that does not exist yet is still the pool: setup creates the product folders
            missing = string(tempname);
            tc.verifyEqual(NonGravForces.dataPool(missing), missing);
            tc.verifyFalse(isfolder(missing), 'resolving a name creates nothing');
        end
    end

    methods (Test, TestTags = {'metamorphic'})
        function dataPoolMetamorphic(tc)
            % (1) HOW THE ROOT IS NAMED CHANGES NOTHING: 'DataPool', DATAPOOL and READING_DATA_POOL
            % naming the same folder give setup the same root and the same records. (2) Every data
            % file setup keeps lies below that root, in the product folder of the contract, and none
            % below the toolbox. (3) Without any name, every data step is skipped -- with every Fetch
            % setting "ifMissing", so a fallback to the toolbox would have had every chance -- and
            % the toolbox's data folder is left exactly as it was.
            pool = string(tempname);
            mkdir(pool);
            tc.addTeardown(@() rmdir(pool, 's'));
            % placeholders where the contract puts each product (presence is by exact file; the
            % kernel's is then refused by its digest, which proves it was looked for in SPICE)
            day = fullfile(pool, "GRACEFO", "L1B", "RL04", "2020", "2020-01-01");
            ceres = ["rsut" "rlut" "rsdt"] + "_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_200003-202509.nc";
            files = [fullfile(pool, "EOP", "finals2000A.all"), fullfile(pool, "SpaceWeather", "SW-All.csv"), ...
                fullfile(pool, "IGRF", "igrf14coeffs.txt"), fullfile(day, ["GNV1B" "SCA1B" "MAG1B"] + "_2020-01-01_C_04.txt"), ...
                fullfile(pool, "CERES", "EBAF", ceres), fullfile(pool, "SPICE", "de440s.bsp")];
            for f = files
                if ~isfolder(fileparts(f)), mkdir(fileparts(f)); end
                fid = fopen(f, 'w'); fwrite(fid, 'placeholder'); fclose(fid);
            end
            never = {'Verbose', false, 'FetchMsis', "never", 'FetchData', "never", 'FetchKernels', "never", ...
                'FetchGroops', "never", 'FetchPool', "never"};
            ids = ["eop" "spaceweather" "igrf" "l1b-reference" "ceres" "de440"];
            tc.environment("", "");
            a = setupNonGravForces(never{:}, 'DataPool', pool);
            tc.environment(pool, "");
            b = setupNonGravForces(never{:});
            tc.environment("", pool);
            c = setupNonGravForces(never{:});
            tc.verifyEqual([a.dataPool b.dataPool c.dataPool], [pool pool pool]);
            tc.verifyEqual(b.data.items, a.data.items, 'DATAPOOL and DataPool give the same records');
            tc.verifyEqual(c.data.items, a.data.items, 'READING_DATA_POOL and DataPool give the same records');
            recs = a.data.items(ismember([a.data.items.id], ids));
            tc.verifyEqual(numel(recs), 1 + 1 + 1 + 1 + 3 + 1, 'one record per file, the Level-1B day one record');
            root = string(fileparts(fileparts(mfilename('fullpath'))));
            for r = recs(ismember([recs.id], ids(1:5)))
                tc.verifyEqual(r.action, "present", r.name + " was not found where the contract puts it");
                kept = splitlines(r.file);
                tc.verifyTrue(all(startsWith(kept, pool)), r.name + " is not below the pool");
                tc.verifyFalse(any(startsWith(kept, root)), r.name + " is below the toolbox");
            end
            tc.verifyEqual(extractBefore(recs([recs.id] == "l1b-reference").file, newline), fullfile(day, "GNV1B_2020-01-01_C_04.txt"));
            k = recs([recs.id] == "de440");
            tc.verifyTrue(startsWith(k.action, "failed: the file on disk has SHA-256"), ...
                "the placeholder in <pool>/SPICE was not the file the kernel step judged: " + k.action);
            % (3) no name at all
            before = localListing(fullfile(root, 'data'));
            tc.environment("", "");
            d = setupNonGravForces('Verbose', false, 'FetchMsis', "never", 'FetchGroops', "never", ...
                'FetchData', "ifMissing", 'FetchKernels', "ifMissing", 'FetchPool', "ifMissing");
            tc.verifyEqual(d.dataPool, "");
            skip = "skipped: DATAPOOL is not set -> set the DATAPOOL environment variable or pass 'DataPool'";
            tc.verifyTrue(all([d.data.items.action] == skip), strjoin(unique([d.data.items.action]), " | "));
            tc.verifyTrue(all([d.data.items.file] == ""));
            tc.verifyTrue(all(ismember(ids, [d.data.items.id])), 'every data source says it was skipped');
            tc.verifyEqual(localListing(fullfile(root, 'data')), before, 'the toolbox data folder changed');
        end
    end

    methods (Access = private)
        function environment(tc, datapool, reading)
            % Set DATAPOOL and READING_DATA_POOL; the first call of a test registers their restoration.
            if isempty(tc.saved)
                tc.saved = {getenv('DATAPOOL'), getenv('READING_DATA_POOL')};
                saved = tc.saved;
                tc.addTeardown(@() setenv('DATAPOOL', saved{1}));
                tc.addTeardown(@() setenv('READING_DATA_POOL', saved{2}));
            end
            setenv('DATAPOOL', char(datapool));
            setenv('READING_DATA_POOL', char(reading));
        end
    end

    properties (Access = private)
        saved = {}
    end
end

function L = localListing(folder)
%LOCALLISTING Every file below folder with its size and date: what "left exactly as it was" compares.
d = dir(fullfile(folder, '**', '*'));
d = d(~[d.isdir]);
L = sort(string(fullfile({d.folder}, {d.name})) + " " + string([d.bytes]) + " " + string([d.datenum]));
end
