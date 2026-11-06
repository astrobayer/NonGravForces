classdef tNgfFetchCatalog < NgfTestCase
%TNGFFETCHCATALOG Tests of NonGravForces.fetchCatalog and of setupNonGravForces' 'Items' against it.
%   The catalogue owns the ids, the spans and the expansion of a selection; setup fetches what it
%   expands. The metamorphic test compares the two per id, without touching the network
%   (every Fetch setting "never"), so that neither can drift from the other unnoticed.
% Created 2026-11-06 by Claude Opus 5.5 (Anthropic).
    properties (Constant)
        % A fixed 'Today': "all" and "in the future" are judged against it, so the counts below
        % do not move with the calendar.
        Today = datetime(2026, 11, 6)
    end

    methods (Test, TestTags = {'contract'})
        function fetchCatalogContract(tc)
            [T, plan] = NonGravForces.fetchCatalog('Today', tc.Today);
            tc.verifyEqual(width(T), 20);
            tc.verifyEqual(string(T.Properties.VariableNames), ["id" "group" "description" "setting" "parameter" ...
                "destination" "unit" "first" "last" "ongoing" "download_MB" "disk_MB" "sizeBasis" "licence" ...
                "licenceText" "notice" "source" "selected" "count" "totalDisk_MB"]);
            tc.verifyEqual(numel(unique(T.id)), height(T), 'ids must be unique');
            tc.verifyClass(T.first, 'datetime');
            tc.verifyTrue(all(ismember(T.licence, ["open" "academic" "unstated" "unreviewed"])));
            tc.verifyTrue(all(ismember(T.destination, ["toolbox" "pool"])));
            tc.verifyTrue(all(T.disk_MB > 0 & T.download_MB > 0));
            tc.verifyTrue(all(T.selected), 'no Items selects every id');
            % THE LICENCE CLASS IS WHAT THIRD-PARTY-NOTICES.md RECORDS: every named heading exists,
            % and only "unreviewed" has none (a class without a record would be a declaration, R3).
            here = fileparts(mfilename('fullpath'));
            notices = string(splitlines(fileread(fullfile(fileparts(here), 'THIRD-PARTY-NOTICES.md'))));
            headings = strtrim(extractAfter(notices(startsWith(notices, "## ") | startsWith(notices, "### ")), " "));
            named = T.notice(T.notice ~= "");
            for k = 1:numel(named)
                tc.verifyTrue(any(headings == named(k)), "THIRD-PARTY-NOTICES.md has no heading """ + named(k) + """");
            end
            tc.verifyEqual(T.notice == "", T.licence == "unreviewed");
            % THE SETTING IS ONE OF setupNonGravForces' OWN OPTIONS, read from its arguments block.
            setupText = fileread(fullfile(fileparts(here), 'setupNonGravForces.m'));
            for s = unique(T.setting)'
                tc.verifyTrue(contains(setupText, "options." + s + " (1,1) string"), s + " is not an option of setupNonGravForces");
            end
            tc.verifyEqual(sort(string(fieldnames(plan)))', sort(["PoolDays" "PoolMonths" "SwarmMonths" "DensityMonths" ...
                "Items" "toleos" "itsgSpans" "swarmFirst" "l1bFirst" "lastComplete" "today"]));
            % every documented identifier fires
            f = @(varargin) NonGravForces.fetchCatalog('Today', tc.Today, varargin{:});
            tc.verifyErrorId(@() f('PoolDays', "2020/01/01"), "NonGravForces:fetchCatalog:inputSyntax");
            tc.verifyErrorId(@() f('PoolDays', "2021-02-29"), "NonGravForces:fetchCatalog:inputSyntax");
            tc.verifyErrorId(@() f('PoolMonths', "2020-03..2020-01"), "NonGravForces:fetchCatalog:inputRange");
            tc.verifyErrorId(@() f('PoolMonths', "2017-12"), "NonGravForces:fetchCatalog:inputSpan");
            tc.verifyErrorId(@() f('SwarmMonths', "2013-12"), "NonGravForces:fetchCatalog:inputSpan");
            tc.verifyErrorId(@() f('DensityMonths', "2026-08"), "NonGravForces:fetchCatalog:inputSpan");
            tc.verifyErrorId(@() f('PoolDays', "2026-11-07"), "NonGravForces:fetchCatalog:inputSpan");
            tc.verifyErrorId(@() f('Items', "nope"), "NonGravForces:fetchCatalog:inputItems");
            % setup refuses a selected id whose setting is "never", before any network action
            tc.verifyErrorId(@() setupNonGravForces('Verbose', false, 'Items', "groops-grid", ...
                'FetchMsis', "never", 'FetchData', "never", 'FetchKernels', "never", 'FetchPool', "never"), ...
                "setupNonGravForces:itemsNever");
            % determinism
            tc.verifyEqual(NonGravForces.fetchCatalog('Today', tc.Today), T);
        end
    end

    methods (Test, TestTags = {'reference'})
        function fetchCatalogReference(tc)
            % The spans against an authority not written here: the HAPI headers of the committed
            % TOLEOS and Swarm excerpts (startDate, stopDate), as the server returned them.
            [~, plan] = NonGravForces.fetchCatalog('Today', tc.Today);
            here = fileparts(mfilename('fullpath'));
            files = dir(fullfile(fileparts(here), 'data', 'samples', 'toleos', '*.csv'));
            tc.assertEqual(numel(files), 3, 'one committed excerpt per TOLEOS dataset');
            for k = 1:numel(files)
                h = fileread(fullfile(files(k).folder, files(k).name));
                ds = string(regexp(h, '"x_dataset":\s*"([^"]+)"', 'tokens', 'once'));
                t0 = string(regexp(h, '"startDate":\s*"([^"Z]+)Z?"', 'tokens', 'once'));
                t1 = string(regexp(h, '"stopDate":\s*"([^"Z]+)Z?"', 'tokens', 'once'));
                row = plan.toleos([plan.toleos.dataset] == ds);
                tc.assertNumElements(row, 1, ds + " is not in the catalogue");
                tc.verifyEqual(row.start, t0, ds + " start");
                tc.verifyEqual(row.stop, t1, ds + " stop");
            end
            files = dir(fullfile(fileparts(here), 'data', 'samples', 'swarm', '*.csv'));
            tc.assertNotEmpty(files);
            for k = 1:numel(files)
                h = fileread(fullfile(files(k).folder, files(k).name));
                t0 = string(regexp(h, '"startDate":\s*"(\d{4}-\d{2}-\d{2})', 'tokens', 'once'));
                tc.verifyEqual(plan.swarmFirst, datetime(t0, 'InputFormat', 'yyyy-MM-dd'), files(k).name);
            end
        end
    end

    methods (Test, TestTags = {'precision'})
        function fetchCatalogPrecision(tc)
            % Counts derived by hand for one selection; each line says how.
            T = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', "2020-01-01..2020-01-31", ...
                'PoolMonths', "2020-01", 'SwarmMonths', "2019-12", 'DensityMonths', "2008-01");
            c = @(id) T.count(T.id == id);
            tc.verifyEqual(c("itsg"), 31 * 2 * 3);    % 31 days, two satellites, three products
            tc.verifyEqual(c("dst"), 2);              % 2020-01 and the month before it
            tc.verifyEqual(c("f30"), 2);              % the adjusted and the absolute file
            tc.verifyEqual(c("swarm"), 2 * 3);        % 2019-12 and 2020-01, three satellites
            tc.verifyEqual(c("toleos"), 1 + 2);       % 2020-01: GRACE-FO 1; 2008-01: CHAMP and GRACE-A
            tc.verifyEqual(c("pool-l1b"), 31);
            tc.verifyEqual([c("pool-spaceweather") c("pool-eop")], [1 1]);
            tc.verifyEqual([c("ceres") c("groops-syn1deg") c("groops-grid") c("de440")], [3 2 2 1]);
            tc.verifyEqual(T.totalDisk_MB(T.id == "itsg"), 186 * T.disk_MB(T.id == "itsg"));
            % "all" months: 2018-06..2026-10 on Today = 2026-11-06 is 7 + 7*12 + 10 = 101 months
            [T, plan] = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolMonths', "all");
            tc.verifyNumElements(plan.PoolMonths, 101);
            tc.verifyEqual([plan.PoolMonths(1) plan.PoolMonths(end)], ["2018-06" "2026-10"]);
            % days 2018-06-01..2026-10-31: 2018 214, 2019..2025 7 years with two leap days, 2026 304
            tc.verifyEqual(T.count(T.id == "pool-l1b"), 214 + 7 * 365 + 2 + 304);
            % "all" reference days: TU Graz's two spans and not the gap between them
            [~, plan] = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', "all");
            tc.verifyEqual(plan.PoolDays([1 end]), ["2002-04-04" "2026-10-31"]);
            tc.verifyFalse(any(plan.PoolDays == "2018-01-01"), 'the GRACE / GRACE-FO gap is not published');
            tc.verifyTrue(any(plan.PoolDays == "2017-06-29") && any(plan.PoolDays == "2018-06-01"));
        end
    end

    methods (Test, TestTags = {'robustness'})
        function fetchCatalogRobustness(tc)
            % no selection: the dated sources count nothing, the undated keep their fixed files
            [T, plan] = NonGravForces.fetchCatalog('Today', tc.Today);
            dated = ["dst" "f30" "itsg" "swarm" "toleos" "pool-l1b" "pool-spaceweather" "pool-eop"];
            tc.verifyEqual(T.count(ismember(T.id, dated)), zeros(numel(dated), 1));
            tc.verifySize(plan.PoolDays, [1 0]);
            % whitespace, duplicates, overlapping ranges and a single day inside a range unify
            [~, a] = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', [" 2020-01-01..2020-01-05 ", ...
                "2020-01-03..2020-01-07", "2020-01-04", "2020-01-07"]);
            tc.verifyEqual(a.PoolDays, "2020-01-0" + (1:7));
            % a leap day is a day
            [~, b] = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', "2020-02-29");
            tc.verifyEqual(b.PoolDays, "2020-02-29");
            % a zoned Today is the same day
            [~, z] = NonGravForces.fetchCatalog('Today', datetime(2026, 11, 6, 12, 0, 0, 'TimeZone', 'UTC'));
            tc.verifyEqual(z.today, tc.Today);
            % Items: the rest stays in the table, unselected and costing nothing
            T = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', "2020-01-01", 'Items', ["itsg" "dst"]);
            tc.verifyEqual(sort(T.id(T.selected))', ["dst" "itsg"]);
            tc.verifyEqual(T.totalDisk_MB(~T.selected), zeros(nnz(~T.selected), 1));
            tc.verifyGreaterThan(T.count(T.id == "f30"), 0, 'count is what the selection implies, selected or not');
            % an empty string array element is not a date
            tc.verifyErrorId(@() NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', ""), ...
                "NonGravForces:fetchCatalog:inputSyntax");
        end
    end

    methods (Test, TestTags = {'metamorphic'})
        function fetchCatalogMetamorphic(tc)
            % (1) a range is the list of its days, and splitting it changes nothing
            [Ta, a] = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', "2020-01-30..2020-02-02");
            [Tb, b] = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', ["2020-02-02" "2020-01-31" "2020-01-30" "2020-02-01"]);
            [Tc, c] = NonGravForces.fetchCatalog('Today', tc.Today, 'PoolDays', ["2020-01-30..2020-01-31" "2020-02-01..2020-02-02"]);
            tc.verifyEqual(b, a); tc.verifyEqual(c, a);
            tc.verifyEqual(Tb, Ta); tc.verifyEqual(Tc, Ta);
            % (2) THE CATALOGUE AND setupNonGravForces AGREE, per id: the count here is the number of
            % distinct sources setup lists for the same selection. Every Fetch setting is "never",
            % so nothing reaches the network; DataPool is a fresh empty folder, so nothing is present.
            pool = string(tempname);
            mkdir(pool);
            tc.addTeardown(@() rmdir(pool, 's'));
            sel = {'PoolDays', ["2008-01-01" "2020-01-01..2020-01-03"], 'PoolMonths', "2020-01", ...
                'SwarmMonths', "2019-12", 'DensityMonths', "2008-02"};
            T = NonGravForces.fetchCatalog(sel{:});
            info = setupNonGravForces('Verbose', false, 'FetchMsis', "never", 'FetchData', "never", ...
                'FetchKernels', "never", 'FetchGroops', "never", 'FetchPool', "never", 'DataPool', pool, sel{:});
            recs = [rmfield(info.sources, {'dir', 'message'}), rmfield(info.data.items, {'file'})];
            tc.verifyTrue(all(ismember([recs.id], T.id)), 'every record carries a catalogue id');
            for k = 1:height(T)
                got = numel(unique([recs([recs.id] == T.id(k)).source]));
                tc.verifyEqual(got, T.count(k), "setup and the catalogue disagree on " + T.id(k));
            end
            tc.verifyEqual(info.plan.PoolDays, ["2008-01-01" "2020-01-01" "2020-01-02" "2020-01-03"]);
            % (3) Items narrows setup to the selection and says so for the rest. The one source
            % selected is NRLMSISE-00, which the gate has already fetched (ciGate builds from it), so
            % "ifMissing" finds it present and the network stays untouched; every data item is left
            % out, and an absent one must say why in the words of Items, not of a Fetch setting.
            info = setupNonGravForces('Verbose', false, 'DataPool', pool, 'Items', "nrlmsise00", 'PoolDays', "2020-01-01");
            tc.verifyTrue(ismember(info.sources([info.sources.id] == "nrlmsise00").action, ["present" "skipped"]));
            actions = [info.data.items.action];
            tc.verifyTrue(all(startsWith(actions, ["present" "missing (not selected by Items)" "failed"])), strjoin(unique(actions), " | "));
            tc.verifyTrue(any(actions == "missing (not selected by Items)"), 'the empty pool holds nothing, so something is reported missing');
        end
    end
end
