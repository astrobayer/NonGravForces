classdef tNgfIo < NgfTestCase
%TNGFIO Tests of the io layer (LP-2): readACC1B, readSCA1B, readGNV1B, readTHR1B, readMAS1B, readMAG1B, igrfCoefficients, importCERES, massModel, spaceWeather.
%
%   Oracle: the real 2020-01-01 GRACE-C RL04 products in data/samples (first 600
%   records of the 1 Hz products; THR1B and MAS1B complete), the IERS
%   finals2000A excerpt in data/eop, and the CelesTrak / SET excerpts in
%   data/spaceweather (HANDOVER Part I, 2026-09-06). Independent oracles are
%   the raw file text (parsed here with regexp, not with gracefoReading), the
%   IAU 1982 GMST polynomial, the IERS leap-second count and JPL D-56935.
%   importCERES: the CERES EBAF 4.2.1 excerpt in data/ceres (December 2019 and
%   January 2020, values unchanged) against tools/pyxcheck/ceres_check.py, a
%   standard-library netCDF-3 reader written from the format specification.
%   Every bound carries the number it was measured at on 2026-09-06.
%
%   No test skips when a fixture is missing: it fails (a skip is `incomplete`
%   and not evidence).
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
    properties
        f      % fixture paths
        eop    % referenceSystem.readEop of data/eop/finals2000A excerpt
        tmp    % scratch folder for corrupted copies
    end
    methods (TestClassSetup)
        function loadFixtures(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            S = fullfile(root, 'data', 'samples'); W = fullfile(root, 'data', 'spaceweather');
            tc.f = struct( ...
                'act', fullfile(S, 'ACT1B_2020-01-01_C_04.txt'), 'sca', fullfile(S, 'SCA1B_2020-01-01_C_04.txt'), ...
                'gnv', fullfile(S, 'GNV1B_2020-01-01_C_04.txt'), 'thr', fullfile(S, 'THR1B_2020-01-01_C_04.txt'), ...
                'mas', fullfile(S, 'MAS1B_2020-01-01_C_04.txt'), 'mag', fullfile(S, 'MAG1B_2020-01-01_C_04.txt'), ...
                'igrf', fullfile(root, 'data', 'igrf', 'igrf14coeffs.txt'), ...
                'rsut', fullfile(root, 'data', 'ceres', 'rsut_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc'), ...
                'rlut', fullfile(root, 'data', 'ceres', 'rlut_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc'), ...
                'rsdt', fullfile(root, 'data', 'ceres', 'rsdt_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc'), ...
                'ceresCells', fullfile(root, 'tools', 'pyxcheck', 'ceres_cells.csv'), ...
                'sw', fullfile(W, 'SW-All_2019-09_2020-03.csv'), ...
                'solfsmy', fullfile(W, 'SOLFSMY_2019-09_2020-03.TXT'), 'dtc', fullfile(W, 'DTCFILE_2019-09_2020-03.TXT'), ...
                'eop', fullfile(root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            for n = string(fieldnames(tc.f))'
                tc.assertTrue(isfile(tc.f.(n)), sprintf('fixture %s is missing: %s', n, tc.f.(n)));
            end
            tc.eop = referenceSystem.readEop(tc.f.eop);
            tc.tmp = tempname; mkdir(tc.tmp);
            tc.addTeardown(@() rmdir(tc.tmp, 's'));
        end
    end
    methods
        function p = corrupt(tc, src, name, editFcn)
            %CORRUPT Copy an L1B/text fixture with editFcn(lines, kFirstData) applied; return the path.
            lines = splitlines(string(fileread(src)));
            k = find(startsWith(lines, "# End of YAML header"), 1) + 1;
            lines = editFcn(lines, k);
            p = fullfile(tc.tmp, name);
            fid = fopen(p, 'w'); fprintf(fid, '%s\n', lines); fclose(fid);
        end
        function tok = rawRecord(tc, src, k)
            %RAWRECORD Whitespace tokens of the k-th data record, parsed independently of gracefoReading.
            lines = splitlines(string(fileread(src)));
            k0 = find(startsWith(lines, "# End of YAML header"), 1);
            tok = split(strtrim(lines(k0 + k)));
            tc.assertGreaterThan(numel(tok), 3);
        end
        function tr = range(~, s0, s1)
            %RANGE decimal-MJD TimeRange for seconds s0..s1 of 2020-01-01 (GPS day 58849).
            tr = [58849 + s0 / 86400, 58849 + s1 / 86400];
        end
        function f = ceres(tc)
            %CERES The three EBAF excerpt files, in the order rsut rlut rsdt.
            f = string({tc.f.rsut, tc.f.rlut, tc.f.rsdt});
        end
        function q = speedRatio(~, fcn, base)
            %SPEEDRATIO median over 5 repeats of (fcn time / base time), after one warm-up of each.
            fcn(); base();
            r = zeros(1, 5);
            for k = 1:5
                t1 = tic; fcn(); a = toc(t1); t2 = tic; base(); b = toc(t2); r(k) = a / b;
            end
            q = median(r);
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function importCERESContract(tc)
            [cer, aux] = NonGravForces.io.importCERES(tc.ceres());
            tc.verifyEqual(cer.product, "CERES-EBAF-4-2-1"); tc.verifyEqual(cer.version, "4.2.1");
            tc.verifySize(cer.rsut, [180 360 2]); tc.verifySize(cer.rlut, [180 360 2]); tc.verifySize(cer.rsdt, [180 360 2]);
            tc.verifyClass(cer.rsut, 'single'); tc.verifySize(cer.time, [1 2]); tc.verifySize(cer.mjd, [1 2]);
            tc.verifyEqual(cer.lat_deg([1 end]), [-89.5 89.5]); tc.verifyEqual(cer.lon_deg([1 end]), [0.5 359.5]);
            tc.verifyEqual(cer.latitude, "geodetic"); tc.verifyEqual(cer.referenceHeight_m, 20e3);
            tc.verifyTrue(startsWith(cer.doi, "10.5067/")); tc.verifyTrue(contains(cer.license, "Creative Commons"));
            tc.verifyEqual(aux.files, tc.ceres()); tc.verifyEqual(aux.timeIndex, [1 2]); tc.verifyEqual(aux.provenance.unit, "io.importCERES");
            C = @NonGravForces.io.importCERES; f = tc.ceres();
            tc.verifyErrorId(@() C(f(1:2)), 'NonGravForces:importCERES:inputFiles');
            tc.verifyErrorId(@() C(f([1 1 3])), 'NonGravForces:importCERES:inputFiles');
            tc.verifyErrorId(@() C([f(1:2) string(fullfile(tc.tmp, 'nothere.nc'))]), 'NonGravForces:importCERES:dataMissing');
            other = fullfile(tc.tmp, 'other_rsdt.nc'); copyfile(tc.f.rsdt, other, 'f'); fileattrib(other, '+w');
            ncwriteatt(other, '/', 'source_id', 'NOT-CERES');
            tc.verifyErrorId(@() C([f(1:2) string(other)]), 'NonGravForces:importCERES:dataFormat');
            tc.verifyErrorId(@() C(f, 'Period', [datetime(2020, 1, 5) datetime(2020, 1, 1)]), 'NonGravForces:importCERES:inputPeriod');
            tc.verifyErrorId(@() C(f, 'Period', [datetime(2019, 12, 1) datetime(2020, 1, 1)]), 'NonGravForces:importCERES:inputPeriod');
            tc.verifyErrorId(@() C(f, 'Period', [datetime(2020, 1, 1) datetime(2020, 2, 1)]), 'NonGravForces:importCERES:inputPeriod');
        end
        function igrfCoefficientsContract(tc)
            [c, aux] = NonGravForces.io.igrfCoefficients(tc.f.igrf);
            tc.verifyEqual(c.nMax, 13);
            tc.verifyEqual(c.svMax, 8);
            tc.verifyEqual(numel(c.epochs), 26);
            tc.verifyEqual(c.epochs(1), 1900);
            tc.verifyEqual(c.epochs(end), 2025);
            tc.verifySize(c.gTable, [14 14 26]);
            tc.verifySize(c.svG, [14 14]);
            tc.verifyEmpty(c.g, 'without an Epoch there is no interpolated field to return');
            tc.verifyEqual(aux.generation, 14);
            tc.verifyEqual(aux.numCoefficients, 195);
            tc.verifyEqual(aux.validTo, 2030);
            c2 = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2020.0);
            tc.verifySize(c2.g, [14 14]);
            tc.verifyEqual(c2.year, 2020);

            I = @NonGravForces.io.igrfCoefficients;
            tc.verifyErrorId(@() I(fullfile(tc.tmp, 'nothere.txt')), 'NonGravForces:igrfCoefficients:dataMissing');
            tc.verifyErrorId(@() I(tc.f.igrf, 'Epoch', 1899), 'NonGravForces:igrfCoefficients:inputRange');
            tc.verifyErrorId(@() I(tc.f.igrf, 'Epoch', 2031), 'NonGravForces:igrfCoefficients:inputRange');
            tc.verifyErrorId(@() I(tc.f.mag), 'NonGravForces:igrfCoefficients:dataFormat');
        end
        function readMAG1BContract(tc)
            [ep, fld, aux] = NonGravForces.io.readMAG1B(tc.f.mag);
            K = numel(ep.sec);
            tc.verifyEqual(K, 600);
            tc.verifyEqual(ep.scale, "GPS");
            tc.verifySize(fld.B_srf_T, [3 K]);
            tc.verifySize(fld.torquerBranch_A, [6 K]);
            tc.verifySize(fld.torquerNet_A, [3 K]);
            tc.verifySize(fld.magnitude_T, [1 K]);
            tc.verifyEqual(aux.torquerBranchNames, ["torque1A","torque2A","torque3A","torque1B","torque2B","torque3B"]);
            tc.verifySize(aux.calibration.magField, [3 K]);
            tc.verifySize(aux.qualityBits, [8 K]);
            tc.verifyEqual(aux.product, "MAG1B");
            tc.verifyEqual(aux.satellite, "C");
            tc.verifyErrorId(@() NonGravForces.io.readMAG1B(tc.f.mas), 'NonGravForces:readMAG1B:dataProduct');
            tc.verifyErrorId(@() NonGravForces.io.readMAG1B(tc.f.mag, 'TimeRange', [58851 58852]), 'NonGravForces:readMAG1B:dataEmpty');
        end
        function readACC1BContract(tc)
            [ep, a, aux] = NonGravForces.io.readACC1B(tc.f.act);
            tc.verifySize(a, [3 600]); tc.verifyEqual(string(ep.scale), "GPS"); tc.verifySize(ep.day, [600 1]);
            tc.verifyEqual(aux.product, "ACT1B"); tc.verifyEqual(aux.satellite, "C"); tc.verifyEqual(aux.numRecords, 600);
            tc.verifyEqual(aux.sampling_s, 1); tc.verifySize(aux.angular_rad_s2, [3 600]); tc.verifySize(aux.residual_m_s2, [3 600]);
            tc.verifyClass(aux.qualflgBits, 'logical'); tc.verifySize(aux.qualflgBits, [600 8]);
            tc.verifyEqual(aux.provenance.unit, "io.readACC1B"); tc.verifyEqual(aux.provenance.files, string(tc.f.act));
            tc.verifyErrorId(@() NonGravForces.io.readACC1B(tc.f.sca), 'NonGravForces:readACC1B:dataProduct');
            tc.verifyErrorId(@() NonGravForces.io.readACC1B(tc.f.act, 'TimeRange', [58850 58851]), 'NonGravForces:readACC1B:dataEmpty');
            tc.verifyError(@() NonGravForces.io.readACC1B(fullfile(tc.tmp, 'nope.txt')), 'MATLAB:validators:mustBeFile');
        end
        function readSCA1BContract(tc)
            [ep, q, aux] = NonGravForces.io.readSCA1B(tc.f.sca);
            tc.verifySize(q, [4 600]); tc.verifySize(ep.sec, [600 1]); tc.verifyEqual(string(ep.scale), "GPS");
            tc.verifyEqual(aux.frame, "BODY->GCRS"); tc.verifyEqual(aux.convention.quaternion, "scalarFirst");
            tc.verifyEqual(aux.convention.product, "hamilton"); tc.verifySize(aux.qualRss, [1 600]); tc.verifySize(aux.scaId, [1 600]);
            tc.verifyTrue(all(q(1, :) >= 0), 'canonical sign q0 >= 0');
            tc.verifyEqual(aux.product, "SCA1B"); tc.verifyEqual(aux.provenance.unit, "io.readSCA1B");
            tc.verifyErrorId(@() NonGravForces.io.readSCA1B(tc.f.act), 'NonGravForces:readSCA1B:dataProduct');
        end
        function readGNV1BContract(tc)
            [ep, r, v, aux] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop);
            tc.verifySize(r, [3 600]); tc.verifySize(v, [3 600]); tc.verifySize(ep.day, [600 1]);
            tc.verifyEqual(aux.frame, "GCRS"); tc.verifyEqual(aux.eopSeries, "finals2000A"); tc.verifyEqual(aux.rotation.path, "direct");
            tc.verifySize(aux.sigma_r_m, [3 600]); tc.verifySize(aux.sigma_v_m_s, [3 600]);
            tc.verifyGreaterThan(min(vecnorm(r)), 6.7e6); tc.verifyLessThan(max(vecnorm(r)), 7.1e6);
            [~, rI, ~, auxI] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Frame', "ITRS");
            tc.verifyEqual(auxI.frame, "ITRS"); tc.verifyEmpty(auxI.rotation); tc.verifyEqual(auxI.eopSeries, "");
            tc.verifyNotEqual(rI, r, 'ITRS and GCRS must differ');
            tc.verifyErrorId(@() NonGravForces.io.readGNV1B(tc.f.gnv), 'NonGravForces:readGNV1B:inputEop');
            tc.verifyErrorId(@() NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', struct('mjd', 1)), 'NonGravForces:readGNV1B:inputEop');
            tc.verifyErrorId(@() NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', fullfile(tc.tmp, 'no.eop')), 'NonGravForces:readGNV1B:inputEop');
            tc.verifyErrorId(@() NonGravForces.io.readGNV1B(tc.f.sca, 'Eop', tc.eop), 'NonGravForces:readGNV1B:dataProduct');
            tc.verifyError(@() NonGravForces.io.readGNV1B(tc.f.gnv, 'Frame', "TEME"), 'MATLAB:validators:mustBeMember');
        end
        function readTHR1BContract(tc)
            [T, aux] = NonGravForces.io.readTHR1B(tc.f.thr);
            tc.verifyEqual(T.Properties.VariableNames, timeSeriesPreparation.qc.eventTable().Properties.VariableNames);
            tc.verifyEqual(height(T), aux.numRecords - aux.numRecordsNoFiring);
            tc.verifyTrue(all(T.class == "impulse") && all(T.source == "THR1B") && all(T.axis == "all") && all(T.confidence == 1));
            tc.verifyTrue(all(isnan(T.amplitude)) && all(isnan(T.sigma)));
            tc.verifyEqual(T.tStart.TimeZone, 'UTC'); tc.verifyTrue(all(T.tEnd >= T.tStart));
            tc.verifyEqual(numel(aux.thruster), 14); tc.verifySize(aux.onTime_s, [14 aux.numRecords]); tc.verifySize(aux.fired, [14 aux.numRecords]);
            tc.verifyEqual(numel(aux.epoch.day), aux.numRecords); tc.verifyEqual(aux.product, "THR1B");
            tc.verifyErrorId(@() NonGravForces.io.readTHR1B(tc.f.mas), 'NonGravForces:readTHR1B:dataProduct');
            tc.verifyErrorId(@() NonGravForces.io.readTHR1B(tc.f.thr, 'TimeRange', [58851 58852]), 'NonGravForces:readTHR1B:dataEmpty');
        end
        function readMAS1BContract(tc)
            [ep, tank, aux] = NonGravForces.io.readMAS1B(tc.f.mas);
            tc.verifyEqual(numel(ep.day), 24); tc.verifySize(tank.total_kg, [1 24]); tc.verifySize(tank.tank_kg, [2 24]);
            tc.verifySize(tank.thruster_kg, [2 24]); tc.verifyTrue(all(isnan(tank.thruster_kg(:))), 'thruster columns are absent in the 2020 files');
            tc.verifyTrue(tank.present.gas_mass_tnk1 && tank.present.gas_mass_tnk2 && ~tank.present.mass_tnk);
            tc.verifyTrue(all(tank.total_kg > 0 & tank.total_kg < 100));
            tc.verifySize(aux.prodFlagBits, [24 8]); tc.verifyEqual(aux.product, "MAS1B");
            tc.verifyErrorId(@() NonGravForces.io.readMAS1B(tc.f.thr), 'NonGravForces:readMAS1B:dataProduct');
        end
        function massModelContract(tc)
            m = NonGravForces.io.massModel("C");
            tc.verifyEqual(m.kind, "constant"); tc.verifyEqual(m.satellite, "C"); tc.verifyEqual(m.span, [-Inf Inf]);
            tc.verifyClass(m.at, 'function_handle'); tc.verifyEqual(m.provenance.unit, "io.massModel");
            ep = referenceSystem.mjdToEpoch([58849.1; 58849.5; 58850], "GPS");
            tc.verifySize(m.at(ep), [1 3]);
            [epM, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            v = NonGravForces.io.massModel("D", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg);
            tc.verifyEqual(v.kind, "mas1b"); tc.verifyEqual(v.scale, "GPS"); tc.verifySize(v.at(ep), [1 3]);
            tc.verifyEqual(v.span, [epM.day(1) + epM.sec(1)/86400, epM.day(end) + epM.sec(end)/86400]);
            tc.verifyError(@() NonGravForces.io.massModel("E"), 'MATLAB:validators:mustBeMember');
            tc.verifyErrorId(@() NonGravForces.io.massModel("C", 'Anchor_kg', 30), 'NonGravForces:massModel:inputAnchor');
            tc.verifyErrorId(@() NonGravForces.io.massModel("C", 'Mode', "mas1b"), 'NonGravForces:massModel:inputMas1b');
            tc.verifyErrorId(@() NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg(1:3)), 'NonGravForces:massModel:inputMas1b');
            tc.verifyErrorId(@() m.at(struct('x', 1)), 'NonGravForces:massModel:inputEpoch');
        end
        function spaceWeatherContract(tc)
            sw = NonGravForces.io.spaceWeather(tc.f.sw);
            tc.verifyEqual(sw.format, "celestrak"); tc.verifyEqual(numel(sw.mjd), 213); tc.verifyEqual(sw.span, [58727 58939]);
            tc.verifySize(sw.index.ap3h, [213 8]); tc.verifySize(sw.index.kp3h, [213 8]); tc.verifyClass(sw.index.dataType, 'string');
            d = sw.query(referenceSystem.mjdToEpoch([58849.5; 58850.2], "UTC"));
            tc.verifySize(d.f107, [1 2]); tc.verifySize(d.apHistory, [7 2]); tc.verifySize(d.kp, [1 2]); tc.verifySize(d.dataType, [1 2]);
            sj = NonGravForces.io.spaceWeather([string(tc.f.dtc) string(tc.f.solfsmy)]);
            tc.verifyEqual(sj.format, "jb2008"); tc.verifyEqual(numel(sj.mjd), 213); tc.verifySize(sj.index.dtc, [213 24]);
            dj = sj.query(referenceSystem.mjdToEpoch(58849.5, "UTC"));
            tc.verifyEqual(sort(string(fieldnames(dj))), sort(["f10";"f10b";"s10";"s10b";"m10";"m10b";"y10";"y10b";"dtc"]));
            sj2 = NonGravForces.io.spaceWeather([string(tc.f.solfsmy) string(tc.f.dtc)], 'Format', "jb2008");
            tc.verifyEqual(sj2.mjd, sj.mjd, 'file order must not matter');
            tc.verifyErrorId(@() NonGravForces.io.spaceWeather(tc.f.dtc), 'NonGravForces:spaceWeather:inputFiles');
            tc.verifyErrorId(@() NonGravForces.io.spaceWeather(tc.f.sw, 'Format', "jb2008"), 'NonGravForces:spaceWeather:inputFormat');
            tc.verifyErrorId(@() sw.query(struct('x', 1)), 'NonGravForces:spaceWeather:inputEpoch');
            tc.verifyError(@() NonGravForces.io.spaceWeather(fullfile(tc.tmp, 'nope.csv')), 'MATLAB:validators:mustBeFile');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function importCERESReference(tc)
            % Oracle: tools/pyxcheck/ceres_check.py reads the same netCDF-3 files with its own parser (T9). Cell values
            % at eight cells per field and month, chosen at the corners, the equator, both hemispheres and the longitude
            % seam, pin the axis order and orientation: they must agree BITWISE (one float, printed with 17 digits).
            % The area-weighted global means are the whole grid in one number each.
            cer = NonGravForces.io.importCERES(tc.ceres());
            opts = setvartype(detectImportOptions(tc.f.ceresCells), {'variable', 'lat_deg', 'lon_deg'}, 'string');
            T = readtable(tc.f.ceresCells, opts);                % "global" rows keep their label
            cells = T(T.lat_deg ~= "global", :);
            tc.assertEqual(height(cells), 48, 'the oracle holds 8 cells per field and month');
            worst = 0;
            for k = 1:height(cells)
                i = find(cer.lat_deg == double(cells.lat_deg(k))); j = find(cer.lon_deg == double(cells.lon_deg(k)));
                got = double(cer.(cells.variable(k))(i, j, cells.month(k)));
                worst = max(worst, abs(got - cells.value(k)));
            end
            tc.verifyBoundAndRecord("importCERES cells vs the independent netCDF-3 reader (48 values)", worst, 0, "W/m^2", "the same float read twice: bitwise");
            tc.verifyEqual(cer.mjd, juliandate(datetime(2000, 3, 1) + days(unique(cells.time_days)')) - 2400000.5, 'the month nodes are the files'' own times');
            glob = T(T.lat_deg == "global", :);
            w = sind(cer.lat_deg' + 0.5) - sind(cer.lat_deg' - 0.5);
            worstG = 0;
            for k = 1:height(glob)
                x = double(cer.(glob.variable(k))(:, :, glob.month(k)));
                worstG = max(worstG, abs(sum(x .* w, 'all') / (360 * sum(w)) / glob.value(k) - 1));
            end
            tc.verifyBoundAndRecord("importCERES global means vs the oracle (6 means)", worstG, 1e-13, "-", "64800 products summed in another order; rounding only");
        end
        function igrfCoefficientsReference(tc)
            % THE ORACLE IS THE PUBLISHED FILE, read a second way. These four numbers are
            % transcribed from igrf14coeffs.txt by eye; if a column or a degree/order index
            % moved, they are the first thing that fails.
            c = NonGravForces.io.igrfCoefficients(tc.f.igrf);
            i2020 = find(c.epochs == 2020); i2025 = find(c.epochs == 2025);
            tc.verifyBoundAndRecord("g_1^0 at 2020.0 against the published table", ...
                abs(c.gTable(2, 1, i2020) - (-29403.41)), 0, "nT", "the axial dipole, read off the file");
            tc.verifyBoundAndRecord("g_1^1 at 2020.0 against the published table", ...
                abs(c.gTable(2, 2, i2020) - (-1451.37)), 0, "nT", "the equatorial dipole");
            tc.verifyBoundAndRecord("h_1^1 at 2020.0 against the published table", ...
                abs(c.hTable(2, 2, i2020) - 4653.35), 0, "nT", "h is a separate table and must not be confused with g");
            tc.verifyBoundAndRecord("g_1^0 at 2025.0 against the published table", ...
                abs(c.gTable(2, 1, i2025) - (-29350.0)), 0, "nT", "the most recent definitive epoch");
            tc.verifyBoundAndRecord("secular variation of g_1^0", ...
                abs(c.svG(2, 1) - 12.6), 0, "nT/yr", "positive: the axial dipole is weakening, since g_1^0 is negative");
            % And the physics the numbers describe: the dipole dominates by an order of magnitude.
            dipole = sqrt(c.gTable(2,1,i2025)^2 + c.gTable(2,2,i2025)^2 + c.hTable(2,2,i2025)^2);
            quad = sqrt(sum(c.gTable(3,:,i2025).^2 + c.hTable(3,:,i2025).^2));
            tc.verifyBoundAndRecord("quadrupole over dipole at 2025.0", quad / dipole, 0.3, "-", ...
                "the Earth's field is dipole-dominated; a ratio near one would mean the degrees were transposed");
        end
        function readMAG1BReference(tc)
            % THE ORACLE IS THE FILE, read independently: the first data line is parsed here by
            % hand and compared against what the reader produced, so a column shifted by one
            % cannot pass. The published header names MfvX/Y/Z_RAW as columns 5 to 7 in
            % microtesla and the six torquer currents as columns 8 to 13 in milliampere.
            txt = splitlines(string(fileread(tc.f.mag)));
            first = find(startsWith(txt, "# End of YAML header"), 1) + 1;
            v = double(split(strtrim(txt(first))));
            [~, fld] = NonGravForces.io.readMAG1B(tc.f.mag);
            tc.verifyBoundAndRecord("readMAG1B field vs the first data line, microtesla to tesla", ...
                max(abs(fld.B_srf_T(:, 1) - v(5:7) * 1e-6)), 0, "T", ...
                "one multiplication by 1e-6; anything but zero means a column moved");
            tc.verifyBoundAndRecord("readMAG1B branch currents vs the first data line, mA to A", ...
                max(abs(fld.torquerBranch_A(:, 1) - v(8:13) * 1e-3)), 0, "A", ...
                "one multiplication by 1e-3, six columns in the published order");
            % And the physics: the Earth's field at GRACE-FO altitude is tens of microtesla.
            tc.verifyBoundAndRecord("measured field magnitude over the excerpt", ...
                max(fld.magnitude_T) * 1e6, 100, "uT", ...
                "the Earth's field at 490 km reaches about 50 microtesla; a bound of 100 catches a misread column without asserting the geomagnetic model");
            tc.verifyGreaterThan(max(fld.magnitude_T) * 1e6, 10, ...
                'a field magnitude below 10 microtesla over five minutes would mean the columns are not the field');
        end
        function readACC1BReference(tc)
            % Oracle: the file text itself, parsed independently (regexp tokens of the first data record).
            [~, a, aux] = NonGravForces.io.readACC1B(tc.f.act);
            tok = tc.rawRecord(tc.f.act, 1);          % gps_time GRACEFO_id lin_x lin_y lin_z ang_x ang_y ang_z res_x res_y res_z qualflg
            tc.verifyEqual(a(:, 1), str2double(tok(3:5)));
            tc.verifyEqual(aux.angular_rad_s2(:, 1), str2double(tok(6:8)));
            tc.verifyEqual(aux.residual_m_s2(:, 1), str2double(tok(9:11)));
            tc.verifyEqual(aux.units.lin_accl_x, 'm/s2');   % JPL D-56935: SI in the header
        end
        function readSCA1BReference(tc)
            % Physics oracle: the body +Z axis (graceFO "Nadir" plate normal) rotated
            % into the GCRS must point to nadir (-r/|r| from GNV1B). Measured 2026-09-06:
            % 0.84 deg mean, 0.86 deg max when conjugated; 96 deg mean when NOT conjugated.
            % Bound 3 deg: the GRACE-FO pitch offset and attitude deadband are ~1 deg.
            [~, q] = NonGravForces.io.readSCA1B(tc.f.sca);
            [~, r] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop);
            rhat = r ./ vecnorm(r);
            ez = squeeze(pagemtimes(referenceSystem.rotationFromQuaternion(q), [0; 0; 1]));
            ang = acosd(max(-1, min(1, sum(ez .* -rhat, 1))));
            tc.verifyBoundAndRecord("SCA1B body +Z vs nadir (conjugated)", max(ang), 3, "deg", "measured 0.86; 96 without the conjugation");
            % the same statistic without the conjugation must be far off, or this test proves nothing
            qNo = referenceSystem.invertQuaternion(q);
            ezNo = squeeze(pagemtimes(referenceSystem.rotationFromQuaternion(qNo), [0; 0; 1]));
            tc.verifyGreaterThan(mean(acosd(max(-1, min(1, sum(ezNo .* -rhat, 1))))), 30, 'the unconjugated sense must fail this test');
        end
        function readGNV1BReference(tc)
            % Oracle: Greenwich mean sidereal time from the IAU 1982 polynomial (Aoki et al.
            % 1982; Meeus 12.4), evaluated here without referenceSystem. The GCRS-ITRS
            % longitude difference of a point equals the hour angle of the equinox at
            % Greenwich up to precession-nutation since J2000 (~0.28 deg in 2020) and the
            % CIO/equinox difference. Measured 2026-09-06: -0.19 deg after the 18 s GPS-UTC
            % offset is applied; bound 0.5 deg. A wrong rotation sense gives ~160 deg.
            [~, r] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop);
            [~, rI] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Frame', "ITRS");
            T = (58849 + 2400000.5 - 2451545.0) / 36525;                       % centuries J2000 -> 2020-01-01 0h
            gmst0 = 24110.54841 + 8640184.812866 * T + 0.093104 * T^2 - 6.2e-6 * T^3;   % s, at 0h UT1
            ut1 = -18;                                                           % record 1: GPS 0 s = UTC 23:59:42 (dUT1 -0.18 s neglected: 0.0007 deg)
            gmst_deg = mod((gmst0 + 1.00273790935 * ut1) / 86400 * 360, 360);
            dlon = mod(rad2deg(atan2(r(2, 1), r(1, 1)) - atan2(rI(2, 1), rI(1, 1))), 360);
            tc.verifyBoundAndRecord("GNV1B GCRS-ITRS longitude vs IAU1982 GMST", abs(dlon - gmst_deg), 0.5, "deg", "precession since J2000 dominates (measured 0.19)");
        end
        function readTHR1BReference(tc)
            % Oracle: the raw first record (time_intg time_frac time_ref id 14 counts 14 on_times ...).
            [T, aux] = NonGravForces.io.readTHR1B(tc.f.thr);
            tok = tc.rawRecord(tc.f.thr, 1);
            tc.verifyEqual(aux.onTime_s(:, 1), str2double(tok(19:32)) / 1000);
            tc.verifyEqual(aux.count(:, 1), str2double(tok(5:18)));
            tc.verifyEqual(aux.thruster(1), "att_ctrl_1_1"); tc.verifyEqual(aux.thruster(13), "orb_ctrl_1");
            % GPS -> UTC: 18 leap seconds in 2020 (IERS Bulletin C 58); the tag is GPS seconds past 2000-01-01T12:00.
            t0 = datetime(2000, 1, 1, 12, 0, 0, 'TimeZone', 'UTC') + seconds(str2double(tok(1)) + str2double(tok(2)) * 1e-6) - seconds(18);
            tc.verifyLessThan(abs(seconds(T.tStart(1) - t0)), 1e-6, 'first event epoch vs raw GPS tag minus 18 s');
            tc.verifyEqual(height(T), 471, 'the complete day has 471 records, all with a firing');
        end
        function readMAS1BReference(tc)
            % Oracle: the values printed in the first record of the real file (also pinned by gracefoReading's realdata test).
            [~, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            tc.verifyEqual(tank.tank_kg(:, 1), [15.25970067419259; 15.23406422398089], 'AbsTol', 1e-12);
            tc.verifyLessThan(abs(tank.total_kg(1) - 30.5), 0.5, 'GRACE-FO carries ~30 kg of cold gas; not 601 kg');
        end
        function massModelReference(tc)
            % Oracle: JPL D-56935 masses (GRACE-C 601.214 kg, GRACE-D 601.209 kg).
            ep = referenceSystem.mjdToEpoch(58849.5, "GPS");
            tc.verifyEqual(NonGravForces.io.massModel("C").at(ep), 601.214);
            tc.verifyEqual(NonGravForces.io.massModel("D").at(ep), 601.209);
        end
        function spaceWeatherReference(tc)
            % Oracle: the raw CSV / TXT lines parsed here by hand (split on comma / whitespace).
            lines = splitlines(string(fileread(tc.f.sw)));
            row = @(d) split(lines(startsWith(lines, d)), ",");
            d0 = row("2020-01-01"); dm1 = row("2019-12-31"); dm2 = row("2019-12-30");
            sw = NonGravForces.io.spaceWeather(tc.f.sw);
            q = sw.query(referenceSystem.mjdToEpoch(58849.5, "UTC"));            % 2020-01-01 12:00 UTC, 3-h slot 5
            tc.verifyEqual(q.f107, str2double(dm1(25)), 'F10.7_OBS of the previous day');
            tc.verifyEqual(q.f107a, str2double(d0(28)), 'F10.7_OBS_CENTER81 of the day');
            tc.verifyEqual(q.ap, str2double(d0(21)), 'AP_AVG');
            tc.verifyEqual(q.kp, str2double(d0(8)) / 10, 'KP5');
            ap = @(r) str2double(r(13:20))';                                     % AP1..AP8
            hist = [str2double(d0(21)); ap(d0) * [0 0 0 0 1 0 0 0]'; ap(d0) * [0 0 0 1 0 0 0 0]'; ap(d0) * [0 0 1 0 0 0 0 0]'; ap(d0) * [0 1 0 0 0 0 0 0]'; ...
                    mean([ap(dm1) * [0 1 1 1 1 1 1 1]' / 7 * 7 / 8 * 8; 0]) * 0 + (sum(ap(dm1) .* [0 1 1 1 1 1 1 1]) + ap(d0) * [1 0 0 0 0 0 0 0]') / 8; ...
                    (sum(ap(dm2) .* [0 1 1 1 1 1 1 1]) + ap(dm1) * [1 0 0 0 0 0 0 0]') / 8];
            tc.verifyEqual(q.apHistory, hist, 'NRLMSISE-00 AP array: [daily, now, -3h, -6h, -9h, mean(-12..-33h), mean(-36..-57h)]');
            tc.verifyEqual(q.dataType, string(d0(27)));
            % JB2008: SOLFSMY line "2020 1": F10 S10 M10 Y10 with lags 1/1/2/5 d; DTC "DTC 2020 1" hourly.
            sl = splitlines(string(fileread(tc.f.solfsmy))); s1 = split(strtrim(sl(startsWith(strtrim(sl), "2020   1 "))));
            dl = splitlines(string(fileread(tc.f.dtc)));     dt = split(strtrim(dl(startsWith(dl, "DTC 2020   1 "))));
            sj = NonGravForces.io.spaceWeather([string(tc.f.solfsmy) string(tc.f.dtc)]);
            e = @(mjd) referenceSystem.mjdToEpoch(mjd, "UTC");
            tc.verifyEqual(sj.query(e(58850.5)).f10, str2double(s1(4)), 'F10 one day later');
            tc.verifyEqual(sj.query(e(58850.5)).s10b, str2double(s1(7)), 'S81c one day later');
            tc.verifyEqual(sj.query(e(58851.5)).m10, str2double(s1(8)), 'M10 two days later');
            tc.verifyEqual(sj.query(e(58854.5)).y10, str2double(s1(10)), 'Y10 five days later');
            tc.verifyEqual(sj.query(e(58849 + 0.5/24)).dtc, str2double(dt(4)), 'DTC hour 1');
            tc.verifyEqual(sj.query(e(58849 + 7.5/24)).dtc, str2double(dt(11)), 'DTC hour 8');
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function importCERESPrecision(tc)
            % The month nodes are exact: the file says 7228.25 and 7259 days after 2000-03-01 00:00 UTC, which are
            % 2019-12-15 06:00 and 2020-01-15 00:00, MJD 58832.25 and 58863 (the CF time units parsed, not assumed).
            cer = NonGravForces.io.importCERES(tc.ceres());
            tc.verifyEqual(cer.time, [datetime(2019, 12, 15, 6, 0, 0) datetime(2020, 1, 15)]);
            tc.verifyEqual(cer.mjd, [58832.25 58863]);
            % the permutation to lat x lon x time is exact against a direct read of the lon x lat x time variable
            raw = ncread(tc.f.rlut, 'rlut');
            tc.verifyEqual(cer.rlut, permute(raw, [2 1 3]), 'lat x lon x time, values untouched');
        end
        function igrfCoefficientsPrecision(tc)
            % Interpolation is one linear blend, so the two ends must be EXACT and the midpoint
            % must be the mean. Anything else means the weight or the bracketing index is wrong.
            c20 = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2020.0);
            c25 = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2025.0);
            cMid = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2022.5);
            tbl = NonGravForces.io.igrfCoefficients(tc.f.igrf);
            i20 = find(tbl.epochs == 2020); i25 = find(tbl.epochs == 2025);
            tc.verifyBoundAndRecord("interpolating to a published epoch returns that epoch", ...
                max(abs(c20.g - tbl.gTable(:, :, i20)), [], 'all'), 0, "nT", ...
                "a weight of zero must reproduce the node, not approximate it");
            tc.verifyBoundAndRecord("the midpoint is the mean of the two nodes", ...
                max(abs(cMid.g - 0.5 * (tbl.gTable(:, :, i20) + tbl.gTable(:, :, i25))), [], 'all'), ...
                1e-11, "nT", "one linear blend at w = 0.5; the bound is the rounding of numbers of order 3e4");
            % Past the last epoch the secular variation extrapolates, and it does so linearly.
            c27 = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2027.0);
            tc.verifyBoundAndRecord("extrapolation is the last epoch plus SV times dt", ...
                max(abs(c27.g - (c25.g + 2 * tbl.svG)), [], 'all'), 1e-11, "nT", ...
                "two years of secular variation, applied once");
        end
        function readMAG1BPrecision(tc)
            % The net current is the sum of two branches, and the sum is only meaningful because
            % the branches never overlap. That is asserted here rather than trusted: it was
            % measured over the full day and it must hold on whatever excerpt is present.
            [~, fld] = NonGravForces.io.readMAG1B(tc.f.mag);
            A = fld.torquerBranch_A(1:3, :); B = fld.torquerBranch_A(4:6, :);
            tc.verifyEqual(nnz(A ~= 0 & B ~= 0), 0, ...
                'an A and a B branch non-zero on the same axis would make the net sum a modelling choice, not an addition');
            tc.verifyBoundAndRecord("torquer net equals A plus B exactly", ...
                max(abs(fld.torquerNet_A - (A + B)), [], 'all'), 0, "A", ...
                "one addition, quoted back");
            tc.verifyTrue(all(A(:) >= 0), 'the A branches carry the positive direction and are never negative');
            tc.verifyTrue(all(B(:) <= 0), 'the B branches carry the negative direction and are never positive');
        end
        function readACC1BPrecision(tc)
            % Time-tag round trip is exact: epoch -> [day sec] equals gracefoReading's time bitwise (bound 0).
            [ep] = NonGravForces.io.readACC1B(tc.f.act);
            t = gracefoReading.readL1b(tc.f.act);
            back = referenceSystem.epochTo(ep, "daySecond");
            tc.verifyBoundAndRecord("ACC1B time tag round trip", max(abs(back(:) - t(:))), 0, "s", "exact by construction");
        end
        function readSCA1BPrecision(tc)
            % invertQuaternion normalises: |norm(q) - 1| at the eps level (measured 2.2e-16); bound 1e-15 = 4.5 eps.
            [~, q] = NonGravForces.io.readSCA1B(tc.f.sca);
            tc.verifyBoundAndRecord("SCA1B unit norm of q_body2gcrs", max(abs(vecnorm(q) - 1)), 1e-15, "-", "measured 2.2e-16");
        end
        function readGNV1BPrecision(tc)
            % The rotated velocity must contain the frame-rate term: a five-point central
            % difference of the GCRS positions (1 Hz) reproduces it to 4.4e-5 m/s (measured;
            % dominated by the GNV1B position noise ~mm and the stencil truncation at 1 s).
            % Bound 1e-3 m/s: 20x the measurement, 5 orders below the 470 m/s error of a
            % dropped Rdot term.
            [~, r, v] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop);
            h = 1;
            vfd = (-r(:, 5:end) + 8 * r(:, 4:end-1) - 8 * r(:, 2:end-3) + r(:, 1:end-4)) / (12 * h);
            tc.verifyBoundAndRecord("GNV1B velocity vs 5-point FD of GCRS position", max(vecnorm(vfd - v(:, 3:end-2))), 1e-3, "m/s", "measured 4.4e-5; 4.7e2 without Rdot");
        end
        function readTHR1BPrecision(tc)
            % Event duration equals the longest on-time of the record exactly (integer ms / 1000, one division).
            [T, aux] = NonGravForces.io.readTHR1B(tc.f.thr);
            dur = seconds(T.tEnd - T.tStart)';
            expect = max(aux.onTime_s(:, any(aux.fired, 1)), [], 1);
            tc.verifyBoundAndRecord("THR1B event duration vs max on-time", max(abs(dur - expect)), 1e-6, "s", "datetime resolution ~1e-7 s near 2020");
        end
        function readMAS1BPrecision(tc)
            % total is the plain sum tank1 + tank2 (same operands, same order): bitwise.
            [~, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            tc.verifyBoundAndRecord("MAS1B total vs tank1+tank2", max(abs(tank.total_kg - sum(tank.tank_kg, 1))), 0, "kg", "bitwise");
        end
        function massModelPrecision(tc)
            % At the MAS1B sample epochs the model equals anchor + gas - gas(1): interp1 at a node is
            % exact, the two additions round at 1e-13 kg (600 kg * eps); bound 1e-9 kg leaves 4 decades.
            [epM, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            m = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg);
            expect = 601.214 + tank.total_kg - tank.total_kg(1);
            tc.verifyBoundAndRecord("massModel mas1b at the sample epochs", max(abs(m.at(epM) - expect)), 1e-9, "kg", "rounding 1e-13");
            tc.verifyTrue(all(m.at(epM) > 600.9 & m.at(epM) < 601.3));
        end
        function spaceWeatherPrecision(tc)
            % kp is the file integer / 10 (the only arithmetic): 10*kp is an integer to 4 eps.
            sw = NonGravForces.io.spaceWeather(tc.f.sw);
            ep = referenceSystem.daySecondToEpoch([repmat(58849, 8, 1), (0:7)' * 10800 + 1], "UTC");
            d = sw.query(ep);
            tc.verifyBoundAndRecord("spaceWeather 10*kp integer", max(abs(10 * d.kp - round(10 * d.kp))), 1e-14, "-", "one division");
            tc.verifyEqual(d.kp, sw.index.kp3h(sw.mjd == 58849, :), 'the eight 3-h slots of the day in order');
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function importCERESVectorization(tc)
            % A period that sits exactly on one node loads that month alone; each such single-month import is the
            % corresponding page of the whole import, bitwise.
            f = tc.ceres();
            cer = NonGravForces.io.importCERES(f);
            for m = 1:2
                one = NonGravForces.io.importCERES(f, 'Period', [cer.time(m) cer.time(m)]);
                tc.verifyEqual(numel(one.time), 1);
                tc.verifyEqual(one.rsut, cer.rsut(:, :, m)); tc.verifyEqual(one.rlut, cer.rlut(:, :, m)); tc.verifyEqual(one.rsdt, cer.rsdt(:, :, m));
            end
        end
        function igrfCoefficientsVectorization(tc)
            % The tables are read once and indexed; asking for many epochs must give the same
            % numbers as asking for each, because there is no state between calls.
            years = [1905.0, 1963.7, 2000.0, 2019.25, 2028.0];
            tbl = NonGravForces.io.igrfCoefficients(tc.f.igrf);
            worst = 0;
            for y = years
                c = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', y);
                tc.verifyEqual(c.year, y);
                tc.verifyEqual(c.epochs, tbl.epochs, 'the tables must not depend on the epoch asked for');
                worst = max(worst, max(abs(c.gTable(:) - tbl.gTable(:))));
            end
            tc.verifyBoundAndRecord("the tables are identical whatever epoch is requested", worst, 0, "nT", ...
                "interpolation reads the tables; it may not modify them");
        end
        function readMAG1BVectorization(tc)
            % A TimeRange is a filter, not a different reader: the selected records must be the
            % same numbers the full read produced, at the same epochs.
            [ep, fld] = NonGravForces.io.readMAG1B(tc.f.mag);
            mjd = referenceSystem.epochTo(ep, "mjd");
            cut = [mjd(1) - 1e-9, mjd(200) + 1e-9];
            [ep2, fld2] = NonGravForces.io.readMAG1B(tc.f.mag, 'TimeRange', cut);
            n = numel(ep2.sec);
            tc.verifyLessThanOrEqual(n, numel(ep.sec));
            tc.verifyGreaterThan(n, 0);
            tc.verifyBoundAndRecord("a TimeRange selects records without changing them", ...
                max(abs(fld2.B_srf_T - fld.B_srf_T(:, 1:n)), [], 'all'), 0, "T", ...
                "a filter that alters a value is not a filter");
        end
        function readACC1BVectorization(tc)
            % A TimeRange excerpt equals the slice of the full read, bitwise (leading and interior).
            [ep, a, aux] = NonGravForces.io.readACC1B(tc.f.act);
            [e1, a1, x1] = NonGravForces.io.readACC1B(tc.f.act, 'TimeRange', tc.range(0, 9.5));
            tc.verifyEqual(a1, a(:, 1:10)); tc.verifyEqual(e1.sec, ep.sec(1:10)); tc.verifyEqual(x1.angular_rad_s2, aux.angular_rad_s2(:, 1:10));
            [e2, a2] = NonGravForces.io.readACC1B(tc.f.act, 'TimeRange', tc.range(300, 309.5));
            tc.verifyEqual(a2, a(:, 301:310)); tc.verifyEqual(e2.sec, ep.sec(301:310));
        end
        function readSCA1BVectorization(tc)
            [~, q] = NonGravForces.io.readSCA1B(tc.f.sca);
            [~, q2] = NonGravForces.io.readSCA1B(tc.f.sca, 'TimeRange', tc.range(300, 309.5));
            tc.verifyEqual(q2, q(:, 301:310));
        end
        function readGNV1BVectorization(tc)
            % The rotation of a batch must equal the rotation of a subset, bitwise (measured 0 difference).
            [~, r, v] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop);
            [~, r1, v1] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop, 'TimeRange', tc.range(0, 9.5));
            [~, r2, v2] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop, 'TimeRange', tc.range(300, 309.5));
            tc.verifyEqual(r1, r(:, 1:10)); tc.verifyEqual(v1, v(:, 1:10));
            tc.verifyEqual(r2, r(:, 301:310)); tc.verifyEqual(v2, v(:, 301:310));
        end
        function readTHR1BVectorization(tc)
            % Whole day equals the two half days concatenated (isequaln: amplitude and sigma are NaN by design).
            T = NonGravForces.io.readTHR1B(tc.f.thr);
            T1 = NonGravForces.io.readTHR1B(tc.f.thr, 'TimeRange', [58849, 58849.5]);
            T2 = NonGravForces.io.readTHR1B(tc.f.thr, 'TimeRange', [58849.5 + 1e-9, 58850]);
            tc.verifyTrue(isequaln([T1; T2], T)); tc.verifyGreaterThan(height(T1), 0); tc.verifyGreaterThan(height(T2), 0);
        end
        function readMAS1BVectorization(tc)
            [ep, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            [e2, t2] = NonGravForces.io.readMAS1B(tc.f.mas, 'TimeRange', [58849 + 5/24, 58849 + 10/24]);
            k = ep.sec >= 5 * 3600 & ep.sec <= 10 * 3600;
            tc.verifyEqual(e2.sec, ep.sec(k)); tc.verifyEqual(t2.total_kg, tank.total_kg(k')); tc.verifyEqual(t2.tank_kg, tank.tank_kg(:, k'));
        end
        function massModelVectorization(tc)
            [epM, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            m = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg);
            ep = referenceSystem.daySecondToEpoch([repmat(58849, 7, 1), (1000:5000:31000)'], "GPS");
            all7 = m.at(ep);
            one = arrayfun(@(k) m.at(referenceSystem.daySecondToEpoch([58849, ep.sec(k)], "GPS")), 1:7);
            tc.verifyEqual(all7, one);
        end
        function spaceWeatherVectorization(tc)
            sw = NonGravForces.io.spaceWeather(tc.f.sw); sj = NonGravForces.io.spaceWeather([string(tc.f.solfsmy) string(tc.f.dtc)]);
            mjd = [58849.01; 58849.5; 58849.99; 58860.3; 58900.7];
            for s = {sw, sj}
                dAll = s{1}.query(referenceSystem.mjdToEpoch(mjd, "UTC"));
                for k = 1:numel(mjd)
                    d1 = s{1}.query(referenceSystem.mjdToEpoch(mjd(k), "UTC"));
                    for fn = string(fieldnames(dAll))'
                        tc.verifyEqual(dAll.(fn)(:, k), d1.(fn), sprintf('%s field %s epoch %d', s{1}.format, fn, k));
                    end
                end
            end
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function importCERESRobustness(tc)
            f = tc.ceres();
            cer = NonGravForces.io.importCERES(f);
            tc.verifyEqual(NonGravForces.io.importCERES(f([3 1 2])), cer, 'the files are recognised by variable_id, not by position');
            % a fill value in the loaded months is refused and named, not propagated as 1e20 W/m^2
            holed = fullfile(tc.tmp, 'holed_rsut.nc'); copyfile(tc.f.rsut, holed, 'f'); fileattrib(holed, '+w');
            ncwrite(holed, 'rsut', single(1e20), [10 20 2]);
            tc.verifyErrorId(@() NonGravForces.io.importCERES([string(holed) f(2:3)]), 'NonGravForces:importCERES:dataMissing');
            % a units attribute that is not W m-2 is refused
            milli = fullfile(tc.tmp, 'milli_rlut.nc'); copyfile(tc.f.rlut, milli, 'f'); fileattrib(milli, '+w');
            ncwriteatt(milli, 'rlut', 'units', 'mW m-2');
            tc.verifyErrorId(@() NonGravForces.io.importCERES([f(1) string(milli) f(3)]), 'NonGravForces:importCERES:dataFormat');
            % time units with a time of day are parsed, not truncated to midnight (MATLAB's regexp dropped those tokens)
            % ... and a file whose time axis disagrees with the other two is refused
            shifted = string(fullfile(tc.tmp, ["a_rsut.nc" "b_rlut.nc" "c_rsdt.nc"]));
            for k = 1:3, copyfile(f(k), shifted(k), 'f'); fileattrib(shifted(k), '+w'); end
            ncwriteatt(shifted(3), 'time', 'units', 'days since 2000-03-01 12:30:00');
            tc.verifyErrorId(@() NonGravForces.io.importCERES([f(1:2) shifted(3)]), 'NonGravForces:importCERES:dataFormat');
            for k = 1:2, ncwriteatt(shifted(k), 'time', 'units', 'days since 2000-03-01 12:30:00'); end
            c2 = NonGravForces.io.importCERES(shifted);
            tc.verifyEqual(c2.time - cer.time, repmat(hours(12) + minutes(30), 1, 2), 'the time of day of the origin is honoured');
        end
        function igrfCoefficientsRobustness(tc)
            % A truncated row must be refused with the mnemonic that names the cause.
            bad = fullfile(tc.tmp, 'igrf_truncated.txt');
            L = splitlines(string(fileread(tc.f.igrf)));
            k = find(startsWith(strtrim(L), "g  1  0"), 1);
            p = split(strtrim(L(k)));
            L(k) = strjoin(p(1:end-3), " ");
            writelines(L, bad);
            tc.verifyErrorId(@() NonGravForces.io.igrfCoefficients(bad), 'NonGravForces:igrfCoefficients:dataFormat');
            % An impossible degree/order pair is refused too, and not silently placed somewhere.
            bad2 = fullfile(tc.tmp, 'igrf_badorder.txt');
            L2 = splitlines(string(fileread(tc.f.igrf)));
            p2 = split(strtrim(L2(k)));
            p2(3) = "5";                                   % m = 5 with n = 1
            L2(k) = strjoin(p2, " ");
            writelines(L2, bad2);
            tc.verifyErrorId(@() NonGravForces.io.igrfCoefficients(bad2), 'NonGravForces:igrfCoefficients:dataColumn');
            % And the unaltered file must pass, or neither refusal proves anything (rule T6).
            tc.verifyWarningFree(@() NonGravForces.io.igrfCoefficients(tc.f.igrf));
            % The upper triangle and the zeroth row are zero and stay in the array.
            c = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2020.0);
            tc.verifyEqual(c.g(1, :), zeros(1, 14), 'degree zero has no coefficient');
            tc.verifyEqual(c.h(:, 1), zeros(14, 1), 'h has no m = 0 term');
            tc.verifyEqual(c.g(2, 3:end), zeros(1, 12), 'm may not exceed n');
        end
        function readMAG1BRobustness(tc)
            % A corrupted copy must be refused with the mnemonic that names the cause, not with
            % whatever MATLAB happens to throw first.
            bad = fullfile(tc.tmp, 'MAG1B_2020-01-01_C_04.txt');
            txt = string(fileread(tc.f.mag));
            % A field of 900 microtesla: physically impossible at any altitude, and exactly what
            % reading the wrong column would produce.
            lines = splitlines(txt);
            k = find(startsWith(lines, "# End of YAML header"), 1) + 1;
            parts = split(strtrim(lines(k)));
            parts(5) = "900.0";
            lines(k) = strjoin(parts, " ");
            writelines(lines, bad);
            tc.verifyErrorId(@() NonGravForces.io.readMAG1B(bad), 'NonGravForces:readMAG1B:dataRange');
            % The same file unaltered must pass, or the refusal above proves nothing (rule T6).
            tc.verifyWarningFree(@() NonGravForces.io.readMAG1B(tc.f.mag));
            % Quality bits are reported even when every one of them is clear, which is the case
            % on this day: an empty report and an all-clear report are different statements.
            [~, ~, aux] = NonGravForces.io.readMAG1B(tc.f.mag);
            tc.verifySize(aux.qualityBits, [8 600]);
        end
        function readACC1BRobustness(tc)
            % A quality mask naming a bit no record sets removes nothing; two files of different products refuse (foreign id).
            [~, a] = NonGravForces.io.readACC1B(tc.f.act);
            [~, a2, aux2] = NonGravForces.io.readACC1B(tc.f.act, 'QualityBitMask', 7);
            tc.verifyEqual(a2, a); tc.verifyEqual(aux2.numRecords, 600);
            tc.verifyError(@() NonGravForces.io.readACC1B([string(tc.f.act) string(tc.f.sca)]), 'gracefoReading:readL1bMulti:inconsistentVariables');
            % a header missing lin_accl_z: refuse by name, never by position
            bad = tc.corrupt(tc.f.act, 'ACT1B_nocol.txt', @(L, k) regexprep(L, '^(\s*- )lin_accl_z(:)', '$1lin_accl_z_renamed$2'));
            tc.verifyErrorId(@() NonGravForces.io.readACC1B(bad), 'NonGravForces:readACC1B:dataColumn');
        end
        function readACC1BMultiEmptyRobustness(tc)
            % O-12 discharge. Until gracefoReading 1.0.1, l1bLoad sent one file
            % through readL1b and several through readL1bMulti, because the
            % multi-file reader raised MATLAB:badsubscript on an all-empty
            % TimeRange. dataEmpty was therefore guaranteed for one file and
            % UNREACHABLE for several. The split is gone; this asserts the path
            % that could not be tested before.
            % The same fixture twice is enough: the range excludes every record,
            % so the concatenation is empty and no ordering check is reached.
            two = [string(tc.f.act) string(tc.f.act)];
            tc.verifyErrorId(@() NonGravForces.io.readACC1B(two, 'TimeRange', [58850 58851]), ...
                'NonGravForces:readACC1B:dataEmpty');
            % and one file still refuses identically, through the same path
            tc.verifyErrorId(@() NonGravForces.io.readACC1B(tc.f.act, 'TimeRange', [58850 58851]), ...
                'NonGravForces:readACC1B:dataEmpty');
        end
        function readSCA1BRobustness(tc)
            % A quaternion off the unit sphere is refused, not repaired.
            bad = tc.corrupt(tc.f.sca, 'SCA1B_bad.txt', @(L, k) localScaleToken(L, k, 4, 2));
            tc.verifyErrorId(@() NonGravForces.io.readSCA1B(bad), 'NonGravForces:readSCA1B:dataQuaternion');
            tc.verifyErrorId(@() NonGravForces.io.readSCA1B(tc.f.sca, 'TimeRange', [58850 58851]), 'NonGravForces:readSCA1B:dataEmpty');
        end
        function readGNV1BRobustness(tc)
            % coord_ref other than 'E' is refused; the whole 1962-now C04 series propagates referenceSystem's refusal unwrapped.
            bad = tc.corrupt(tc.f.gnv, 'GNV1B_bad.txt', @(L, k) localSetToken(L, k, 3, "I"));
            tc.verifyErrorId(@() NonGravForces.io.readGNV1B(bad, 'Eop', tc.eop), 'NonGravForces:readGNV1B:dataFrame');
            tc.verifyErrorId(@() NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop, 'TimeRange', [58850 58851]), 'NonGravForces:readGNV1B:dataEmpty');
            [~, r] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Frame', "ITRS", 'Eop', tc.eop);   % Eop ignored, not refused, for ITRS
            tc.verifySize(r, [3 600]);
        end
        function readTHR1BRobustness(tc)
            bad = tc.corrupt(tc.f.thr, 'THR1B_bad.txt', @(L, k) localSetToken(L, k, 3, "U"));
            tc.verifyErrorId(@() NonGravForces.io.readTHR1B(bad), 'NonGravForces:readTHR1B:dataTimeRef');
            % a record with no firing is counted, not listed
            noFire = tc.corrupt(tc.f.thr, 'THR1B_nofire.txt', @(L, k) localZeroTokens(L, k, 19:32));
            [T, aux] = NonGravForces.io.readTHR1B(noFire);
            tc.verifyEqual(aux.numRecordsNoFiring, 1); tc.verifyEqual(height(T), aux.numRecords - 1);
        end
        function readMAS1BRobustness(tc)
            % 150 kg in a tank column is not propellant: refused.
            bad = tc.corrupt(tc.f.mas, 'MAS1B_bad.txt', @(L, k) localSetToken(L, k, 7, "150.0"));
            tc.verifyErrorId(@() NonGravForces.io.readMAS1B(bad), 'NonGravForces:readMAS1B:dataRange');
            tc.verifyErrorId(@() NonGravForces.io.readMAS1B(tc.f.mas, 'TimeRange', [58851 58852]), 'NonGravForces:readMAS1B:dataEmpty');
        end
        function massModelRobustness(tc)
            [epM, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            m = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg);
            % hold within two sample intervals (2 h) beyond the ends; refuse further
            tc.verifyEqual(m.at(referenceSystem.mjdToEpoch(58849, "GPS")), 601.214, 'first hour of the day: held at the first sample');
            tc.verifyEqual(m.at(referenceSystem.mjdToEpoch(58850 + 1.9/24, "GPS")), m.at(epM) * [zeros(23, 1); 1], 'AbsTol', 1e-12);
            tc.verifyErrorId(@() m.at(referenceSystem.mjdToEpoch(58850 + 2.1/24, "GPS")), 'NonGravForces:massModel:inputEpoch');
            tc.verifyErrorId(@() m.at(referenceSystem.mjdToEpoch(58848.5, "GPS")), 'NonGravForces:massModel:inputEpoch');
            g = tank.total_kg; g(5) = NaN;
            tc.verifyErrorId(@() NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', g), 'NonGravForces:massModel:inputMas1b');
            epRev = epM; epRev.sec = flipud(epM.sec);
            tc.verifyErrorId(@() NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epRev, 'GasMass_kg', tank.total_kg), 'NonGravForces:massModel:inputMas1b');
            % a gas series that would push the mass below 500 kg is refused at evaluation
            bigLoss = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg - 200 * (0:23) / 23);
            tc.verifyErrorId(@() bigLoss.at(epM), 'NonGravForces:massModel:numericRange');
        end
        function spaceWeatherRobustness(tc)
            sw = NonGravForces.io.spaceWeather(tc.f.sw);
            e = @(mjd) referenceSystem.mjdToEpoch(mjd, "UTC");
            tc.verifyErrorId(@() sw.query(e(58720)), 'NonGravForces:spaceWeather:dataRange');            % before the table
            tc.verifyErrorId(@() sw.query(e(58727.5)), 'NonGravForces:spaceWeather:dataRange');          % first day: no previous day for F10.7
            tc.verifyErrorId(@() sw.query(e(58940)), 'NonGravForces:spaceWeather:dataRange');            % after the table
            d = sw.query(e(58730.5)); tc.verifyTrue(isfinite(d.f107), 'fourth day: 57 h history available');
            % a blank Kp field in the day queried is refused, never returned as NaN
            lines = splitlines(string(fileread(tc.f.sw))); k = find(startsWith(lines, "2020-01-01")); tok = split(lines(k), ","); tok(8) = ""; lines(k) = strjoin(tok, ",");
            blank = fullfile(tc.tmp, 'SW_blank.csv'); fid = fopen(blank, 'w'); fprintf(fid, '%s\n', lines); fclose(fid);
            swb = NonGravForces.io.spaceWeather(blank);
            tc.verifyErrorId(@() swb.query(e(58849.5)), 'NonGravForces:spaceWeather:dataMissing');     % slot 5 of 2020-01-01
            tc.verifyTrue(isfinite(swb.query(e(58849.1)).kp), 'slot 1 of the same day is intact');
            junk = fullfile(tc.tmp, 'junk.txt'); fid = fopen(junk, 'w'); fprintf(fid, 'hello\nworld\n'); fclose(fid);
            tc.verifyErrorId(@() NonGravForces.io.spaceWeather(junk), 'NonGravForces:spaceWeather:inputFormat');
            sl = splitlines(string(fileread(tc.f.solfsmy))); sl(10) = "  2019 249   garbage";
            badS = fullfile(tc.tmp, 'SOLFSMY_bad.TXT'); fid = fopen(badS, 'w'); fprintf(fid, '%s\n', sl); fclose(fid);
            tc.verifyErrorId(@() NonGravForces.io.spaceWeather([string(badS) string(tc.f.dtc)]), 'NonGravForces:spaceWeather:dataMalformed');
            sj = NonGravForces.io.spaceWeather([string(tc.f.solfsmy) string(tc.f.dtc)]);
            tc.verifyErrorId(@() sj.query(e(58731.5)), 'NonGravForces:spaceWeather:dataRange');          % fifth day: Y10 lag needs day-5
            tc.verifyTrue(isfinite(sj.query(e(58732.5)).y10));
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function importCERESMetamorphic(tc)
            % Any period strictly inside the two nodes needs both and returns the whole import unchanged; the
            % three fields are independent, so the rsut of one import is the rsut of any other.
            f = tc.ceres();
            cer = NonGravForces.io.importCERES(f);
            mid = NonGravForces.io.importCERES(f, 'Period', [datetime(2020, 1, 1) datetime(2020, 1, 2, 12, 0, 0)]);
            tc.verifyEqual(mid, cer);
            zoned = NonGravForces.io.importCERES(f, 'Period', datetime([2020 2020], [1 1], [1 2], 'TimeZone', 'Europe/Berlin'));
            tc.verifyEqual(zoned.time, cer.time, 'a zoned period is compared in UTC');
        end
        function igrfCoefficientsMetamorphic(tc)
            tbl = NonGravForces.io.igrfCoefficients(tc.f.igrf);
            % 1. Interpolation is monotone in the weight: halfway is between the ends, coefficient
            %    by coefficient. A sign error in the weight would break this and nothing else.
            i20 = find(tbl.epochs == 2020); i25 = find(tbl.epochs == 2025);
            a = tbl.gTable(:, :, i20); b = tbl.gTable(:, :, i25);
            mid = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2022.5).g;
            between = (mid >= min(a, b) - 1e-9) & (mid <= max(a, b) + 1e-9);
            tc.verifyTrue(all(between, 'all'), 'every interpolated coefficient must lie between its two nodes');
            % 2. A quarter of the way is a quarter of the difference.
            q = NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2021.25).g;
            tc.verifyBoundAndRecord("a quarter step moves a quarter of the difference", ...
                max(abs((q - a) - 0.25 * (b - a)), [], 'all'), 1e-11, "nT", ...
                "linearity in the weight, checked away from the midpoint where a factor of two would hide");
            % 3. Secular variation exists only where the model publishes it: degree <= 8.
            tc.verifyEqual(nnz(tbl.svG(10:end, :)), 0, 'SV is published to degree 8 and no further');
            tc.verifyGreaterThan(nnz(tbl.svG(2:9, :)), 0, 'and it must be present below that');
        end
        function readMAG1BMetamorphic(tc)
            [~, fld, aux] = NonGravForces.io.readMAG1B(tc.f.mag);
            % 1. The magnitude is the norm of the components, by definition and not by accident.
            tc.verifyBoundAndRecord("magnitude equals the norm of the components", ...
                max(abs(fld.magnitude_T - vecnorm(fld.B_srf_T))), 0, "T", "one vecnorm, quoted back");
            % 2. Unit conversion is a scaling and nothing else: the ratio to the file's numbers is
            %    exactly 1e-6 everywhere, so no offset or clipping crept in.
            txt = splitlines(string(fileread(tc.f.mag)));
            k = find(startsWith(txt, "# End of YAML header"), 1);
            raw = zeros(3, 600);
            for r = 1:600
                v = double(split(strtrim(txt(k + r))));
                raw(:, r) = v(5:7);
            end
            nz = raw ~= 0;
            ratios = fld.B_srf_T(nz) ./ raw(nz);
            tc.verifyBoundAndRecord("microtesla to tesla is a pure scaling", ...
                max(abs(ratios - 1e-6)) / 1e-6, 1e-15, "-", ...
                "every non-zero component divided by its file value gives the same factor");
            % 3. The calibration columns are returned unapplied. On this day they are constant,
            %    which is recorded rather than asserted as a law: a later file may vary them.
            tc.verifyBoundAndRecord("calibration columns are constant over this excerpt", ...
                max(max(aux.calibration.magField, [], 2) - min(aux.calibration.magField, [], 2)), 0, "-", ...
                "measured [0, 65000, -65000] for all 172800 records of the full day; recorded because the values look like sentinels, not scale factors");
        end
        function readACC1BMetamorphic(tc)
            % Split/merge: two adjoining TimeRanges concatenated equal the whole read, bitwise.
            [ep, a, aux] = NonGravForces.io.readACC1B(tc.f.act);
            [e1, a1, x1] = NonGravForces.io.readACC1B(tc.f.act, 'TimeRange', tc.range(0, 299.5));
            [e2, a2, x2] = NonGravForces.io.readACC1B(tc.f.act, 'TimeRange', tc.range(299.5, 600));
            tc.verifyEqual([a1 a2], a); tc.verifyEqual([e1.sec; e2.sec], ep.sec);
            tc.verifyEqual([x1.qualflgBits; x2.qualflgBits], aux.qualflgBits);
        end
        function readSCA1BMetamorphic(tc)
            % Inverting the output once more reproduces the (normalised) file quaternion at the eps level
            % (two normalisations; measured 3.3e-16). Bound 1e-14.
            [~, q] = NonGravForces.io.readSCA1B(tc.f.sca);
            [~, obs] = gracefoReading.readL1b(tc.f.sca);
            qFile = referenceSystem.normalizeQuaternion(obs(:, 1:4)');
            tc.verifyBoundAndRecord("SCA1B double conjugation vs file", max(abs(referenceSystem.invertQuaternion(q) - qFile), [], 'all'), 1e-14, "-", "measured 3.3e-16");
        end
        function readGNV1BMetamorphic(tc)
            % (1) a rotation preserves |r| (measured 4e-16 relative; bound 1e-14);
            % (2) round trip GCRS -> ITRS through applyRotation's inverse reproduces the file
            %     (measured 5e-9 m; bound 1e-6 m, 1000x the double rounding at 7e6 m and 1000x below
            %     the HANDOVER requirement of 1 mm);
            % (3) an Eop struct and the same series as a file give bitwise the same state.
            [ep, r, v] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop);
            [~, rI, vI] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Frame', "ITRS");
            tc.verifyBoundAndRecord("GNV1B |r| GCRS vs ITRS", max(abs(vecnorm(r) - vecnorm(rI)) ./ vecnorm(rI)), 1e-14, "rel", "measured 4e-16");
            [R, Rdot] = referenceSystem.rotationItrsToGcrs(ep, 'Eop', tc.eop);
            [rB, vB] = referenceSystem.applyRotation(R, Rdot, r, v, 'Inverse', true);
            tc.verifyBoundAndRecord("GNV1B round trip GCRS->ITRS position", max(vecnorm(rB - rI)), 1e-6, "m", "measured 5e-9; requirement 1e-3");
            tc.verifyBoundAndRecord("GNV1B round trip GCRS->ITRS velocity", max(vecnorm(vB - vI)), 1e-9, "m/s", "measured 1e-11");
            [~, rF, vF] = NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.f.eop);
            tc.verifyEqual(rF, r); tc.verifyEqual(vF, v);
        end
        function readTHR1BMetamorphic(tc)
            % The event table and aux describe the same records: one event per fired record, in order,
            % and the UTC epochs are the GPS epochs minus the leap-second count.
            [T, aux] = NonGravForces.io.readTHR1B(tc.f.thr);
            fired = any(aux.fired, 1);
            tc.verifyEqual(height(T), nnz(fired));
            gps = referenceSystem.epochTo(aux.epoch, "datetime"); gps.TimeZone = 'UTC';
            dt = seconds(gps(fired) - T.tStart);
            tc.verifyLessThan(max(abs(dt - 18)), 1e-6, 'GPS - UTC = 18 s throughout 2020');
            tc.verifyTrue(issorted(T.tStart));
        end
        function readMAS1BMetamorphic(tc)
            % Split/merge over TimeRange equals the whole read (isequaln: NaN columns by design).
            [ep, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            [e1, t1] = NonGravForces.io.readMAS1B(tc.f.mas, 'TimeRange', [58849, 58849.5]);
            [e2, t2] = NonGravForces.io.readMAS1B(tc.f.mas, 'TimeRange', [58849.5, 58850]);
            tc.verifyEqual([e1.sec; e2.sec], ep.sec);
            tc.verifyTrue(isequaln([t1.tank_kg t2.tank_kg], tank.tank_kg));
            tc.verifyTrue(isequaln([t1.total_kg t2.total_kg], tank.total_kg));
        end
        function massModelMetamorphic(tc)
            % (1) shifting the anchor by D shifts every mass by D (one addition: 1e-12 kg);
            % (2) the epoch scale of the query does not matter (converted inside): bitwise-equal day/sec
            %     after the 18 s GPS-UTC shift, so identical masses to 1e-12 kg.
            [epM, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            m1 = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg);
            m2 = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg, 'Anchor_kg', 601.214 + 10);
            ep = referenceSystem.daySecondToEpoch([repmat(58849, 5, 1), (3600:7200:32400)'], "GPS");
            tc.verifyBoundAndRecord("massModel anchor shift", max(abs(m2.at(ep) - m1.at(ep) - 10)), 1e-12, "kg", "one addition at 600 kg");
            epU = referenceSystem.convertTimeScale(ep, "UTC");
            tc.verifyBoundAndRecord("massModel query scale invariance", max(abs(m1.at(epU) - m1.at(ep))), 1e-12, "kg", "exact 18 s shift");
            c = NonGravForces.io.massModel("D");
            tc.verifyEqual(c.at(ep), c.at(epU));
        end
        function spaceWeatherMetamorphic(tc)
            % Daily indices are step functions: 00:00:00 and 23:59:59 of one day agree on every daily field;
            % the JB2008 lags are relations between query and table: query(d).f10 == table(d-1), etc.
            sw = NonGravForces.io.spaceWeather(tc.f.sw);
            e = @(mjd) referenceSystem.mjdToEpoch(mjd, "UTC");
            a = sw.query(e(58860)); b = sw.query(e(58860 + 86399/86400));
            tc.verifyEqual([a.f107 a.f107a a.ap], [b.f107 b.f107a b.ap]); tc.verifyEqual(a.dataType, b.dataType);
            tc.verifyEqual(a.apHistory(1), b.apHistory(1));
            sj = NonGravForces.io.spaceWeather([string(tc.f.solfsmy) string(tc.f.dtc)]);
            for d = [58800 58849 58900]
                q = sj.query(e(d + 0.3));
                tc.verifyEqual(q.f10, sj.index.f10(sj.mjd == d - 1)); tc.verifyEqual(q.s10b, sj.index.s10b(sj.mjd == d - 1));
                tc.verifyEqual(q.m10, sj.index.m10(sj.mjd == d - 2)); tc.verifyEqual(q.y10b, sj.index.y10b(sj.mjd == d - 5));
                tc.verifyEqual(q.dtc, sj.index.dtc(sj.mjd == d, 8), 'hour index floor(0.3*24)+1 = 8');
            end
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function importCERESSpeed(tc)
            % Everything importCERES adds to three ncread calls (attribute checks, axis comparisons, the permutation and
            % the fill screen) against those bare reads: measured see last_run; budget 3x the 2026-09-15 figure.
            f = tc.ceres();
            base = @() cellfun(@(p, v) ncread(p, v), cellstr(f), {'rsut', 'rlut', 'rsdt'}, 'UniformOutput', false);
            q = tc.speedRatio(@() NonGravForces.io.importCERES(f), base);
            tc.verifyBoundAndRecord("importCERES vs three bare ncread calls", q, 25, "-", "checks and one permutation on 1.2e5 floats per field; measured 8.3 on 2026-09-15, budget 3x");
        end
        function igrfCoefficientsSpeed(tc)
            NonGravForces.io.igrfCoefficients(tc.f.igrf);
            t = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.io.igrfCoefficients(tc.f.igrf, 'Epoch', 2020.0); t(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("igrfCoefficients, 195 rows x 27 columns", 1e3 * median(t), 2000, "ms", ...
                "a smoke ceiling on a 42 kB parse, not a claim about this host");
        end
        function readMAG1BSpeed(tc)
            % Reading is dominated by gracefoReading's parse; what this checks is that nothing
            % quadratic sits in the marshalling above it.
            NonGravForces.io.readMAG1B(tc.f.mag);
            t = zeros(1, 3);
            for k = 1:3
                t0 = tic; NonGravForces.io.readMAG1B(tc.f.mag); t(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("readMAG1B, 600 records", 1e3 * median(t), 5000, "ms", ...
                "a smoke ceiling on a 600-record parse, not a claim about this host");
        end
        function readACC1BSpeed(tc)
            % Reader over raw gracefoReading.readL1b: measured 0.98 (median of 5, warm). Bound 5 = 3x a
            % generous 1.5 expectation for the epoch struct, checks and provenance.
            q = tc.speedRatio(@() NonGravForces.io.readACC1B(tc.f.act), @() gracefoReading.readL1b(tc.f.act));
            tc.verifyBoundAndRecord("readACC1B / readL1b time ratio", q, 5, "-", "measured 0.98");
        end
        function readSCA1BSpeed(tc)
            % Measured 1.46 (norm check + conjugation). Bound 5.
            q = tc.speedRatio(@() NonGravForces.io.readSCA1B(tc.f.sca), @() gracefoReading.readL1b(tc.f.sca));
            tc.verifyBoundAndRecord("readSCA1B / readL1b time ratio", q, 5, "-", "measured 1.46");
        end
        function readGNV1BSpeed(tc)
            % GCRS read over ITRS read (the rotation with Rdot for 600 epochs): measured 1.16 warm, 6.3 cold.
            % Bound 5 after warm-up (3x the 1.5 expectation).
            q = tc.speedRatio(@() NonGravForces.io.readGNV1B(tc.f.gnv, 'Eop', tc.eop), @() NonGravForces.io.readGNV1B(tc.f.gnv, 'Frame', "ITRS"));
            tc.verifyBoundAndRecord("readGNV1B GCRS / ITRS time ratio", q, 5, "-", "measured 1.16 warm");
        end
        function readTHR1BSpeed(tc)
            % Measured 1.04 (event table build is cheap). Bound 5.
            q = tc.speedRatio(@() NonGravForces.io.readTHR1B(tc.f.thr), @() gracefoReading.readL1b(tc.f.thr));
            tc.verifyBoundAndRecord("readTHR1B / readL1b time ratio", q, 5, "-", "measured 1.04");
        end
        function readMAS1BSpeed(tc)
            % Measured 1.47. Bound 5.
            q = tc.speedRatio(@() NonGravForces.io.readMAS1B(tc.f.mas), @() gracefoReading.readL1b(tc.f.mas));
            tc.verifyBoundAndRecord("readMAS1B / readL1b time ratio", q, 5, "-", "measured 1.47");
        end
        function massModelSpeed(tc)
            % at() over one day at 1 Hz against a bare interp1 of the same size: measured 2.24
            % (validation, clamping, the range check). Bound 7 = 3x.
            [epM, tank] = NonGravForces.io.readMAS1B(tc.f.mas);
            m = NonGravForces.io.massModel("C", 'Mode', "mas1b", 'Epoch', epM, 'GasMass_kg', tank.total_kg);
            ep = referenceSystem.daySecondToEpoch([repmat(58849, 86400, 1), (0:86399)'], "GPS");
            tg = epM.day + epM.sec / 86400; tq = ep.day + ep.sec / 86400;
            q = tc.speedRatio(@() m.at(ep), @() interp1(tg, tank.total_kg, min(max(tq, tg(1)), tg(end))));
            tc.verifyBoundAndRecord("massModel.at(86400) / interp1 time ratio", q, 7, "-", "measured 2.24");
        end
        function spaceWeatherSpeed(tc)
            % query of one day at 1 Hz (86400 epochs, 0.05 s measured) against an ismember day lookup of
            % the same size: measured 65 (four day rows, 20 history picks, five NaN scans). Bound 200 = 3x.
            sw = NonGravForces.io.spaceWeather(tc.f.sw);
            ep = referenceSystem.daySecondToEpoch([repmat(58849, 86400, 1), (0:86399)'], "UTC");
            base = @() localBaseline(sw.mjd, ep);
            q = tc.speedRatio(@() sw.query(ep), base);
            tc.verifyBoundAndRecord("spaceWeather.query(86400) / ismember time ratio", q, 200, "-", "measured 65");
        end
    end

    % ============================================== io.importGroopsGrid (LP-E WP-E1)
    methods (Static)
        function f = writeGroopsGrid(folder, name, type, lonDeg, latDeg, values, mjd)
            %WRITEGROOPSGRID Write a GROOPS binary file in the decoded layout, for the tests.
            %   THE TESTS WRITE THEIR OWN FIXTURE because the real grids are 1.25 GB and live in
            %   DATAPOOL. A round trip through this writer exercises exactly what could be
            %   misdecoded: the length-prefixed header, the padding, the count-prefixed vectors,
            %   the time nodes and the column-major value block.
            f = fullfile(folder, name);
            fid = fopen(f, 'w', 'l');
            closer = onCleanup(@() fclose(fid));
            fwrite(fid, 'B', 'char');
            for str = [string(type), "20200123"]
                fwrite(fid, strlength(str), 'uint64');
                fwrite(fid, char(str), 'char');
            end
            fwrite(fid, repmat(' ', 1, 4 + 2 * (type == "griddedDataTimeSeries")), 'char');
            isSeries = type == "griddedDataTimeSeries";
            if isSeries
                fwrite(fid, [0 numel(mjd) 1 1], 'uint64');
            else
                fwrite(fid, 1, 'uint64');
            end
            fwrite(fid, [6378137.0 6356752.314140], 'double');
            nLon = numel(lonDeg); nLat = numel(latDeg);
            % EXACT band solid angle, not a cosine times a width: the reference test asserts the
            % sum is 4 pi, and a first-order area would make the READER look wrong for a fault
            % that was the writer's (measured 2.6e-2 relative, 2026-09-17).
            dPhi = deg2rad(180 / nLat); dLam = 2 * pi / nLon;
            phi = deg2rad(latDeg(:));
            area = repmat(dLam * (sin(phi + dPhi / 2) - sin(phi - dPhi / 2)), 1, nLon);
            for v = {deg2rad(lonDeg(:)), deg2rad(latDeg(:)), 6378137 + zeros(nLat, 1), reshape(area', [], 1)}
                fwrite(fid, numel(v{1}), 'uint64');
                fwrite(fid, v{1}, 'double');
            end
            fwrite(fid, 1, 'uint64'); fwrite(fid, 1, 'double');
            if isSeries
                for k = 1:numel(mjd)
                    fwrite(fid, mod(mjd(k), 1), 'double');
                    fwrite(fid, floor(mjd(k)), 'int64');
                end
                fwrite(fid, 0, 'uint64');
            else
                fwrite(fid, size(values, 3), 'uint64');
            end
            for k = 1:size(values, 3)                 % column-major: one whole grid at a time
                g = values(:, :, k)';                 % longitude varies fastest
                fwrite(fid, g(:), 'double');
            end
        end
        function [f, vals] = groopsFixture(folder, name, type)
            %GROOPSFIXTURE A 8 x 4 grid with values that are their own index, so a transposed or
            % wrongly ordered read is visible in the number and not only in a statistic.
            lon = -157.5:45:157.5; lat = 67.5:-45:-67.5;
            n = numel(lat) * numel(lon);
            nk = 3;
            vals = reshape(1:(n * nk), numel(lon), numel(lat), nk);
            vals = permute(vals, [2 1 3]);
            mjd = [];
            if type == "griddedDataTimeSeries", mjd = 58000 + (0:nk) * 30; end
            f = tNgfIo.writeGroopsGrid(folder, name, type, lon, lat, vals, mjd);
        end
    end
    methods (Test, TestTags = {'contract'})
        function importGroopsGridContract(tc)
            d = tc.tmp;
            [f, want] = tNgfIo.groopsFixture(d, 'g.dat', "griddedData");
            [g, aux] = NonGravForces.io.importGroopsGrid(f);
            tc.verifyEqual(g.type, "griddedData"); tc.verifyEqual(g.version, "20200123");
            tc.verifyFalse(g.isTimeSeries); tc.verifyEmpty(g.mjd);
            tc.verifySize(g.values, [4 8 3]); tc.verifySize(g.area_sr, [4 8]);
            tc.verifyEqual(numel(g.lon_deg), 8); tc.verifyEqual(numel(g.lat_deg), 4);
            tc.verifyEqual(g.values, want, 'the values must come back exactly as written, in the same order');
            tc.verifyEqual(aux.pointCount, 32); tc.verifyEqual(aux.gridCount, 3);
            tc.verifyEqual(aux.ellipsoid.a, 6378137.0);
            tc.verifyEqual(aux.provenance.unit, "importGroopsGrid");
            [s, sx] = NonGravForces.io.importGroopsGrid(tNgfIo.groopsFixture(d, 's.dat', "griddedDataTimeSeries"));
            tc.verifyTrue(s.isTimeSeries);
            tc.verifyEqual(numel(s.mjd), 4, 'a series carries one more node than it has grids');
            tc.verifyEqual(sx.gridCount, 3);
            tc.verifyEqual(s.mjd, 58000 + (0:3) * 30);
        end
    end
    methods (Test, TestTags = {'reference'})
        function importGroopsGridReference(tc)
            % The oracle is the format itself, asserted by the three facts that make a misreading
            % visible; each is checked on a file whose contents are known independently.
            d = tc.tmp;
            [f, want] = tNgfIo.groopsFixture(d, 'r.dat', "griddedDataTimeSeries");
            [g, aux] = NonGravForces.io.importGroopsGrid(f);
            tc.verifyBoundAndRecord("GROOPS grid: solid angles sum to 4 pi", ...
                abs(aux.areaSum_sr / (4 * pi) - 1), 1e-12, "-", "the writer lays a cosine-weighted grid on the sphere; the reader must not rescale it");
            tc.verifyEqual(numel(g.mjd) - 1, aux.gridCount, 'nodes minus one is the grid count');
            tc.verifyEqual(g.values, want, 'column-major: one whole grid at a time, longitude fastest');
            % Latitude decreases and longitude increases, as both real products do.
            tc.verifyTrue(all(diff(g.lat_deg) < 0) && all(diff(g.lon_deg) > 0));
        end
    end
    methods (Test, TestTags = {'precision'})
        function importGroopsGridPrecision(tc)
            % float64 in, float64 out: a read that rescales or single-casts shows here and nowhere else.
            d = tc.tmp;
            lon = -135:90:135; lat = 45:-90:-45;
            v = zeros(2, 4); v(:) = [pi; exp(1); 1/3; 1e-300; 1e300; -pi; realmin; 0.1];
            f = tNgfIo.writeGroopsGrid(d, 'p.dat', "griddedData", lon, lat, v, []);
            g = NonGravForces.io.importGroopsGrid(f);
            tc.verifyBoundAndRecord("GROOPS grid values survive the round trip", ...
                max(abs(g.values(:) - v(:))), 0, "-", "float64 written and read unchanged: bitwise or the reader is lossy");
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function importGroopsGridVectorization(tc)
            % Reading a subrange must give exactly what reading everything gives for those grids:
            % the seek arithmetic is the only difference and it must change no number.
            d = tc.tmp;
            f = tNgfIo.groopsFixture(d, 'v.dat', "griddedDataTimeSeries");
            all3 = NonGravForces.io.importGroopsGrid(f);
            [one, aux] = NonGravForces.io.importGroopsGrid(f, 'TimeRange', [58061 58062]);
            tc.verifyEqual(aux.gridsRead, 3, 'MJD 58061 falls in the third month');
            tc.verifyEqual(one.values, all3.values(:, :, 3), 'a seeked grid equals the same grid read in bulk');
            [two, aux2] = NonGravForces.io.importGroopsGrid(f, 'TimeRange', [58035 58065]);
            tc.verifyEqual(aux2.gridsRead, [2 3]);
            tc.verifyEqual(two.values, all3.values(:, :, 2:3));
            % .mjd IS TRIMMED TO THE GRIDS RETURNED, so a seeked struct says which months it holds.
            % It used to carry the file's whole node list beside a subset of values, and
            % env.earthRadiance's ceresClimatology backend could not tell one grid from another
            % (2413 nodes beside 2 grids, 2026-09-17). The invariant the reference test asserts --
            % nodes are one more than grids -- now holds after a seek too, which is the point.
            tc.verifyEqual(one.mjd, all3.mjd(3:4), 'the nodes bounding the one grid returned');
            tc.verifyEqual(two.mjd, all3.mjd(2:4));
            tc.verifyEqual(numel(one.mjd), size(one.values, 3) + 1);
            tc.verifyEqual(numel(two.mjd), size(two.values, 3) + 1);
            % and the FILE's own span stays reportable, as a fact about the file rather than these values
            tc.verifyEqual(aux.mjdFile, all3.mjd([1 end]));
            tc.verifyEqual(aux2.mjdFile, all3.mjd([1 end]));
        end
    end
    methods (Test, TestTags = {'robustness'})
        function importGroopsGridRobustness(tc)
            d = tc.tmp;
            f = tNgfIo.groopsFixture(d, 'x.dat', "griddedData");
            tc.verifyErrorId(@() NonGravForces.io.importGroopsGrid(fullfile(d, 'absent.dat')), 'NonGravForces:importGroopsGrid:inputFile');
            % An HTML error page saved as .dat is the failure the magic check exists for.
            html = fullfile(d, 'html.dat');
            fid = fopen(html, 'w'); fwrite(fid, '<!DOCTYPE html><html>404</html>', 'char'); fclose(fid);
            tc.verifyErrorId(@() NonGravForces.io.importGroopsGrid(html), 'NonGravForces:importGroopsGrid:dataFormat');
            % A truncated file: the payload no longer fills it exactly.
            cut = fullfile(d, 'cut.dat');
            % BYTES, not text. fileread decodes, and a decoded binary is corrupted in ways that have
            % nothing to do with the truncation this means to test.
            fi = fopen(f, 'r'); b = fread(fi, inf, '*uint8'); fclose(fi);
            fo = fopen(cut, 'w'); fwrite(fo, b(1:end-12), 'uint8'); fclose(fo);
            tc.verifyErrorId(@() NonGravForces.io.importGroopsGrid(cut), 'NonGravForces:importGroopsGrid:dataFormat');
            % TimeRange on a file with no time, and a range outside the series.
            tc.verifyErrorId(@() NonGravForces.io.importGroopsGrid(f, 'TimeRange', [1 2]), 'NonGravForces:importGroopsGrid:inputRange');
            s = tNgfIo.groopsFixture(d, 'y.dat', "griddedDataTimeSeries");
            tc.verifyErrorId(@() NonGravForces.io.importGroopsGrid(s, 'TimeRange', [1 2]), 'NonGravForces:importGroopsGrid:inputRange');
            tc.verifyErrorId(@() NonGravForces.io.importGroopsGrid(s, 'TimeRange', [58100 58000]), 'NonGravForces:importGroopsGrid:inputRange');
            tc.verifyErrorId(@() NonGravForces.io.importGroopsGrid(s, 'TimeRange', 58000), 'NonGravForces:importGroopsGrid:inputRange');
            % An unknown GROOPS type is refused BY NAME rather than misread.
            other = fullfile(d, 'o.dat');
            fid = fopen(other, 'w', 'l');
            fwrite(fid, 'B', 'char'); fwrite(fid, 10, 'uint64'); fwrite(fid, 'instrument', 'char');
            fwrite(fid, 8, 'uint64'); fwrite(fid, '20200123', 'char'); fclose(fid);
            tc.verifyErrorId(@() NonGravForces.io.importGroopsGrid(other), 'NonGravForces:importGroopsGrid:unsupportedType');
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function importGroopsGridMetamorphic(tc)
            % Writing what was read must reproduce the file: the layout is then self-consistent,
            % which no single-direction test can establish.
            d = tc.tmp;
            f = tNgfIo.groopsFixture(d, 'm1.dat', "griddedDataTimeSeries");
            g = NonGravForces.io.importGroopsGrid(f);
            again = tNgfIo.writeGroopsGrid(d, 'm2.dat', "griddedDataTimeSeries", ...
                g.lon_deg, g.lat_deg, g.values, g.mjd);
            h = NonGravForces.io.importGroopsGrid(again);
            tc.verifyEqual(h.values, g.values); tc.verifyEqual(h.mjd, g.mjd);
            tc.verifyEqual(h.lon_deg, g.lon_deg); tc.verifyEqual(h.lat_deg, g.lat_deg);
            tc.verifyBoundAndRecord("GROOPS grid: read-write-read is a fixed point", ...
                max(abs(h.values(:) - g.values(:))), 0, "-", "bitwise; the two files describe the same grid");
        end
    end
    methods (Test, TestTags = {'speed'})
        function importGroopsGridSpeed(tc)
            % SEEKING IS THE POINT. One month of a many-month series must not cost what the whole
            % series costs, or a 1.25 GB file is unusable for a one-day arc.
            d = tc.tmp;
            % Big enough that the VALUES dominate: 180 x 90 x 60 grids is 7.8 MB, against 130 kB
            % for one grid. At the first attempt the fixture was so small that the header decided the
            % time and the ratio measured noise (2.4, meaninglessly).
            lon = -178:2:180; lat = 89:-2:-89;
            nk = 60;
            v = rand(numel(lat), numel(lon), nk);
            big = tNgfIo.writeGroopsGrid(d, 'spBig.dat', "griddedDataTimeSeries", lon, lat, v, 58000 + (0:nk) * 30);
            small = tNgfIo.writeGroopsGrid(d, 'spSmall.dat', "griddedDataTimeSeries", lon, lat, v(:, :, 1:2), 58000 + (0:2) * 30);
            cleanupBig = onCleanup(@() delete(big, small));
            one = @(f) NonGravForces.io.importGroopsGrid(f, 'TimeRange', [58000 58001]);
            one(big); one(small);
            tBig = zeros(1, 5); tSmall = zeros(1, 5);
            for k = 1:5
                t0 = tic; one(big);   tBig(k) = toc(t0);
                t0 = tic; one(small); tSmall(k) = toc(t0);
            end
            tc.verifyBoundAndRecord("GROOPS grid: one grid of 60 / one grid of 2", median(tBig) / median(tSmall), 3, "-", ...
                "seeking means the cost of one grid does not grow with the file: 30x the payload must not cost 30x. A reader that ignored TimeRange would land near 30");
        end
    end

end

% ------------------------------------------------------------------------- local helpers
function L = localSetToken(L, k, j, value)
tok = split(strtrim(L(k))); tok(j) = string(value); L(k) = strjoin(tok, " ");
end
function L = localScaleToken(L, k, j, factor)
tok = split(strtrim(L(k))); tok(j) = string(factor * str2double(tok(j))); L(k) = strjoin(tok, " ");
end
function L = localZeroTokens(L, k, jj)
tok = split(strtrim(L(k))); tok(jj) = "0"; L(k) = strjoin(tok, " ");
end
function x = localBaseline(mjd, ep)
[~, r] = ismember(ep.day, mjd); s = floor(ep.sec / 10800) + 1; x = r + s;
end
