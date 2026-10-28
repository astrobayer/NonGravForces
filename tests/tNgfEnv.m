classdef tNgfEnv < NgfTestCase
%TNGFENV Tests of the env layer (LP-3): density, eclipseFactor, earthRadiance, radiationSources, and buildMex.
%
%   Oracles (none built here): the NRLMSISE-00 DOCUMENTATION §7 table (17
%   published cases, parsed by this class from whichever release is present in
%   mex/nrlmsise00 -- the oracle travels with the source, so tracking upstream
%   cannot silently invalidate it); NRL's msis2.1_test_in.txt and
%   msis2.1_test_ref_dp.txt (200 cases, read from the NRLMSIS 2.1 package the
%   msis21 MEX was built from, for the same reason);
%   tools/pyxcheck/eclipse_check.py (pure-Python disc-overlap integration);
%   tools/pyxcheck/knocke_check.py (pure-Python 0.25 deg cap integration and an
%   independent restatement of Knocke's Table 1); the real 2020-01-01 GRACE-C
%   state from data/samples. Every bound carries the number measured on
%   2026-09-07.
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
    properties
        root
        msisTable    % 17x11 published values (cgs): HE O N2 O2 AR RHO H N ANOM TINF TG
        msisInput    % struct of 1x17 inputs from nrlmsise-00_test.c
        msis21In     % 200x9 NRL inputs: iyd sec alt glat glong stl f107a f107 Ap
        msis21Ref    % 200x11 NRL outputs (cgs): He O N2 O2 Ar rho H N O* NO T
        cer          % io.importCERES of the EBAF excerpt in data/ceres (2019-12, 2020-01)
        igrf         % io.igrfCoefficients of data/igrf/igrf14coeffs.txt at 2020.0
        eclipse      % oracle table
        knocke       % oracle table
        au
        cLight
    end
    methods (TestClassSetup)
        function loadOracles(tc)
            tc.root = fileparts(fileparts(mfilename('fullpath')));
            % The oracle travels WITH the source: section 7 of the package's own DOCUMENTATION
            % is its published testing output, so whichever NRLMSISE-00 release is present is
            % checked against its own table. Since 2026-09-09 that source is not in this
            % repository (O-34) and the run needs it fetched first.
            docPath = fullfile(tc.root, 'mex', 'nrlmsise00', 'DOCUMENTATION');
            tc.assertTrue(isfile(docPath), ...
                sprintf(['the NRLMSISE-00 source is not at %s. It is not redistributed with this ' ...
                         'toolbox (THIRD-PARTY-NOTICES.md): run setupNonGravForces, which fetches ' ...
                         'it at the pinned commit the first time it runs.'], docPath));
            L = splitlines(string(fileread(docPath)));
            rows = L(~cellfun(@isempty, regexp(L, '^\d\.\d{6}E[+-]\d{2} ', 'once')));
            tc.assertEqual(numel(rows), 17, 'DOCUMENTATION §7 must hold the 17 published rows');
            tc.msisTable = cell2mat(cellfun(@(s) sscanf(s, '%f')', cellstr(rows), 'UniformOutput', false));
            % inputs: the base case and the deviations, transcribed from nrlmsise-00_test.c (published driver)
            n = 17; in.doy = 172*ones(1,n); in.sec = 29000*ones(1,n); in.alt_km = 400*ones(1,n); in.lat = 60*ones(1,n); in.lon = -70*ones(1,n);
            in.lst = 16*ones(1,n); in.f107a = 150*ones(1,n); in.f107 = 150*ones(1,n); in.ap = 4*ones(1,n);
            in.doy(2) = 81; in.sec(3) = 75000; in.alt_km(3) = 1000; in.alt_km(4) = 100; in.alt_km(11) = 0; in.alt_km(12) = 10; in.alt_km(13) = 30;
            in.alt_km(14) = 50; in.alt_km(15) = 70; in.alt_km(17) = 100; in.lat(5) = 0; in.lon(6) = 0; in.lst(7) = 4; in.f107a(8) = 70; in.f107(9) = 180; in.ap(10) = 40;
            tc.msisInput = in;
            tc.eclipse = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'eclipse_oracle.csv'), 'VariableNamingRule', 'preserve');
            tc.knocke  = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'knocke_oracle.csv'), 'VariableNamingRule', 'preserve');
            tc.au = NonGravForces.const().au;
            tc.cLight = referenceSystem.iersConstants().speedOfLight;
            tc.assertTrue(NonGravForces.buildMex('Check', true, 'Target', "msis").present, 'the MSIS MEX must be built before this class runs (ciGate builds it)');
            % NRLMSIS 2.1: the oracle is read from the package the MEX was built from (buildMex reports where).
            r21 = NonGravForces.buildMex('Check', true, 'Target', "msis21");
            tc.assertTrue(r21.upToDate, 'the NRLMSIS 2.1 MEX must be built and current before this class runs (setupNonGravForces fetches the package, ciGate builds it)');
            tc.msis21In = readmatrix(fullfile(r21.sourceDir, 'msis2.1_test_in.txt'), 'FileType', 'text', 'NumHeaderLines', 1);
            tc.msis21Ref = readmatrix(fullfile(r21.sourceDir, 'msis2.1_test_ref_dp.txt'), 'FileType', 'text', 'NumHeaderLines', 1);
            tc.assertSize(tc.msis21In, [200 9], 'msis2.1_test_in.txt holds 200 cases of 9 inputs');
            tc.assertSize(tc.msis21Ref, [200 20], 'msis2.1_test_ref_dp.txt repeats the 9 inputs and prints 11 outputs');
            tc.assertEqual(tc.msis21Ref(:, 1:2), tc.msis21In(:, 1:2), 'the reference rows are the input rows, in order');
            tc.msis21Ref = tc.msis21Ref(:, 10:20);
            % DTM2020 (2026-10-27): the oracle is CNES's own benchmark (Mod_DTM2020F107Kp_Benchmark.f90), transcribed below.
            tc.assertTrue(NonGravForces.buildMex('Check', true, 'Target', "dtm2020").upToDate, ...
                'the DTM2020 MEX must be built and current before this class runs (setupNonGravForces fetches CNES''s files, ciGate builds them)');
            tc.igrf = NonGravForces.io.igrfCoefficients(fullfile(tc.root, 'data', 'igrf', 'igrf14coeffs.txt'), 'Epoch', 2020.0);
            tc.cer = NonGravForces.io.importCERES(string(fullfile(tc.root, 'data', 'ceres', ...
                ["rsut" "rlut" "rsdt"] + "_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_201912-202001.nc")));
        end
    end
    methods
        function [ep, lla, drv, lst, ref] = dtmBenchmark(~)
            %DTMBENCHMARK CNES's four DTM2020 operational benchmark cases (Mod_DTM2020F107Kp_Benchmark.f90,
            % S. Bruinsma 22/10/20): day 180, latitude and longitude 0, local time 3.1415 rad, Kp 3 and its
            % 24 h mean 3; F10.7 = its 81-day mean = 80 or 180; altitude 300 or 800 km. Order as printed:
            % (80, 300), (180, 300), (80, 800), (180, 800). ref: tz [K], tinf [K], ro [g/cm^3], as printed.
            ep = referenceSystem.daySecondToEpoch([59394 * ones(4, 1), 43200 * ones(4, 1)], "UTC");   % 2021-06-29: day 180
            lla = [zeros(2, 4); [300 300 800 800] * 1e3];
            f = [80 180 80 180];
            drv = struct('f107', f, 'f107bar', f, 'kp3', 3, 'kp24', 3);
            lst = 3.1415 * 12 / pi;
            ref = struct('tz', [843.243 1312.266 844.099 1334.874], 'tinf', [844.099 1334.875 844.099 1334.875], ...
                         'ro', [0.94719e-14 0.29053e-13 0.26600e-17 0.39104e-16]);
        end
        function [rho, sp, T] = msisCases(tc, idx, anomalous)
            % Evaluate the published cases idx (1-based) with the published (inconsistent) local time.
            in = tc.msisInput;
            ep = referenceSystem.daySecondToEpoch([51544 + in.doy(idx)' - 1, in.sec(idx)'], "UTC");   % year 2000: irrelevant to the model
            lla = [deg2rad(in.lat(idx)); deg2rad(in.lon(idx)); in.alt_km(idx) * 1e3];
            if all(idx >= 16)
                ap = 100 * ones(7, numel(idx));                    % cases 16, 17: the ap array (all 100), switch 9 = -1
            else
                ap = in.ap(idx);
            end
            d = struct('f107', in.f107(idx), 'f107a', in.f107a(idx), 'ap', ap);
            [rho, sp, T] = NonGravForces.env.density("msis00", ep, lla, [], 'Drivers', d, 'AnomalousOxygen', anomalous, 'LocalSolarTime_h', in.lst(idx));
        end
        function [rho, sp, T] = msis21Cases(tc, idx, varargin)
            % NRL's test cases idx through env.density. The test program reads its inputs as real(4)
            % and widens them, so the same single-precision values are handed over here; the year is
            % irrelevant to the model (only the day of year enters) and 2001 has no day 366.
            X = double(single(tc.msis21In(idx, :)));
            doy = mod(X(:, 1), 1000);
            ep = referenceSystem.daySecondToEpoch([51910 + doy - 1, X(:, 2)], "UTC");
            lla = [deg2rad(X(:, 4))'; deg2rad(X(:, 5))'; X(:, 3)' * 1e3];
            d = struct('f107', X(:, 8)', 'f107a', X(:, 7)', 'ap', X(:, 9)');
            [rho, sp, T] = NonGravForces.env.density("msis21", ep, lla, [], 'Drivers', d, varargin{:});
        end
        function c = syntheticCeres(tc, fields)
            %SYNTHETICCERES The excerpt's grid and months with caller-defined fields: fields(lat, lon, month) -> [rsut rlut rsdt].
            c = tc.cer;
            [LON, LAT] = meshgrid(c.lon_deg, c.lat_deg);
            for m = 1:numel(c.mjd)
                [a, b, d] = fields(LAT, LON, m);
                c.rsut(:, :, m) = single(a); c.rlut(:, :, m) = single(b); c.rsdt(:, :, m) = single(d);
            end
        end
        function clm = syntheticClimatology(~, fAlb, fLw, nLat, nLon, mjd)
            %SYNTHETICCLIMATOLOGY Two io.importGroopsGrid-shaped structs with caller-defined fields.
            %   THE LATITUDES ARE GEOCENTRIC AND THE RADII FOLLOW THEM, because that is what the real
            %   file holds and what localCheckClimatology measures: the grid is laid out uniformly in
            %   GEODETIC latitude and the file stores the geocentric latitude of those same points.
            %   Building the fixture the other way round would make the unit look wrong for a fault
            %   that was the fixture's -- the 2.5e-2 solid-angle lesson of WP-E1, one layer up.
            ell = referenceSystem.ellipsoidParameters("WGS84");
            latGeod = 90 - 180 * ((1:nLat) - 0.5) / nLat;                  % uniform GEODETIC centres
            latGeoc = atand((1 - ell.e2) * tand(latGeod));                 % what the file stores
            lon = -180 + 360 * ((1:nLon) - 0.5) / nLon;
            th = deg2rad(latGeoc(:));
            radius = ell.a * ell.b ./ sqrt((ell.a * sin(th)).^2 + (ell.b * cos(th)).^2) + 30000.1;
            [LON, LAT] = meshgrid(lon, latGeoc);
            T = numel(mjd) - 1;
            va = zeros(nLat, nLon, T); vl = zeros(nLat, nLon, T);
            for m = 1:T
                va(:, :, m) = fAlb(LAT, LON, m);
                vl(:, :, m) = fLw(LAT, LON, m);
            end
            dPhi = pi / nLat; dLam = 2 * pi / nLon;
            phi = deg2rad(latGeoc(:));
            area = repmat(dLam * (sin(phi + dPhi / 2) - sin(phi - dPhi / 2)), 1, nLon);
            mk = @(v) struct('type', "griddedDataTimeSeries", 'version', "test", ...
                'lon_deg', lon, 'lat_deg', latGeoc, 'radius_m', radius(:)', 'area_sr', area, ...
                'values', v, 'mjd', mjd, 'isTimeSeries', true);
            clm = struct('albedo', mk(va), 'longwave', mk(vl));
        end
        function R = eyes(~, K)
            R = repmat(eye(3), 1, 1, K);
        end
        function [info, rG, ep] = hwmAt(tc, eop, doy, sec, glat, glon, alt_km, varargin)
            % env.wind "hwm14" at 1995 day-of-year/UT and geodetic points given in degrees and km (NRL's check
            % driver inputs, single precision as the driver declares them); the GCRS position goes through the
            % same rotation the unit applies, so the model sees the geodetic point back to rounding.
            n = numel(glat);
            X = double(single([repmat(sec, 1, n / numel(sec)); glat; glon; alt_km]));
            ep = referenceSystem.daySecondToEpoch([repmat(49718 + doy - 1, n, 1), X(1, :)'], "UTC");      % MJD 49718 = 1995-01-01
            rI = referenceSystem.geodeticToCartesian([X(2, :)' X(3, :)' X(4, :)' * 1e3], 'Ellipsoid', "WGS84", 'AngleUnit', "deg")';
            R = referenceSystem.rotationItrsToGcrs(ep, 'Eop', eop);
            rG = reshape(pagemtimes(R, reshape(rI, 3, 1, n)), 3, n);
            [~, info] = NonGravForces.env.wind("hwm14", ep, rG, zeros(3, n), 'Eop', eop, varargin{:});
        end
        function [ep, r, q, eop, v] = sampleState(tc)
            eop = referenceSystem.readEop(fullfile(tc.root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            [ep, r, v] = NonGravForces.io.readGNV1B(fullfile(tc.root, 'data', 'samples', 'GNV1B_2020-01-01_C_04.txt'), 'Eop', eop);
            [~, q] = NonGravForces.io.readSCA1B(fullfile(tc.root, 'data', 'samples', 'SCA1B_2020-01-01_C_04.txt'));
        end
        function [Ea, Ei, src] = knockeCase(tc, k, grid, backend)
            K = tc.knocke;
            ep = referenceSystem.datetimeToEpoch(datetime(1981, 12, 22) + days(K.days_since_t0(k)), "UTC");
            r = [K.rx(k); K.ry(k); K.rz(k)]; sun = tc.au * [K.sunx(k); K.suny(k); K.sunz(k)];
            opt = NonGravForces.Options("standard", "erp.backend", backend, "srp.eclipse", "cylindrical");
            if ~isempty(grid), opt = opt.set("erp.grid", grid); end
            src = NonGravForces.env.radiationSources(ep, r, [1; 0; 0; 0], opt, 'Sun', sun);
            Mt = (src.M - 1) / 2;
            Ea = sum(src.P(2:1+Mt)) * tc.cLight; Ei = sum(src.P(2+Mt:end)) * tc.cLight;
        end
        function q = speedRatio(~, fcn, base)
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
        function windContract(tc)
            [ep, r, ~, eop, v] = tc.sampleState();
            [vrel, info] = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
            tc.verifySize(vrel, size(v)); tc.verifyEqual(info.model, "corotation");
            tc.verifySize(info.corotation_m_s, size(v)); tc.verifyEqual(info.K, size(v, 2));
            tc.verifyFalse(info.windGiven);
            [vn, iN] = NonGravForces.env.wind("none", ep, r, v);
            tc.verifyEqual(vn, v, 'the none model returns the inertial velocity untouched');
            tc.verifyEqual(iN.corotation_m_s, zeros(size(v)));
            w = @(varargin) NonGravForces.env.wind(varargin{:});
            % hwm14: the co-rotation plus the model wind, and its inputs refused wherever nothing reads them
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            [vh, ih] = w("hwm14", ep, r, v, 'Eop', eop, 'SpaceWeather', sw);
            tc.verifySize(ih.hwm14.northEast_m_s, [2 size(v, 2)]); tc.verifySize(ih.hwm14.gcrs_m_s, size(v)); tc.verifySize(ih.hwm14.geodetic, size(v));
            tc.verifyTrue(ih.hwm14.disturbed); tc.verifySize(ih.hwm14.ap, [1 size(v, 2)]); tc.verifyEqual(ih.model, "hwm14");
            tc.verifyEqual(vh, v - ih.corotation_m_s - ih.hwm14.gcrs_m_s, 'vrel is v minus the co-rotation minus the model wind');
            tc.verifyEmpty(info.hwm14, 'no model wind for corotation');
            tc.verifyErrorId(@() w("hwm14", ep, r, v, 'Eop', eop), 'NonGravForces:wind:inputAp');
            tc.verifyErrorId(@() w("hwm14", ep, r, v, 'Eop', eop, 'SpaceWeather', sw, 'Ap', 4), 'NonGravForces:wind:inputAp');
            tc.verifyErrorId(@() w("hwm14", ep, r, v, 'Eop', eop, 'Ap', -1), 'NonGravForces:wind:inputAp');
            tc.verifyErrorId(@() w("hwm14", ep, r, v, 'Eop', eop, 'Disturbed', false, 'Ap', 4), 'NonGravForces:wind:inputAp');
            tc.verifyErrorId(@() w("corotation", ep, r, v, 'Eop', eop, 'Ap', 4), 'NonGravForces:wind:inputAp');
            tc.verifyErrorId(@() w("hwm14", ep, r, v, 'Disturbed', false), 'NonGravForces:wind:inputEop');
            tc.verifyErrorId(@() w("corotation", ep, r, v), 'NonGravForces:wind:inputEop');
            tc.verifyErrorId(@() w("corotation", ep, r(:, 1:2), v, 'Eop', eop), 'NonGravForces:wind:inputShape');
            tc.verifyErrorId(@() w("none", ep, r, v, 'Wind_m_s', zeros(3, 2)), 'NonGravForces:wind:inputShape');
            bad = r; bad(1) = NaN;
            tc.verifyErrorId(@() w("none", ep, bad, v), 'NonGravForces:wind:numericNaN');
            tc.verifyErrorId(@() w("none", ep, r, v, 'Wind_m_s', nan(size(v))), 'NonGravForces:wind:numericNaN');
            tc.verifyError(@() w("hwm07", ep, r, v), 'MATLAB:validators:mustBeMember');
        end
        function densityDtm2020Contract(tc)
            [ep, lla, drv, lst] = tc.dtmBenchmark();
            [rho, sp, T, info] = NonGravForces.env.density("dtm2020", ep, lla, [], 'Drivers', drv, 'LocalSolarTime_h', lst);
            tc.verifySize(rho, [1 4]); tc.verifySize(sp, [8 4]); tc.verifySize(T, [2 4]); tc.verifyTrue(all(rho > 0));
            tc.verifyEqual(info.backend, "dtm2020"); tc.verifyEqual(info.drivers.source, "caller"); tc.verifyFalse(info.anomalousOxygen);
            tc.verifyEqual(sp([5 8], :), zeros(2, 4), 'DTM2020 has no argon and no anomalous oxygen');
            tc.verifyEqual(info.dtmInput.doy, 180 * ones(1, 4)); tc.verifyEqual(info.dtmInput.alt_km, [300 300 800 800]);
            % from a celestrak spaceWeather struct: the query's DTM fields, unchanged
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            e2 = referenceSystem.mjdToEpoch([58849.5; 58850.25], "UTC"); l2 = [0 0.5; 1 2; 490e3 495e3];
            [~, ~, ~, i2] = NonGravForces.env.density("dtm2020", e2, l2, sw);
            q = sw.query(e2);
            tc.verifyEqual([i2.drivers.f107; i2.drivers.f107bar; i2.drivers.kp3; i2.drivers.kp24], [q.f107; q.f107Last81; q.kpDelayed3h; q.kpMean24h]);
            tc.verifyEqual(i2.drivers.source, "spaceWeather");
            % refusals, each paired with the call that must succeed
            e1 = referenceSystem.mjdToEpoch(58849.5, "UTC");
            d1 = struct('f107', 80, 'f107bar', 80, 'kp3', 3, 'kp24', 3);
            tc.verifyErrorId(@() NonGravForces.env.density("dtm2020", e1, [0; 0; 149.9e3], [], 'Drivers', d1), 'NonGravForces:density:inputRange');
            tc.verifyTrue(NonGravForces.env.density("dtm2020", e1, [0; 0; 150e3], [], 'Drivers', d1) > 0, '150 km is DTM2020 and evaluates');
            tc.verifyErrorId(@() NonGravForces.env.density("dtm2020", e1, [0; 0; 490e3], [], 'Drivers', struct('f107', 80, 'f107a', 80, 'ap', 4)), 'NonGravForces:density:inputDrivers');
            tc.verifyErrorId(@() NonGravForces.env.density("dtm2020", e1, [0; 0; 490e3], []), 'NonGravForces:density:inputDrivers');
            tc.verifyErrorId(@() NonGravForces.env.density("dtm2020", e1, [0; 0; 490e3], [], 'Drivers', setfield(d1, 'kp3', NaN)), 'NonGravForces:density:numericNaN'); %#ok<SFLD>
        end
        function densityContract(tc)
            ep = referenceSystem.mjdToEpoch([58849.5; 58849.6], "UTC"); lla = [0 0.1; 0 0.2; 490e3 500e3];
            d = struct('f107', [70 71], 'f107a', [71 71], 'ap', [3 4]);
            [rho, sp, T, info] = NonGravForces.env.density("msis00", ep, lla, [], 'Drivers', d);
            tc.verifySize(rho, [1 2]); tc.verifySize(sp, [8 2]); tc.verifySize(T, [2 2]); tc.verifyTrue(all(rho > 0));
            tc.verifyEqual(info.backend, "msis00"); tc.verifyEqual(info.drivers.source, "caller"); tc.verifyTrue(info.anomalousOxygen);
            tc.verifyEqual(info.provenance.unit, "env.density"); tc.verifyEqual(info.msisInput.alt_km, [490 500]);
            p = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);
            [rhoE, spE, TE, infoE] = NonGravForces.env.density("exponential", ep, lla, [], 'Exponential', p);
            tc.verifySize(rhoE, [1 2]); tc.verifyTrue(all(isnan(spE(:)))); tc.verifyTrue(all(isnan(TE(:)))); tc.verifyEqual(infoE.drivers.source, "none");
            tc.verifyError(@() NonGravForces.env.density("nrlmsis", ep, lla, []), 'MATLAB:validators:mustBeMember');
            tc.verifyErrorId(@() NonGravForces.env.density("msis00", ep, lla(1:2, :), [], 'Drivers', d), 'NonGravForces:density:inputShape');
            tc.verifyErrorId(@() NonGravForces.env.density("msis00", ep, lla, []), 'NonGravForces:density:inputDrivers');
            tc.verifyErrorId(@() NonGravForces.env.density("exponential", ep, lla, []), 'NonGravForces:density:inputExponential');
            [rho21, sp21, T21, info21] = NonGravForces.env.density("msis21", ep, lla, [], 'Drivers', d);
            tc.verifySize(rho21, [1 2]); tc.verifySize(sp21, [8 2]); tc.verifySize(T21, [2 2]); tc.verifyTrue(all(rho21 > 0) && all(sp21(:) > 0));
            tc.verifyEqual(info21.backend, "msis21"); tc.verifyEqual(info21.drivers.source, "caller"); tc.verifyTrue(info21.anomalousOxygen);
            tc.verifyEqual(info21.msisInput.alt_km, [490 500]); tc.verifyFalse(isfield(info21.msisInput, 'lst_h'), 'msis21 forms local time itself');
            tc.verifyErrorId(@() NonGravForces.env.density("msis21", ep, lla, [], 'Drivers', d, 'LocalSolarTime_h', [16 16]), 'NonGravForces:density:inputOption');
            tc.verifyErrorId(@() NonGravForces.env.density("msis21", ep, lla, []), 'NonGravForces:density:inputDrivers');
            tc.verifyError(@() NonGravForces.env.density("jb2008", ep, lla, [], 'Drivers', d), 'MATLAB:validators:mustBeMember');
        end
        function magneticFieldContract(tc)
            r = [[6861.2e3; 0; 0] [0; 5e6; 4.8e6]];
            [B, info] = NonGravForces.env.magneticField(tc.igrf, r);
            tc.verifySize(B, [3 2]); tc.verifyEqual(info.itrs_T, B); tc.verifySize(info.northEastDown_T, [3 2]);
            tc.verifyEqual(info.intensity_T, vecnorm(B)); tc.verifyEqual(info.maxDegree, 13); tc.verifyEqual(info.year, 2020.0);
            tc.verifyTrue(all(info.intensity_T > 2e-5 & info.intensity_T < 6e-5), 'tens of microtesla at LEO');
            tc.verifyEqual(info.provenance.unit, "env.magneticField");
            M = @(varargin) NonGravForces.env.magneticField(varargin{:});
            tables = NonGravForces.io.igrfCoefficients(fullfile(tc.root, 'data', 'igrf', 'igrf14coeffs.txt'));
            tc.verifyErrorId(@() M(tables, r), 'NonGravForces:magneticField:inputCoefficients');
            tc.verifyErrorId(@() M(tc.igrf, r(1:2, :)), 'NonGravForces:magneticField:inputShape');
            tc.verifyErrorId(@() M(tc.igrf, r, 'Rotation', eye(2)), 'NonGravForces:magneticField:inputShape');
            tc.verifyErrorId(@() M(tc.igrf, r, 'MaxDegree', 14), 'NonGravForces:magneticField:inputRange');
            tc.verifyErrorId(@() M(tc.igrf, r, 'MaxDegree', 0), 'NonGravForces:magneticField:inputRange');
            tc.verifyErrorId(@() M(tc.igrf, [r zeros(3, 1)]), 'NonGravForces:magneticField:inputRange');
            tc.verifyErrorId(@() M(tc.igrf, [NaN; 0; 7e6]), 'NonGravForces:magneticField:numericNaN');
            tc.verifyErrorId(@() M(tc.igrf, r, 'Rotation', nan(3)), 'NonGravForces:magneticField:numericNaN');
        end
        function eclipseFactorContract(tc)
            % at 490 km the penumbra spans 67.97..68.51 deg from the anti-Sun point (b -/+ a); 68.25 is inside it
            r = [-(6378137 + 490e3) * cosd([0 68.25 90]); 0 0 0; (6378137 + 490e3) * sind([0 68.25 90])]; s = [tc.au; 0; 0];
            [nu, info] = NonGravForces.env.eclipseFactor(r, s);
            tc.verifySize(nu, [1 3]); tc.verifyEqual(nu(1), 0); tc.verifyTrue(nu(2) > 0 && nu(2) < 1); tc.verifyEqual(nu(3), 1);
            tc.verifyEqual(info.model, "conical"); tc.verifyEqual(info.radiusEarth_m, referenceSystem.ellipsoidParameters("WGS84").a); tc.verifyEmpty(info.moon);
            nuY = NonGravForces.env.eclipseFactor(r, s, 'Model', "cylindrical"); tc.verifyEqual(nuY, [0 1 1]);
            [nuM, infoM] = NonGravForces.env.eclipseFactor(r, s, 'Moon', [0; 384400e3; 0]); tc.verifyEqual(nuM, nu); tc.verifySize(infoM.moon, [1 3]);
            tc.verifyEqual(NonGravForces.env.eclipseFactor(r, repmat(s, 1, 3)), nu, 'a 3xK Sun equals the broadcast 3x1');
            tc.verifyErrorId(@() NonGravForces.env.eclipseFactor(r(1:2, :), s), 'NonGravForces:eclipseFactor:inputShape');
            tc.verifyErrorId(@() NonGravForces.env.eclipseFactor(r, [s s]), 'NonGravForces:eclipseFactor:inputShape');
            tc.verifyErrorId(@() NonGravForces.env.eclipseFactor(r, s, 'Model', "solaarsCF"), 'NonGravForces:eclipseFactor:inputEclipticNormal');
            tc.verifyErrorId(@() NonGravForces.env.eclipseFactor(r, s, 'Model', "solaarsCF", 'EclipticNormal', [0; 0; 2]), 'NonGravForces:eclipseFactor:inputEclipticNormal');
            [nuS, infoS] = NonGravForces.env.eclipseFactor(r, s, 'Model', "solaarsCF", 'EclipticNormal', [0; 0; 1]);
            tc.verifySize(nuS, [1 3]); tc.verifyEqual(nuS(1), 0); tc.verifyTrue(nuS(2) > 0 && nuS(2) < 1); tc.verifyEqual(nuS(3), 1);
            ell = referenceSystem.ellipsoidParameters("WGS84");
            tc.verifyEqual(infoS.solaars.clampedCount, 0); tc.verifyEqual(infoS.solaars.oblatenessScale, ell.a / ell.b);
            tc.verifyError(@() NonGravForces.env.eclipseFactor(r, s, 'Model', "flat"), 'MATLAB:validators:mustBeMember');
        end
        function earthRadianceContract(tc)
            ep = referenceSystem.mjdToEpoch([58849.5; 58850.5], "UTC"); r = [[6.88e6; 0; 0] [0; 0; 7e6]];
            t = NonGravForces.env.earthRadiance("knocke", ep, [12 24], r);
            tc.verifyEqual(t.M, 288); tc.verifySize(t.n, [3 288 2]); tc.verifySize(t.area_m2, [288 2]); tc.verifySize(t.albedo, [288 2]); tc.verifySize(t.emissivity, [288 2]);
            tc.verifySize(t.lat_rad, [288 2]); tc.verifySize(t.capAngle_rad, [1 2]); tc.verifyEqual(t.grid, [12 24]); tc.verifyEqual(t.backend, "knocke");
            tc.verifyTrue(all(t.albedo(:) > 0 & t.albedo(:) < 1) && all(t.emissivity(:) > 0 & t.emissivity(:) < 1));
            p = NonGravForces.env.earthRadiance("pointSource", ep, [], r);
            tc.verifyEqual(p.M, 0); tc.verifyEqual(p.albedo, [0.34 0.34]); tc.verifyEqual(p.emissivity, [0.68 0.68]);
            tc.verifyErrorId(@() NonGravForces.env.earthRadiance("knocke", ep, [0 24], r), 'NonGravForces:earthRadiance:inputGrid');
            tc.verifyErrorId(@() NonGravForces.env.earthRadiance("knocke", ep, [12 24], r(:, 1)), 'NonGravForces:earthRadiance:inputShape');
            tc.verifyErrorId(@() NonGravForces.env.earthRadiance("knocke", ep, [12 24], [r(:, 1) [1e6; 0; 0]]), 'NonGravForces:earthRadiance:inputRange');
            % ceres: the same cap tiles on the 20 km reference sphere, the measured fields in place of Knocke's
            [tC, rC] = deal(NonGravForces.env.earthRadiance("ceres", ep, [12 24], r, 'Ceres', tc.cer, 'RotationItrsToGcrs', tc.eyes(2)), r);
            tc.verifyEqual(tC.M, 288); tc.verifySize(tC.albedo, [288 2]); tc.verifySize(tC.exitanceLw_W_m2, [288 2]); tc.verifyEmpty(tC.emissivity);
            tc.verifyEqual(tC.radius_m, referenceSystem.ellipsoidParameters("WGS84").a + 20e3);
            tc.verifyEqual(tC.ceres.product, "CERES-EBAF-4-2-1"); tc.verifySize(tC.ceres.months, [2 2]); tc.verifySize(tC.ceres.timeWeight, [1 2]);
            tc.verifyTrue(all(tC.exitanceLw_W_m2(:) > 100 & tC.exitanceLw_W_m2(:) < 400) && all(tC.albedo(:) >= 0));
            tc.verifyEmpty(t.exitanceLw_W_m2); tc.verifyEmpty(p.exitanceLw_W_m2);
            E = @(varargin) NonGravForces.env.earthRadiance(varargin{:});
            tc.verifyErrorId(@() E("ceres", ep, [12 24], rC), 'NonGravForces:earthRadiance:inputCeres');
            tc.verifyErrorId(@() E("ceres", ep, [12 24], rC, 'Ceres', tc.cer), 'NonGravForces:earthRadiance:inputCeres');
            tc.verifyErrorId(@() E("ceres", ep, [12 24], rC, 'Ceres', tc.cer, 'RotationItrsToGcrs', tc.eyes(3)), 'NonGravForces:earthRadiance:inputCeres');
            tc.verifyErrorId(@() E("ceres", ep, [12 24], rC, 'Ceres', rmfield(tc.cer, 'rsdt'), 'RotationItrsToGcrs', tc.eyes(2)), 'NonGravForces:earthRadiance:inputCeres');
            tc.verifyErrorId(@() E("knocke", ep, [12 24], rC, 'Ceres', tc.cer), 'NonGravForces:earthRadiance:inputCeres');
            tc.verifyErrorId(@() E("ceres", referenceSystem.mjdToEpoch([58870; 58849.5], "UTC"), [12 24], rC, 'Ceres', tc.cer, 'RotationItrsToGcrs', tc.eyes(2)), 'NonGravForces:earthRadiance:dataRange');
        end
        function radiationSourcesDe440Reference(tc)
            % THE PRECISE PRESET'S OWN DEFAULT, evaluated at last. Until 2026-09-16 the only test
            % naming rs.de440 asserted that it REFUSES without kernels, so a shipped preset's default
            % sun.source had never produced a number (HANDOVER O-7, O-15).
            %
            % Oracle: JPL's DE440 through referenceSystem's MICE provider -- an ephemeris this project
            % did not compute. The claim is NOT that the two agree to rounding; they are different
            % ephemerides. It is that referenceSystem's analytic Sun is accurate enough for radiation
            % pressure, which is what justifies rs.analytic being a preset value at all.
            kernel = fullfile(tc.root, 'data', 'kernels', 'de440s.bsp');
            tc.assertTrue(isfile(kernel), ...
                'data/kernels/de440s.bsp is missing: run setupNonGravForces (FetchKernels ifMissing). This test does NOT skip -- a preset default that nothing evaluates is exactly how rs.de440 stayed unrun for two levels.');
            eph = referenceSystem.readEphemeris(string(kernel));
            [ep, r, q] = tc.sampleState();
            ep5 = referenceSystem.daySecondToEpoch([ep.day(1:5) ep.sec(1:5)], "GPS"); r5 = r(:, 1:5); q5 = q(:, 1:5);
            [srcD, infoD] = NonGravForces.env.radiationSources(ep5, r5, q5, ...
                NonGravForces.Options("standard", "sun.source", "rs.de440"), 'Ephemeris', eph);
            [srcA, infoA] = NonGravForces.env.radiationSources(ep5, r5, q5, NonGravForces.Options("standard"));
            tc.verifyEqual(infoD.sunSource, "rs.de440", 'the source list must report which ephemeris produced it');
            tc.verifyEqual(infoA.sunSource, "rs.analytic");
            % Direction, in GCRS so the comparison does not run through the attitude.
            d = infoD.sunUnit_gcrs; a = infoA.sunUnit_gcrs;
            arcsec = max(atan2(vecnorm(cross(d, a)), dot(d, a))) * 180 / pi * 3600;
            tc.verifyBoundAndRecord("DE440 vs the analytic Sun: direction", arcsec, 0.05, "arcsec", ...
                "measured 0.0031 arcsec 2026-09-16; 1 arcsec of Sun direction moves a panel's cosine by 5e-6");
            tc.verifyBoundAndRecord("DE440 vs the analytic Sun: distance", ...
                max(abs(infoD.sunDistance_m ./ infoA.sunDistance_m - 1)), 1e-7, "-", ...
                "measured 2.0e-9 (0.3 km) 2026-09-16");
            % What the toolbox actually emits: the solar radiation pressure of source 1.
            tc.verifyBoundAndRecord("DE440 vs the analytic Sun: solar radiation pressure", ...
                max(abs(srcD.P(1, :) ./ srcA.P(1, :) - 1)), 1e-7, "-", ...
                "twice the relative distance difference (1/r^2); measured 4.0e-9, far below the 1e-3 uncertainty of the solar constant itself");
        end
        function radiationSourcesContract(tc)
            [ep, r, q] = tc.sampleState();
            ep5 = referenceSystem.daySecondToEpoch([ep.day(1:5) ep.sec(1:5)], "GPS"); r5 = r(:, 1:5); q5 = q(:, 1:5);
            [src, info] = NonGravForces.env.radiationSources(ep5, r5, q5, NonGravForces.Options("standard"));
            tc.verifyEqual(src.M, 1 + 2 * 288); tc.verifySize(src.dir, [3 577 5]); tc.verifySize(src.P, [577 5]); tc.verifySize(src.band, [1 577]);
            tc.verifyEqual(src.band(1), "vis"); tc.verifyEqual(nnz(src.band == "vis"), 289); tc.verifyEqual(nnz(src.band == "ir"), 288);
            tc.verifySize(src.nu, [1 5]); tc.verifyEqual(src.frame, "body"); tc.verifyTrue(all(src.P(:) >= 0));
            tc.verifyLessThan(max(abs(vecnorm(src.dir(:, 1, :)) - 1)), 1e-14, 'unit Sun direction');
            tc.verifyEqual(info.sunSource, "rs.analytic"); tc.verifyEqual(info.erpBackend, "knocke"); tc.verifyEqual(info.provenance.unit, "env.radiationSources");
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep5, r5(:, 1:4), q5, NonGravForces.Options("standard")), 'NonGravForces:radiationSources:inputShape');
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep5, r5, q5(1:3, :), NonGravForces.Options("standard")), 'NonGravForces:radiationSources:inputShape');
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep5, r5, q5, NonGravForces.Options("standard", "sun.source", "rs.de440")), 'NonGravForces:radiationSources:depKernel');
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep5, r5, q5, NonGravForces.Options("standard"), 'Sun', [tc.au; 0; 0]), 'NonGravForces:radiationSources:inputShape');
            tc.verifyError(@() NonGravForces.env.radiationSources(ep5, r5, q5, struct()), 'MATLAB:validation:UnableToConvert');
            % ceres: 1 + 2 x tiles sources as for knocke; the CERES inputs are refused where nothing reads them
            oC = NonGravForces.Options("standard", "erp.backend", "ceres");
            sC = NonGravForces.env.radiationSources(ep5, r5, q5, oC, 'Ceres', tc.cer, 'RotationItrsToGcrs', tc.eyes(5));
            tc.verifyEqual(sC.M, 1 + 2 * 288); tc.verifyTrue(all(sC.P(:) >= 0)); tc.verifyEqual(sC.meta.erpBackend, "ceres");
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep5, r5, q5, oC), 'NonGravForces:earthRadiance:inputCeres');
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep5, r5, q5, NonGravForces.Options("standard", "erp.backend", "none"), 'Ceres', tc.cer), 'NonGravForces:radiationSources:inputCeres');
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep5, r5, q5, NonGravForces.Options("standard"), 'Ceres', tc.cer, 'RotationItrsToGcrs', tc.eyes(5)), 'NonGravForces:earthRadiance:inputCeres');
        end
        function buildMexContract(tc)
            % 'Target' names one gateway and the report is scalar; the default "all" returns one
            % element per target and is asserted at the end of this method.
            rep = NonGravForces.buildMex('Check', true, 'Target', "msis");
            tc.verifyTrue(rep.present); tc.verifyTrue(rep.upToDate); tc.verifyFalse(rep.built); tc.verifyEqual(rep.compiler, "");
            tc.verifyEqual(strlength(rep.sourceDigest), 64); tc.verifyEqual(rep.sourceDigest, rep.binaryDigest);
            tc.verifyTrue(all(isfile(rep.sources))); tc.verifyEqual(numel(rep.sources), 4);
            tc.verifyTrue(endsWith(rep.file, "msisGateway." + mexext));
            tc.verifyEqual(rep.target, "msis");
            all2 = NonGravForces.buildMex('Check', true);
            tc.verifyEqual([all2.target], ["msis" "msis21" "hwm14" "dtm2020" "rhs" "kernel"]);
            rh = all2(3);
            % 7 = the gateway, the shared chdir-guarded loader ngfHwm14Load.c (which the rhs target
            % links too), NRL's hwm14.f90, this toolbox's shim, and the three data files.
            tc.verifyTrue(rh.upToDate); tc.verifyEqual(numel(rh.sources), 7); tc.verifyTrue(endsWith(rh.file, "hwm14Gateway." + mexext));
            tc.verifyTrue(all(isfile(fullfile(fileparts(rh.file), ["hwm123114.bin" "dwm07b104i.dat" "gd2qd.dat"]))), 'the three data files are copied beside the gateway');
            % msis21: the gateway, NRL's seven Fortran files, the shim, the parameter file and the licence, all digested
            r21 = all2(2);
            tc.verifyTrue(r21.upToDate); tc.verifyEqual(numel(r21.sources), 11); tc.verifyTrue(endsWith(r21.file, "msis21Gateway." + mexext));
            tc.verifyEqual(nnz(endsWith(r21.sources, ".F90")), 8); tc.verifyTrue(any(endsWith(r21.sources, "msis21.parm")));
            tc.verifyTrue(all(isfile(fullfile(fileparts(r21.file), ["msis21.parm" "nrlmsis2.1_license..txt"]))), 'the parameter file and NRL''s licence are copied beside the gateway');
            % rhs (LP-C WP-C9): the one target that links everything -- its own C, the kernel's,
            % NRLMSISE-00's, and the NRLMSIS 2.1 and HWM14 Fortran -- and the one that compiles
            % against a source it does not own, gjToolbox's csrc/gj_native.h.
            % dtm2020 (2026-10-27): the gateway, CNES's subroutine file, the shim, the coefficient file and
            % CNES's licence, all digested; the coefficients and the licence travel beside the binary.
            rd = all2(4);
            tc.verifyTrue(rd.upToDate); tc.verifyEqual(numel(rd.sources), 5); tc.verifyTrue(endsWith(rd.file, "dtm2020Gateway." + mexext));
            tc.verifyTrue(all(isfile(fullfile(fileparts(rd.file), ["DTM_2020_F107_Kp.dat" "LICENSE_MCM_DTM.txt"]))), 'the coefficient file and CNES''s licence are copied beside the gateway');
            rr = all2(5);
            tc.verifyTrue(rr.upToDate); tc.verifyTrue(endsWith(rr.file, "rhsGateway." + mexext));
            tc.verifyEqual(numel(rr.sources), 19);                    % 9 C/H + 8 msis21 F90 + 2 hwm14 F90
            tc.verifyTrue(any(endsWith(rr.sources, "gj_native.h")), 'the frozen ABI header takes part in the digest, so a new ABI rebuilds this target');
            tc.verifyTrue(endsWith(rr.sourceDir, "csrc"), 'the rhs target''s third-party source is odeSuite''s csrc');
            % ... and it must resolve WITHOUT the environment variable, which is the branch CI
            % depends on and the only one nothing had ever exercised: every run on the owner's host
            % set NGF_GJTOOLBOX, so a stray '..' in the folder walk sent buildMex to the parent of
            % odeSuite and four green gates in a row said nothing (LP-C, 2026-09-16).
            saved = getenv('NGF_GJTOOLBOX');
            restore = onCleanup(@() setenv('NGF_GJTOOLBOX', saved));
            setenv('NGF_GJTOOLBOX', '');
            fromPath = NonGravForces.buildMex('Check', true, 'Target', "rhs");
            tc.verifyEqual(fromPath.sourceDir, rr.sourceDir, 'the path walk must reach the same csrc the environment variable names');
            tc.verifyTrue(fromPath.sourcesPresent, 'gj_native.h must be found through the gjToolbox on the path alone');
            clear restore
            tc.verifyEqual(all2(6).sourceDir, "", 'the kernel has no third-party source');
            tc.verifyError(@() NonGravForces.buildMex('Check', true, 'Target', "nosuch"), 'MATLAB:validators:mustBeMember');
        end
    end

    % ===================================================================== reference
    methods (Test, TestTags = {'reference'})
        function windReference(tc)
            % An independent route to the same number: the nominal rigid rotation omega x r,
            % with omega from referenceSystem's own IERS constant and the position in ITRS.
            % The two agree to the polar motion and length-of-day terms the nominal rate omits.
            [ep, r, ~, eop, v] = tc.sampleState();
            [~, info] = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
            R = referenceSystem.rotationItrsToGcrs(ep, 'Eop', eop);
            K = size(r, 2);
            rI = reshape(pagemtimes(pagetranspose(R), reshape(r, 3, 1, K)), 3, K);
            om = referenceSystem.iersConstants().earthRotationRate;
            nominalI = [-om * rI(2, :); om * rI(1, :); zeros(1, K)];
            nominalG = reshape(pagemtimes(R, reshape(nominalI, 3, 1, K)), 3, K);
            rel = max(vecnorm(info.corotation_m_s - nominalG)) / max(vecnorm(nominalG));
            tc.verifyBoundAndRecord("co-rotation vs a nominal omega x r", rel, 5e-3, "-", "Rdot carries polar motion and the length-of-day term that a nominal rate omits");
            % The magnitude must follow omega * sqrt(x^2 + y^2) in ITRS at EVERY epoch, which is the
            % rigid-rotation law itself. The arc spans -59.5 to -21.6 deg latitude, so the speed runs
            % from 255.8 to 467.3 m/s: asserting one number would only have described its first epoch.
            nominalSpeed = om * hypot(rI(1, :), rI(2, :));
            tc.verifyBoundAndRecord("co-rotation speed vs omega * sqrt(x^2 + y^2) over the arc", ...
                max(abs(vecnorm(info.corotation_m_s) - nominalSpeed) ./ nominalSpeed), 5e-3, "-", ...
                "the residual is the polar motion and the length-of-day term; measured 2e-6 over 600 epochs spanning 256 to 467 m/s");
        end
        function windHwm14Reference(tc)
            % Oracle: NRL's own check output for HWM14 (Check/gfortran.txt, read from the package the MEX was built from),
            % printed to 1 mm/s, for the height, latitude, local-time, longitude, day-of-year and ap profiles of
            % checkhwm14.f90, whose fixed inputs are transcribed below from that driver. Through THIS unit: epoch,
            % GCRS position, rotation, geodetic conversion and the rotation of the wind back are all in the loop.
            % Quiet columns use 'Disturbed', false; total columns 'Ap'. Rows the unit cannot pose are left out and
            % named: the two poles (north and east undefined by a position), day 0 (no such date), local times
            % whose UT the driver pushes past 86400 s (a date interface normalises them into the next day), and the
            % surface row (a point on the ellipsoid comes back from the geodetic round trip a nanometre below it,
            % where the unit refuses because HWM14 does not return).
            r21 = NonGravForces.buildMex('Check', true, 'Target', "hwm14");
            txt = splitlines(string(fileread(fullfile(r21.sourceDir, 'gfortran.txt'))));
            eop = referenceSystem.readEop(fullfile(referenceSystem.dataFolder(), 'finals2000A.all'));
            prof = struct('title', {"height profile", "latitude profile", "local time profile", "longitude profile", "day of year profile", "magnetic activity profile"});
            worst = 0; rows = 0;
            for p = 1:numel(prof)
                k0 = find(strtrim(txt) == prof(p).title, 1) + 4;
                tab = [];
                while k0 <= numel(txt) && strlength(strtrim(txt(k0))) > 0
                    tab(end+1, :) = double(split(strtrim(txt(k0))))'; %#ok<AGROW>
                    k0 = k0 + 1;
                end
                x = tab(:, 1)'; n = numel(x);
                switch p
                    case 1, doy = 150 * ones(1, n); sec = 43200 * ones(1, n); la = -45 * ones(1, n); lo = -85 * ones(1, n); h = x; ap = 80 * ones(1, n); keep = x > 0;
                    case 2, doy = 305 * ones(1, n); sec = 64800 * ones(1, n); la = x; lo = 30 * ones(1, n); h = 250 * ones(1, n); ap = 48 * ones(1, n); keep = abs(x) < 90;
                    case 3, doy = 75 * ones(1, n); lo = -70 * ones(1, n); sec = double(single((single(x) - single(lo) / 15) * 3600)); la = 45 * ones(1, n); h = 125 * ones(1, n); ap = 30 * ones(1, n); keep = sec < 86400;
                    case 4, doy = 330 * ones(1, n); sec = 21600 * ones(1, n); la = -5 * ones(1, n); lo = x; h = 40 * ones(1, n); ap = 4 * ones(1, n); keep = true(1, n);
                    case 5, doy = x; sec = 75600 * ones(1, n); la = -65 * ones(1, n); lo = -135 * ones(1, n); h = 200 * ones(1, n); ap = 15 * ones(1, n); keep = x > 0;
                    case 6, doy = 280 * ones(1, n); sec = 75600 * ones(1, n); la = 38 * ones(1, n); lo = 125 * ones(1, n); h = 350 * ones(1, n); ap = x; keep = true(1, n);
                end
                for i = find(keep)
                    iq = tc.hwmAt(eop, doy(i), sec(i), la(i), lo(i), h(i), 'Disturbed', false);
                    it = tc.hwmAt(eop, doy(i), sec(i), la(i), lo(i), h(i), 'Ap', ap(i));
                    q = iq.hwm14.northEast_m_s; t = it.hwm14.northEast_m_s;
                    worst = max([worst, abs(q' - tab(i, 2:3)), abs(t' - tab(i, 6:7)), abs((t - q)' - tab(i, 4:5)) / 2]);
                    rows = rows + 1;
                end
            end
            tc.assertGreaterThan(rows, 90, 'the six profiles must be exercised');
            tc.verifyBoundAndRecord("hwm14 through env.wind vs NRL Check/gfortran.txt (" + string(rows) + " rows, m/s)", worst, 6e-4, "m/s", ...
                "half a printed unit (5e-4) plus the model's single-precision response to a geodetic round trip; disturbed = total - quiet carries two roundings, halved");
        end
        function densityDtm2020Reference(tc)
            % CNES's benchmark, every printed number: tz and tinf to 3 decimals, ro to 5 significant digits.
            [ep, lla, drv, lst, ref] = tc.dtmBenchmark();
            [rho, ~, T] = NonGravForces.env.density("dtm2020", ep, lla, [], 'Drivers', drv, 'LocalSolarTime_h', lst);
            tc.verifyBoundAndRecord("dtm2020 vs CNES benchmark: tz (4 cases)", max(abs(T(2, :) - ref.tz)), 1e-3, "K", ...
                "printed to 3 decimals (half a unit 5e-4); first run 2026-10-27");
            tc.verifyBoundAndRecord("dtm2020 vs CNES benchmark: tinf (4 cases)", max(abs(T(1, :) - ref.tinf)), 1e-3, "K", ...
                "printed to 3 decimals (half a unit 5e-4); first run 2026-10-27");
            tc.verifyBoundAndRecord("dtm2020 vs CNES benchmark: ro, relative (4 cases)", max(abs(rho / 1e3 - ref.ro) ./ ref.ro), 1e-4, "-", ...
                "printed to 5 significant digits (half a unit <= 1.9e-5 relative); first run 2026-10-27");
        end
        function densityReference(tc)
            % The 17 published cases of DOCUMENTATION §7 (printed with %E, 7 significant digits)
            % must be reproduced digit for digit in every column: HE O N2 O2 AR RHO H N ANOM TINF TG.
            [rhoA, spA, TA] = tc.msisCases(1:15, false);
            [rhoB, spB, TB] = tc.msisCases(16:17, false);
            rho = [rhoA rhoB]; sp = [spA spB]; T = [TA TB];
            mine = [sp(1:5, :) / 1e6; rho / 1e3; sp(6:8, :) / 1e6; T];          % back to the reference's cgs
            worst = 0; mismatches = strings(0, 1);
            for i = 1:17
                for j = 1:11
                    a = sprintf('%.6E', mine(j, i)); b = sprintf('%.6E', tc.msisTable(i, j));
                    if ~strcmp(a, b), mismatches(end+1) = sprintf("case %d col %d: %s vs %s", i - 1, j, a, b); end %#ok<AGROW>
                    worst = max(worst, abs(mine(j, i) - tc.msisTable(i, j)) / max(abs(tc.msisTable(i, j)), 1e-300));
                end
            end
            tc.verifyEmpty(mismatches, strjoin(mismatches, newline));
            tc.verifyBoundAndRecord("msis00 vs DOCUMENTATION table (17 cases)", worst, 5e-7, "rel", "printed 7 digits; measured 4.8e-7");
            % NRLMSIS 2.1: NRL's 200 cases. The reference prints sngl() of the double-precision model
            % with E13.4 (densities, cgs) and F8.2 (T), and the package states that -DDBLE reproduces it
            % exactly on any compiler. The claim is therefore "within half a unit of the last printed
            % digit", widened only by the single-precision cast that happened before the printing.
            % NO (column 10) is not returned by this toolbox; missing species (9.999e-38) must be 0.
            [rho21, sp21, T21] = tc.msis21Cases(1:200);
            mine21 = [sp21(1:5, :) * 1e-6; rho21 * 1e-3; sp21(6:8, :) * 1e-6; T21(2, :)]';   % 200x10: He O N2 O2 Ar rho H N O* T
            ref21 = tc.msis21Ref(:, [1:9 11]);
            missing = ref21 < 1e-36;                                           % printed 0.9999E-37
            tc.verifyEqual(mine21(missing), zeros(nnz(missing), 1), 'a species NRLMSIS does not model there is returned as 0');
            tc.verifyGreaterThan(nnz(missing), 100, 'the oracle must exercise the missing-species rule (measured 103 rows print it)');
            e = floor(log10(abs(ref21))) + 1;                                  % 0.dddd E e
            halfUnit = 0.5 * 10 .^ (e - 4);
            halfUnit(:, 10) = 0.005;                                           % F8.2
            allowed = halfUnit + double(eps('single')) * abs(ref21);
            ratio = abs(mine21 - ref21) ./ allowed;
            ratio(missing) = 0;
            tc.verifyBoundAndRecord("msis21 vs NRL msis2.1_test_ref_dp (200 cases, 10 columns), in half last printed digits", max(ratio, [], 'all'), 1, "-", "the package's own statement: double precision reproduces its reference exactly; rounding of the print only");
        end
        function radiationSourcesCeresReference(tc)
            % Oracle: tools/pyxcheck/ceres_check.py integrates the SAME continuous CERES field (geodetic bilinear in space,
            % linear in time, albedo = SWup/SWdown on the instantaneous insolation, longwave as given, on the 20 km
            % reference sphere) over the visible cap by Gauss-Legendre quadrature to < 2e-5, written from a specification
            % by an agent that did not open this toolbox (T9). Eight cases: sub-solar, terminator, polar summer and night,
            % over the pole, the longitude seam, a node time, the night side.
            % THE CLAIM IS IDENTITY OF THE MODEL, and the tile sum converges to it at second order in the ring width. So
            % the Richardson extrapolation of [24 48] and [48 96] must meet the oracle to its own quadrature error, and
            % the ratio of the two tiling errors must be 4; a model that differed would leave a floor and break the ratio.
            T = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'ceres_oracle.csv'), 'TextType', 'string');
            err = zeros(height(T), 4, 2);                        % case x [albSum lwSum albVec lwVec] x grid
            for g = 1:2
                grid = [24 48] * g;
                o = NonGravForces.Options("standard", "erp.backend", "ceres", "erp.grid", grid, "srp.eclipse", "cylindrical");
                for i = 1:height(T)
                    ep = referenceSystem.datetimeToEpoch(datetime(T.epoch_utc(i), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''), "UTC");
                    r = [T.rx_m(i); T.ry_m(i); T.rz_m(i)];
                    rs = [T.sx(i); T.sy(i); T.sz(i)] * tc.au * sqrt(o.srp.tsi / T.S_W_m2(i));   % so that tsi (au/|rs|)^2 = S of the case
                    src = NonGravForces.env.radiationSources(ep, r, [1; 0; 0; 0], o, 'Sun', rs, 'Ceres', tc.cer, 'RotationItrsToGcrs', eye(3));
                    Mt = (src.M - 1) / 2; ia = 2:1+Mt; il = 2+Mt:src.M;
                    scale = T.E_alb_sum(i) + T.E_lw_sum(i);      % the total Earth irradiance: albedo alone can be exactly 0
                    err(i, :, g) = [sum(src.P(ia)) * tc.cLight - T.E_alb_sum(i), sum(src.P(il)) * tc.cLight - T.E_lw_sum(i), ...
                        norm(src.dir(:, ia) * src.P(ia) * tc.cLight - [T.E_alb_x(i); T.E_alb_y(i); T.E_alb_z(i)]), ...
                        norm(src.dir(:, il) * src.P(il) * tc.cLight - [T.E_lw_x(i); T.E_lw_y(i); T.E_lw_z(i)])] / scale;
                end
            end
            tc.verifyBoundAndRecord("ceres Earth sources vs the oracle, tiling error at [24 48] (precise preset)", max(abs(err(:, :, 1)), [], 'all'), 1e-2, "-", ...
                "second-order tiling error of the cap; measured 3.5e-3 of the total Earth irradiance on 2026-09-15");
            rich = (4 * err(:, 1:2, 2) - err(:, 1:2, 1)) / 3;    % the scalar sums, extrapolated to zero ring width
            tc.verifyBoundAndRecord("ceres Earth sources, Richardson-extrapolated tile sum vs the oracle", max(abs(rich), [], 'all'), 1e-4, "-", ...
                "the oracle's own quadrature error (1.7e-5) plus the next term of the tiling series; measured 5.6e-5 on 2026-09-15");
            lw = abs(err(:, 2, :));
            tc.verifyBoundAndRecord("ceres longwave tiling error ratio [24 48] / [48 96], departure from 4", max(abs(lw(:, 1) ./ lw(:, 2) - 4)), 0.5, "-", ...
                "second-order convergence to the same integral: the ratios measured 3.80 to 4.20 on 2026-09-15");
            tc.verifyTrue(any(T.w_t == 0) && any(T.w_t == 1), 'the oracle must exercise both month nodes exactly');
        end
        function magneticFieldReference(tc)
            % Oracle 1: tools/pyxcheck/igrf_check.py, an independent standard-library synthesis written from the IGRF
            % description by an agent that did not open this toolbox (T9), for 12 cases: both poles exactly, the dipole
            % alone, the time rule inside a five-year step and in the secular-variation extrapolation, 1900, and a
            % geostationary radius. Each case takes its coefficients from io.igrfCoefficients at the case's year.
            T = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'igrf_oracle.csv'));
            worst = 0;
            for i = 1:height(T)
                c = NonGravForces.io.igrfCoefficients(fullfile(tc.root, 'data', 'igrf', 'igrf14coeffs.txt'), 'Epoch', T.year(i));
                u = [cosd(T.lat_deg(i)) * cosd(T.lon_deg(i)); cosd(T.lat_deg(i)) * sind(T.lon_deg(i)); sind(T.lat_deg(i))];
                [B, info] = NonGravForces.env.magneticField(c, T.r_km(i) * 1e3 * u, 'MaxDegree', T.N(i));
                want = [T.Bx_nT(i); T.By_nT(i); T.Bz_nT(i)] * 1e-9;
                wantNED = [T.X_nT(i); T.Y_nT(i); T.Z_nT(i)] * 1e-9;
                worst = max(worst, norm(B - want) / norm(want));
                if abs(T.lat_deg(i)) < 90
                    % AT A POLE north and east are not defined by the position: the oracle names them by the case's
                    % longitude, a position on the axis carries none. The Cartesian vector is compared there instead.
                    worst = max(worst, norm(info.northEastDown_T - wantNED) / norm(wantNED));
                end
            end
            tc.verifyBoundAndRecord("magneticField vs the independent synthesis (12 cases, Cartesian; NED off the poles)", worst, 1e-13, "-", ...
                "two routes to the same series of at most 104 terms; measured 8.7e-16 on 2026-09-15");
            % Oracle 2: the GRACE-C magnetometer (MAG1B, raw) over the first five minutes of 2020-01-01 at the GNV1B
            % positions. The model is the main field only and the instrument is uncalibrated, so the claim is the
            % intensity to a few parts per thousand; measured 0.9961 to 0.9985 on 2026-09-15.
            S = fullfile(tc.root, 'data', 'samples');
            [epG, rI] = NonGravForces.io.readGNV1B(fullfile(S, 'GNV1B_2020-01-01_C_04.txt'), 'Frame', "ITRS");
            [epM, fld] = NonGravForces.io.readMAG1B(fullfile(S, 'MAG1B_2020-01-01_C_04.txt'));
            [tf, loc] = ismember(round(epM.sec(:)'), epG.sec(:)');
            tc.assertGreaterThan(nnz(tf), 500, 'MAG1B and GNV1B excerpts must overlap');
            [~, info] = NonGravForces.env.magneticField(tc.igrf, rI(:, loc(tf)));
            q = fld.magnitude_T(tf) ./ info.intensity_T;
            tc.verifyBoundAndRecord("GRACE-C MAG1B intensity / IGRF-14 intensity, largest departure from 1", max(abs(q - 1)), 0.012, "-", ...
                "main field only against a raw magnetometer; measured 3.9e-3 on 2026-09-15, bound 3x");
        end
        function eclipseFactorReference(tc)
            % Oracle: pure-Python numerical disc-overlap integration (eclipse_check.py); cylindrical bitwise.
            E = tc.eclipse; r = [E.rx E.ry E.rz]'; s = [E.sx E.sy E.sz]';
            nuC = NonGravForces.env.eclipseFactor(r, s); nuY = NonGravForces.env.eclipseFactor(r, s, 'Model', "cylindrical");
            % HANDOVER accuracy: nu to 1e-12; the oracle's own truncation is ~1e-12. Bound 1e-11.
            tc.verifyBoundAndRecord("eclipse conical vs Python disc integration", max(abs(nuC - E.nu_conical')), 1e-11, "-", "oracle ~1e-12; atan2 form (acos form measured 6.4e-10)");
            tc.verifyEqual(nuY, E.nu_cylindrical', 'cylindrical is a pure geometric test: bitwise');
            tc.verifyGreaterThan(nnz(nuC > 0 & nuC < 1), 3, 'the oracle must exercise the penumbra');
            % SOLAARS-CF: Table 3.2 against the independent restatement, and the formula against the pure-Python
            % evaluation on a (rR, rE) grid; the satellite is displaced along the ecliptic normal, where the
            % oblateness scaling does not act, so rE is exactly the grid value.
            C = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'solaars_cf_coefficients.csv'), 'TextType', 'string');
            O = readtable(fullfile(tc.root, 'tools', 'pyxcheck', 'solaars_oracle.csv'));
            rS = [-O.rR_Mm' * 1e6; zeros(1, height(O)); O.rE_Mm' * 1e6];
            [nuS, infoS] = NonGravForces.env.eclipseFactor(rS, [tc.au; 0; 0], 'Model', "solaarsCF", 'EclipticNormal', [0; 0; 1]);
            tc.verifyEqual(infoS.solaars.rE_Mm, O.rE_Mm', 'AbsTol', 1e-12); tc.verifyEqual(infoS.solaars.rR_Mm, O.rR_Mm', 'AbsTol', 1e-12);
            tc.verifyBoundAndRecord("SOLAARS-CF vs Python restatement of eq. 3.86", max(abs(nuS - O.nu')), 1e-14, "-", "same coefficients, tanh/exp rounding");
            tc.verifyEqual(height(C), 8); tc.verifyEqual(C.b1(1), 0.1715); tc.verifyEqual(C.b4(8), -0.1318); tc.verifyEqual(C.b3(4), 6.385);
        end
        function earthRadianceReference(tc)
            % Oracle: knocke_coefficients.csv (independent restatement of Table 1) and the formula evaluated here from it.
            Lc = splitlines(strtrim(string(fileread(fullfile(tc.root, 'tools', 'pyxcheck', 'knocke_coefficients.csv'))))); Lc = Lc(2:end);
            kv = split(Lc, ","); names = strtrim(kv(:, 1)); values = strtrim(kv(:, 2));
            val = @(n) str2double(values(names == n));
            ep = referenceSystem.datetimeToEpoch(datetime(1981, 12, 22) + days([0 100 200]), "UTC"); r = repmat([6.88e6; 0; 0], 1, 3);
            t = NonGravForces.env.earthRadiance("knocke", ep, [12 24], r);
            for n = ["a0","a1","a2","e0","e1","e2","period_days"], tc.verifyEqual(t.coefficients.(n), val(n), sprintf('coefficient %s', n)); end
            tc.verifyEqual(t.coefficients.t0, values(names == "t0"));
            x = sin(t.lat_rad); seasonal = cos(2 * pi / val("period_days") * [0 100 200]);
            aExp = val("a0") + val("a1") * (x .* seasonal) + val("a2") * 0.5 * (3 * x.^2 - 1);
            eExp = val("e0") + val("e1") * (x .* seasonal) + val("e2") * 0.5 * (3 * x.^2 - 1);
            tc.verifyBoundAndRecord("Knocke albedo tile values vs Table 1 formula", max(abs(t.albedo(:) - aExp(:))), 1e-15, "-", "same arithmetic from an independent coefficient file");
            tc.verifyBoundAndRecord("Knocke emissivity tile values vs Table 1 formula", max(abs(t.emissivity(:) - eExp(:))), 1e-15, "-", "");
            tc.verifyEqual(t.daysSinceT0, [0 100 200]);
        end
        function radiationSourcesReference(tc)
            % Oracle: knocke_check.py cap integration at 0.25 deg. Errors relative to the longwave flux, per grid.
            % Measured 2026-09-07 (worst over the six zonal cases): [12 24] 2.3e-2, [24 48] 5.4e-3, [48 96] 1.4e-3;
            % first-order convergence in the ring count. Bounds 3x-ish: 4e-2, 1e-2, 3e-3.
            K = tc.knocke; zonal = find(~K.uniform)'; uni = find(K.uniform)';
            bounds = struct('g12', 4e-2, 'g24', 1e-2, 'g48', 3e-3); grids = {[12 24], [24 48], [48 96]}; names = ["g12","g24","g48"];
            for gi = 1:3
                worst = 0;
                for k = zonal
                    [Ea, Ei] = tc.knockeCase(k, grids{gi}, "knocke");
                    worst = max(worst, max(abs(Ea - K.E_albedo_W_m2(k)), abs(Ei - K.E_ir_W_m2(k))) / K.E_ir_W_m2(k));
                end
                tc.verifyBoundAndRecord("ERP cap grid " + mat2str(grids{gi}) + " vs Python integration", worst, bounds.(names(gi)), "rel", "relative to the longwave flux, worst of 6 cases");
            end
            % point source: the far-field Lambert sphere; at 100 R the finite-distance (R/r) terms are 1e-2 and the
            % oracle's rim discretisation ~2e-3. Bound 2e-2 on albedo, 5e-3 on longwave (measured: see last_run).
            far = uni(K.rx(uni).^2 + K.rz(uni).^2 > (50 * 6378137)^2);
            worstA = 0; worstI = 0;
            for k = far
                [Ea, Ei, src] = tc.knockeCase(k, [], "pointSource");
                tc.verifyEqual(src.M, 3);
                worstA = max(worstA, abs(Ea - K.E_albedo_W_m2(k)) / K.E_albedo_W_m2(k)); worstI = max(worstI, abs(Ei - K.E_ir_W_m2(k)) / K.E_ir_W_m2(k));
            end
            tc.verifyBoundAndRecord("pointSource albedo vs uniform sphere at 100 R", worstA, 2e-2, "rel", "far-field law, O(R/r)");
            tc.verifyBoundAndRecord("pointSource longwave vs uniform sphere at 100 R", worstI, 5e-3, "rel", "exact law; oracle rim tiles ~2e-3");
            % SRP: at exactly 1 au and full light P = tsi / c (CODATA c, IAU 2015 B3 tsi 1361).
            ep = referenceSystem.mjdToEpoch(58849.5, "UTC");
            src = NonGravForces.env.radiationSources(ep, [7e6; 0; 0], [1; 0; 0; 0], NonGravForces.Options("standard", "erp.backend", "none"), 'Sun', [tc.au + 7e6; 0; 0]);
            tc.verifyBoundAndRecord("SRP at 1 au vs tsi/c", abs(src.P(1) - 1361 / 299792458) / (1361 / 299792458), 1e-15, "rel", "one division");
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function windPrecision(tc)
            % The corotation branch subtracts exactly Rdot*(R'*r) and nothing else.
            [ep, r, ~, eop, v] = tc.sampleState();
            [vrel, info] = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
            [R, Rdot] = referenceSystem.rotationItrsToGcrs(ep, 'Eop', eop);
            K = size(r, 2);
            rI = reshape(pagemtimes(pagetranspose(R), reshape(r, 3, 1, K)), 3, K);
            want = reshape(pagemtimes(Rdot, reshape(rI, 3, 1, K)), 3, K);
            tc.verifyBoundAndRecord("co-rotation term vs Rdot*(R'*r)", max(abs(info.corotation_m_s - want), [], 'all'), 0, "m/s", "the same two products in the same order: bitwise");
            tc.verifyEqual(vrel, v - want, 'and the relative velocity is that difference, bitwise');
            tc.verifyEqual(info.speed_m_s, vecnorm(vrel));
        end
        function windHwm14Precision(tc)
            % The model wind is HORIZONTAL at the geodetic point: its GCRS vector, rotated back to the ITRS, has no
            % component along the ellipsoid normal and reproduces the north and east components exactly.
            [ep, r, ~, eop, v] = tc.sampleState();
            [~, info] = NonGravForces.env.wind("hwm14", ep, r, v, 'Eop', eop, 'Disturbed', false);
            R = referenceSystem.rotationItrsToGcrs(ep, 'Eop', eop);
            K = size(r, 2);
            wI = reshape(pagemtimes(pagetranspose(R), reshape(info.hwm14.gcrs_m_s, 3, 1, K)), 3, K);
            la = info.hwm14.geodetic(1, :); lo = info.hwm14.geodetic(2, :);
            up = [cos(la) .* cos(lo); cos(la) .* sin(lo); sin(la)];
            north = [-sin(la) .* cos(lo); -sin(la) .* sin(lo); cos(la)]; east = [-sin(lo); cos(lo); zeros(1, K)];
            tc.verifyBoundAndRecord("hwm14 wind along the ellipsoid normal", max(abs(sum(wI .* up, 1))), 1e-12, "m/s", "two rotations of a 100 m/s vector");
            tc.verifyBoundAndRecord("hwm14 north and east recovered from the GCRS vector", ...
                max(abs([sum(wI .* north, 1) - info.hwm14.northEast_m_s(1, :), sum(wI .* east, 1) - info.hwm14.northEast_m_s(2, :)])), 1e-12, "m/s", "the same two rotations");
        end
        function densityDtm2020Precision(tc)
            % DTM3 forms ro as the plain sum of its six partial densities; env.density converts each with the
            % toolbox's molar masses, so n_j m_j summed must return rho to a few ulps.
            [ep, lla, drv, lst] = tc.dtmBenchmark();
            [rho, sp, ~] = NonGravForces.env.density("dtm2020", ep, lla, [], 'Drivers', drv, 'LocalSolarTime_h', lst);
            c = NonGravForces.const();
            back = sum(sp .* (c.speciesMass(:) * 1e3 * c.amu), 1);
            tc.verifyBoundAndRecord("dtm2020: sum of n_j m_j over rho - 1", max(abs(back ./ rho - 1)), 8 * eps, "-", ...
                "six conversions and a sum against the model's own sum: rounding only");
        end
        function densitySpeciesSumPrecision(tc)
            % kernel.dragPanelForce re-forms the mass density from the species and const.speciesMass,
            % while NRLMSISE-00 forms its own total from ITS molecular masses. Its ACCURACY block
            % quotes the difference, so the difference has to be measured here rather than asserted.
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            ep = referenceSystem.mjdToEpoch([58849.2 58849.5 58849.8], "UTC");
            lla = [0 0.5 -0.7; 0 1.0 2.0; 490e3 500e3 480e3];
            [rho, sp] = NonGravForces.env.density("msis00", ep, lla, sw);
            c = NonGravForces.const();
            rel = abs(sum(sp .* (c.speciesMass(:) / c.NA), 1) ./ rho - 1);
            tc.verifyBoundAndRecord("msis00 total vs the species sum over const.speciesMass", max(rel), 1e-3, "-", "the backend's molecular masses differ from CODATA/IUPAC values; measured 5.7e-4 on 2026-09-07 and quoted by kernel.dragPanelForce");
            tc.verifyGreaterThan(max(rel), 1e-5, 'if this ever drops the ACCURACY block of dragPanelForce is stale in the other direction');
            % NRLMSIS 2.1 forms its total as a dot product with ITS masses (msis_constants.F90, CIPM 2007 molar masses
            % over 1e3 N_A, NO weighted 0). With those masses the species must re-form the total to rounding.
            [rho21, sp21] = NonGravForces.env.density("msis21", ep, lla, sw);
            NA = 6.02214076e23;
            m21 = [4.0 31.9988/2 28.0134 31.9988 39.948 1.0 28.0134/2 31.9988/2]' / (1e3 * NA);     % He O N2 O2 Ar H N O*
            tc.verifyBoundAndRecord("msis21 total vs the species sum over the model's own masses", max(abs(sum(sp21 .* m21, 1) ./ rho21 - 1)), 8 * eps, "-", "eight products summed in another order");
            rel21 = abs(sum(sp21 .* (c.speciesMass(:) / c.NA), 1) ./ rho21 - 1);
            tc.verifyBoundAndRecord("msis21 total vs the species sum over const.speciesMass", max(rel21), 1e-3, "-", "the model's molar masses against this toolbox's; the kernel re-forms the density from the latter");
        end
        function densityPrecision(tc)
            % Exponential model: rho(h0) = rho0 exactly; rho(h0 + H) = rho0/e to 4 eps (one exp, one multiply).
            ep = referenceSystem.mjdToEpoch([58849.5; 58849.5], "UTC"); p = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);
            rho = NonGravForces.env.density("exponential", ep, [0 0; 0 0; 450e3 510e3], [], 'Exponential', p);
            tc.verifyEqual(rho(1), 3e-12);
            tc.verifyBoundAndRecord("exponential density one scale height", abs(rho(2) - 3e-12 / exp(1)) / 3e-12, 4 * eps, "rel", "exp and a multiply");
        end
        function magneticFieldPrecision(tc)
            % Degree 1 IS a dipole: B = mu0/(4 pi) [3 (m . rhat) rhat - m] / r^3 with the moment vector
            % m = (4 pi a^3 / mu0) [g11; h11; g10], so B = (a/r)^3 [3 (d . rhat) rhat - d], d = [g11; h11; g10].
            g = tc.igrf.g * 1e-9; h = tc.igrf.h * 1e-9;
            d = [g(2, 2); h(2, 2); g(2, 1)];
            r = [[6861.2e3; 0; 0] [1e6; -3e6; 6e6] [0; 0; -7e6] [42164e3; 1; 2]];
            B = NonGravForces.env.magneticField(tc.igrf, r, 'MaxDegree', 1);
            rh = r ./ vecnorm(r);
            want = (6371.2e3 ./ vecnorm(r)).^3 .* (3 * (d' * rh) .* rh - d);
            tc.verifyBoundAndRecord("magneticField degree 1 vs the closed-form dipole", max(vecnorm(B - want) ./ vecnorm(want)), 1e-14, "-", ...
                "the same three coefficients through two formulas, one of them the exact pole");
        end
        function eclipseFactorPrecision(tc)
            % Exact limits of the conical model: full light 1, umbra 0, annular 1 - b^2/a^2 (Sun larger than the body).
            s = [tc.au; 0; 0]; R = referenceSystem.ellipsoidParameters("WGS84").a; Rs = NonGravForces.const().Rsun;
            rFar = [-3e9; 0; 0];                                     % far behind the Earth: Earth's disc smaller than the Sun's
            nu = NonGravForces.env.eclipseFactor([[7e6; 0; 0] [-6.9e6; 0; 0] rFar], s);
            tc.verifyEqual(nu(1), 1); tc.verifyEqual(nu(2), 0);
            a = asin(Rs / (tc.au + 3e9)); b = asin(R / 3e9);
            tc.verifyBoundAndRecord("eclipse annular limit 1 - b^2/a^2", abs(nu(3) - (1 - b^2 / a^2)), 1e-15, "-", "same asin, one division");
            % SOLAARS-CF limits: deep umbra (rE 6.0 Mm) and clear of the atmosphere (rE 7.0 Mm) at rR 2.56 Mm (490 km orbit):
            % the sigmoids saturate, so nu and 1 - nu are below 1e-12 there; the day side is exactly 1.
            rL = [-2.56e6 -2.56e6 +2.56e6; 0 0 0; 6.0e6 7.0e6 6.4e6];
            nuL = NonGravForces.env.eclipseFactor(rL, s, 'Model', "solaarsCF", 'EclipticNormal', [0; 0; 1]);
            tc.verifyBoundAndRecord("SOLAARS-CF umbra limit", nuL(1), 1e-12, "-", "three saturated sigmoids");
            tc.verifyBoundAndRecord("SOLAARS-CF full-light limit", 1 - nuL(2), 1e-12, "-", "");
            tc.verifyEqual(nuL(3), 1, 'day side');
        end
        function earthRadianceClimatologyPrecision(tc)
            % A field LINEAR in the grid's own latitude is reproduced exactly by bilinear
            % interpolation, so every tile must carry the closed form at the GEOCENTRIC latitude it
            % was read at. The ceres test above forms the GEODETIC latitude for the same reason and
            % the two differ by up to 0.19 deg: if this unit ever converted, this bound would break
            % by four orders of magnitude rather than subtly.
            clm = tc.syntheticClimatology(@(LAT, LON, m) 0.3 + 0.002 * LAT + 0.01 * m, ...
                                          @(LAT, LON, m) 220 + 0 * LAT, 180, 360, [58800 58830 58860]);
            ep = referenceSystem.mjdToEpoch([58810; 58840], "UTC");
            r = [[6.87e6; 1e5; 2e6] [0; 1e3; 6.9e6]];
            t = NonGravForces.env.earthRadiance("ceresClimatology", ep, [12 24], r, 'Climatology', clm, 'RotationItrsToGcrs', tc.eyes(2));
            n = t.n;
            latC = atan2d(squeeze(n(3, :, :)), squeeze(hypot(n(1, :, :), n(2, :, :))));
            lim = max(clm.albedo.lat_deg);
            want = 0.3 + 0.002 * min(max(latC, -lim), lim) + 0.01 * [1 2];      % PIECEWISE constant in time
            tc.verifyBoundAndRecord("ceresClimatology albedo of a latitude-linear field vs its closed form", max(abs(t.albedo - want), [], 'all'), 1e-13, "-", "bilinear reproduces a linear field; rounding only");
            tc.verifyBoundAndRecord("ceresClimatology longwave of a constant field", max(abs(t.exitanceLw_W_m2 - 220), [], 'all'), 1e-12, "W/m^2", "every node equal");
            tc.verifyBoundAndRecord("ceresClimatology tile latitude is GEOCENTRIC", max(abs(rad2deg(t.lat_rad) - latC), [], 'all'), 1e-12, "deg", "the grid's own convention, not converted");
            % the reference level is READ from the file, not assumed
            tc.verifyBoundAndRecord("ceresClimatology reference height measured from the grid", abs(t.climatology.referenceHeight_m - 30000.1), 1e-6, "m", "the constant the synthetic radii were built with");
        end
        function earthRadianceCeresPrecision(tc)
            % A field LINEAR in latitude is reproduced exactly by bilinear interpolation, and one linear in time by the
            % time blend, so every tile must carry the closed form at the geodetic latitude it was read at -- which this
            % test forms itself from the tile normal with tan(phi_d) = tan(phi_c)/(1 - f)^2 (poleward of +-89.5 deg the
            % grid is held constant, and so is the closed form).
            ell = referenceSystem.ellipsoidParameters("WGS84");
            c = tc.syntheticCeres(@(LAT, LON, m) deal(100 + 0.5 * LAT + 10 * m, 200 + LAT .* 0 + 0 * LON, 400 + 0 * LAT));
            ep = referenceSystem.mjdToEpoch([58840; 58850.5], "UTC");
            r = [[6.87e6; 1e5; 2e6] [0; 1e3; 6.9e6]];
            t = NonGravForces.env.earthRadiance("ceres", ep, [12 24], r, 'Ceres', c, 'RotationItrsToGcrs', tc.eyes(2));
            n = t.n;
            latD = atan2d(squeeze(n(3, :, :)), (1 - ell.f)^2 * squeeze(hypot(n(1, :, :), n(2, :, :))));
            wt = ([58840 58850.5] - c.mjd(1)) / (c.mjd(2) - c.mjd(1));
            want = (100 + 0.5 * min(max(latD, -89.5), 89.5) + 10 * (1 + wt)) / 400;
            tc.verifyBoundAndRecord("ceres albedo of a latitude-linear field vs its closed form", max(abs(t.albedo - want), [], 'all'), 1e-13, "-", "bilinear and linear blends reproduce a linear field; rounding of single-precision nodes exactly representable");
            tc.verifyBoundAndRecord("ceres longwave of a constant field", max(abs(t.exitanceLw_W_m2 - 200), [], 'all'), 1e-12, "W/m^2", "a blend of equal values");
            tc.verifyBoundAndRecord("ceres tile latitude vs the geodetic formula", max(abs(rad2deg(t.lat_rad) - latD), [], 'all'), 1e-12, "deg", "the same formula");
        end
        function earthRadiancePrecision(tc)
            % Cap area: sum of the tiles equals 2 pi R^2 (1 - cos alpha_max) (measured 4e-16); tile normals unit.
            ep = referenceSystem.mjdToEpoch([58849.5; 58850.5], "UTC"); r = [[6.88e6; 0; 0] [0; 0; 7e6]];
            t = NonGravForces.env.earthRadiance("knocke", ep, [12 24], r);
            capArea = 2 * pi * t.radius_m^2 * (1 - cos(t.capAngle_rad));
            tc.verifyBoundAndRecord("ERP cap tiles area closure", max(abs(sum(t.area_m2, 1) ./ capArea - 1)), 1e-14, "rel", "measured 4e-16");
            tc.verifyBoundAndRecord("ERP tile normals unit norm", max(abs(vecnorm(t.n) - 1), [], 'all'), 1e-14, "-", "");
        end
        function radiationSourcesPrecision(tc)
            % Anti-solar geometry: every albedo tile is unlit, so the albedo pressure is exactly zero while longwave is not;
            % the Sun direction in the body frame with the identity attitude is the GCRS unit vector to 1e-15.
            ep = referenceSystem.mjdToEpoch(58849.5, "UTC"); r = [-6.88e6; 0; 0]; sun = [tc.au; 0; 0];
            src = NonGravForces.env.radiationSources(ep, r, [1; 0; 0; 0], NonGravForces.Options("standard", "srp.eclipse", "cylindrical"), 'Sun', sun);
            Mt = (src.M - 1) / 2;
            tc.verifyEqual(sum(src.P(2:1+Mt)), 0); tc.verifyGreaterThan(sum(src.P(2+Mt:end)), 0); tc.verifyEqual(src.nu, 0);
            d = (sun - r) / norm(sun - r);
            tc.verifyBoundAndRecord("Sun direction body vs GCRS for identity q", max(abs(src.dir(:, 1) - d)), 1e-15, "-", "identity rotation");
        end
    end

    % ===================================================================== vectorization
    methods (Test, TestTags = {'vectorization'})
        function windVectorization(tc)
            [ep, r, ~, eop, v] = tc.sampleState();
            vrel = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
            for k = [1 3 numel(ep.sec)]
                e1 = referenceSystem.daySecondToEpoch([ep.day(k) ep.sec(k)], ep.scale);
                v1 = NonGravForces.env.wind("corotation", e1, r(:, k), v(:, k), 'Eop', eop);
                tc.verifyEqual(vrel(:, k), v1, sprintf('epoch %d batched vs alone', k));
            end
        end
        function windHwm14Vectorization(tc)
            % Batched equals alone, bitwise: the gateway loops per point and the model keeps no state that changes a value.
            [ep, r, ~, eop, v] = tc.sampleState();
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            idx = [1 7 300 600];
            e4 = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], ep.scale);
            [v4, i4] = NonGravForces.env.wind("hwm14", e4, r(:, idx), v(:, idx), 'Eop', eop, 'SpaceWeather', sw);
            for k = 1:4
                e1 = referenceSystem.daySecondToEpoch([ep.day(idx(k)) ep.sec(idx(k))], ep.scale);
                [v1, i1] = NonGravForces.env.wind("hwm14", e1, r(:, idx(k)), v(:, idx(k)), 'Eop', eop, 'SpaceWeather', sw);
                tc.verifyEqual(i4.hwm14.northEast_m_s(:, k), i1.hwm14.northEast_m_s); tc.verifyEqual(v4(:, k), v1);
            end
        end
        function densityDtm2020Vectorization(tc)
            % The gateway loops over epochs, so a batch equals the same epochs one at a time, bitwise.
            ep = referenceSystem.mjdToEpoch(58849 + (0:5)' / 6, "UTC");
            lla = [deg2rad([-80 -30 0 20 60 89]); deg2rad([0 60 120 180 240 300]); (480:4:500) * 1e3];
            d = struct('f107', 71, 'f107bar', 70, 'kp3', [0 1 2 3 4 5], 'kp24', 2);
            [rho, sp, T] = NonGravForces.env.density("dtm2020", ep, lla, [], 'Drivers', d);
            for k = 1:6
                dk = d; dk.kp3 = d.kp3(k);
                ek = referenceSystem.daySecondToEpoch([ep.day(k) ep.sec(k)], ep.scale);
                [r1, s1, t1] = NonGravForces.env.density("dtm2020", ek, lla(:, k), [], 'Drivers', dk);
                tc.verifyEqual([r1; s1; t1], [rho(k); sp(:, k); T(:, k)], sprintf('epoch %d batched vs alone', k));
            end
        end
        function densityVectorization(tc)
            % The 15 daily-ap published cases in one call equal the 15 single calls, bitwise (the MEX loop is per epoch).
            [rhoAll, spAll, TAll] = tc.msisCases(1:15, true);
            for i = 1:15
                [rho1, sp1, T1] = tc.msisCases(i, true);
                tc.verifyEqual(rhoAll(i), rho1); tc.verifyEqual(spAll(:, i), sp1); tc.verifyEqual(TAll(:, i), T1);
            end
            % msis21: the 200 NRL cases at once equal every tenth case alone (the gateway loops per epoch;
            % MSISCALC's cache of the previous input may only ever skip work, never change a value).
            [rho21, sp21, T21] = tc.msis21Cases(1:200);
            for i = 1:10:200
                [rho1, sp1, T1] = tc.msis21Cases(i);
                tc.verifyEqual(rho21(i), rho1); tc.verifyEqual(sp21(:, i), sp1); tc.verifyEqual(T21(:, i), T1);
            end
        end
        function magneticFieldVectorization(tc)
            r = [[6861.2e3; 0; 0] [1e6; -3e6; 6e6] [0; 0; -7e6] [3e6; 3e6; 5.5e6]];
            M = referenceSystem.rotationFromQuaternion([cos([0.1 0.5 0.9 1.3]); sin([0.1 0.5 0.9 1.3]); zeros(2, 4)]);
            B = NonGravForces.env.magneticField(tc.igrf, r, 'Rotation', M);
            for k = 1:4
                tc.verifyEqual(NonGravForces.env.magneticField(tc.igrf, r(:, k), 'Rotation', M(:, :, k)), B(:, k), sprintf('column %d batched vs alone', k));
            end
        end
        function eclipseFactorVectorization(tc)
            E = tc.eclipse; r = [E.rx E.ry E.rz]'; s = [E.sx E.sy E.sz]';
            nuAll = NonGravForces.env.eclipseFactor(r, s);
            one = arrayfun(@(k) NonGravForces.env.eclipseFactor(r(:, k), s(:, k)), 1:size(r, 2));
            tc.verifyEqual(nuAll, one);
            nuSAll = NonGravForces.env.eclipseFactor(r, s, 'Model', "solaarsCF", 'EclipticNormal', [0; 0; 1]);
            oneS = arrayfun(@(k) NonGravForces.env.eclipseFactor(r(:, k), s(:, k), 'Model', "solaarsCF", 'EclipticNormal', [0; 0; 1]), 1:size(r, 2));
            tc.verifyEqual(nuSAll, oneS);
        end
        function earthRadianceClimatologyVectorization(tc)
            clm = tc.syntheticClimatology(@(LAT, LON, m) 0.3 + 1e-4 * LON + 0.01 * m, ...
                                          @(LAT, LON, m) 200 + 0.1 * LAT, 90, 180, [58800 58830 58860]);
            ep = referenceSystem.mjdToEpoch([58810; 58815; 58840], "UTC");
            r = [[6.87e6; 1e5; 2e6] [1e5; 6.9e6; 0] [0; 1e3; 6.9e6]];
            t = NonGravForces.env.earthRadiance("ceresClimatology", ep, [6 12], r, 'Climatology', clm, 'RotationItrsToGcrs', tc.eyes(3));
            for k = 1:3
                % the epoch is INDEXED, not rebuilt: mjdToEpoch(day + sec/86400) is a different
                % quantity from the struct's own (day, sec) split and the shapes diverge with it
                epk = ep; epk.day = ep.day(k); epk.sec = ep.sec(k);
                one = NonGravForces.env.earthRadiance("ceresClimatology", epk, ...
                    [6 12], r(:, k), 'Climatology', clm, 'RotationItrsToGcrs', tc.eyes(1));
                tc.verifyEqual(one.albedo, t.albedo(:, k), sprintf('epoch %d batched vs alone', k));
                tc.verifyEqual(one.exitanceLw_W_m2, t.exitanceLw_W_m2(:, k));
                tc.verifyEqual(one.climatology.gridIndex, t.climatology.gridIndex(k));
            end
        end
        function earthRadianceCeresVectorization(tc)
            ep = referenceSystem.mjdToEpoch([58835.5; 58849.1; 58860.9], "UTC"); r = [[6.88e6; 0; 0] [0; 3e6; 6.2e6] [4e6; -4e6; 4e6]];
            R = referenceSystem.rotationFromQuaternion([cos([0.1 0.7 2.1]); zeros(2, 3); sin([0.1 0.7 2.1])]);
            t = NonGravForces.env.earthRadiance("ceres", ep, [12 24], r, 'Ceres', tc.cer, 'RotationItrsToGcrs', R);
            for k = 1:3
                epk = ep; epk.day = ep.day(k); epk.sec = ep.sec(k);
                t1 = NonGravForces.env.earthRadiance("ceres", epk, [12 24], r(:, k), 'Ceres', tc.cer, 'RotationItrsToGcrs', R(:, :, k));
                tc.verifyEqual(t.albedo(:, k), t1.albedo); tc.verifyEqual(t.exitanceLw_W_m2(:, k), t1.exitanceLw_W_m2); tc.verifyEqual(t.lat_rad(:, k), t1.lat_rad);
            end
        end
        function earthRadianceVectorization(tc)
            ep = referenceSystem.mjdToEpoch([58849.5; 58900.25; 58950.1], "UTC"); r = [[6.88e6; 0; 0] [0; 0; 7e6] [4e6; 4e6; 4e6]];
            t = NonGravForces.env.earthRadiance("knocke", ep, [12 24], r);
            for k = 1:3
                epk = ep; epk.day = ep.day(k); epk.sec = ep.sec(k);
                t1 = NonGravForces.env.earthRadiance("knocke", epk, [12 24], r(:, k));
                tc.verifyEqual(t.n(:, :, k), t1.n); tc.verifyEqual(t.area_m2(:, k), t1.area_m2); tc.verifyEqual(t.albedo(:, k), t1.albedo); tc.verifyEqual(t.emissivity(:, k), t1.emissivity);
            end
        end
        function radiationSourcesVectorization(tc)
            [ep, r, q] = tc.sampleState();
            idx = [1 100 250 400 600];
            epB = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], "GPS");
            src = NonGravForces.env.radiationSources(epB, r(:, idx), q(:, idx), NonGravForces.Options("standard"));
            for k = 1:5
                ep1 = referenceSystem.daySecondToEpoch([ep.day(idx(k)) ep.sec(idx(k))], "GPS");
                s1 = NonGravForces.env.radiationSources(ep1, r(:, idx(k)), q(:, idx(k)), NonGravForces.Options("standard"));
                tc.verifyEqual(src.P(:, k), s1.P); tc.verifyEqual(src.dir(:, :, k), s1.dir); tc.verifyEqual(src.nu(k), s1.nu);
            end
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function windRobustness(tc)
            [ep, r, ~, eop, v] = tc.sampleState();
            K = size(r, 2);
            w = repmat([10; -20; 5], 1, K);
            [v1, i1] = NonGravForces.env.wind("none", ep, r, v, 'Wind_m_s', w);
            tc.verifyEqual(v1, v - w); tc.verifyTrue(i1.windGiven);
            v2 = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop, 'Wind_m_s', w);
            v0 = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
            tc.verifyEqual(v2, v0 - w, 'the wind is removed after the co-rotation, not instead of it');
            % a satellite on the rotation axis has no co-rotation at all
            e1 = referenceSystem.daySecondToEpoch([ep.day(1) ep.sec(1)], ep.scale);
            R = referenceSystem.rotationItrsToGcrs(e1, 'Eop', eop);
            rPole = R * [0; 0; 6.87e6];
            [~, iP] = NonGravForces.env.wind("corotation", e1, rPole, v(:, 1), 'Eop', eop);
            tc.verifyBoundAndRecord("co-rotation on the rotation axis", norm(iP.corotation_m_s), 1e-2, "m/s", "omega x r vanishes on the axis; what remains is the polar-motion and nutation rate");
        end
        function windHwm14Robustness(tc)
            % Ground level, the model's top (the quiet wind is held constant above its highest level) and a storm ap all
            % evaluate; the current folder is the caller's afterwards (the gateway loads its tables from its own folder).
            eop = referenceSystem.readEop(fullfile(tc.root, 'data', 'eop', 'finals2000A_2019-09_2020-03.all'));
            ep = referenceSystem.mjdToEpoch([58849.5; 58849.5; 58849.5], "UTC");
            rI = referenceSystem.geodeticToCartesian([10 20 1e3; 10 20 600e3; -60 100 450e3], 'Ellipsoid', "WGS84", 'AngleUnit', "deg")';
            R = referenceSystem.rotationItrsToGcrs(ep, 'Eop', eop);
            rG = reshape(pagemtimes(R, reshape(rI, 3, 1, 3)), 3, 3);
            here = pwd;
            [~, info] = NonGravForces.env.wind("hwm14", ep, rG, zeros(3, 3), 'Eop', eop, 'Ap', [4 4 400]);
            tc.verifyEqual(pwd, here, 'the current folder must survive the call');
            tc.verifyTrue(all(isfinite(info.hwm14.northEast_m_s(:))) && all(abs(info.hwm14.northEast_m_s(:)) < 1000));
            e1 = referenceSystem.mjdToEpoch(58849.5, "UTC");
            [~, i600] = NonGravForces.env.wind("hwm14", e1, rG(:, 2) * 1.02, zeros(3, 1), 'Eop', eop, 'Ap', 4);
            tc.verifyTrue(all(isfinite(i600.hwm14.northEast_m_s)), 'above the model''s top it still answers');
            % HWM14 does not return for a negative altitude (measured in a standalone build): refused before the model
            below = R(:, :, 1) * referenceSystem.geodeticToCartesian([10 20 -1], 'Ellipsoid', "WGS84", 'AngleUnit', "deg")';
            tc.verifyErrorId(@() NonGravForces.env.wind("hwm14", e1, below, zeros(3, 1), 'Eop', eop, 'Ap', 4), 'NonGravForces:wind:inputRange');
        end
        function densityDtm2020Robustness(tc)
            % The edges of the drivers and of the position all evaluate to finite positive densities, and a
            % local time of 24 h is 0 h.
            e = referenceSystem.mjdToEpoch(58849.5 * ones(4, 1), "UTC");
            lla = [deg2rad([90 -90 0 0]); 0 0 0 0; 490e3 490e3 150e3 1500e3];
            for kp = [0 9]
                rho = NonGravForces.env.density("dtm2020", e, lla, [], 'Drivers', struct('f107', 65, 'f107bar', 65, 'kp3', kp, 'kp24', kp));
                tc.verifyTrue(all(isfinite(rho) & rho > 0), sprintf('Kp %d at the poles, at 150 km and at 1500 km', kp));
            end
            d = struct('f107', 70, 'f107bar', 70, 'kp3', 2, 'kp24', 2);
            e1 = referenceSystem.mjdToEpoch(58849.5, "UTC");
            r0 = NonGravForces.env.density("dtm2020", e1, [0; 0; 490e3], [], 'Drivers', d, 'LocalSolarTime_h', 0);
            r24 = NonGravForces.env.density("dtm2020", e1, [0; 0; 490e3], [], 'Drivers', d, 'LocalSolarTime_h', 24);
            tc.verifyEqual(r24, r0, 'local time 24 h is 0 h');
            tc.verifyErrorId(@() NonGravForces.env.density("dtm2020", e1, [0; 0; 490e3], [], 'Drivers', d, 'LocalSolarTime_h', [1 2]), 'NonGravForces:density:inputDrivers');
        end
        function densityRobustness(tc)
            ep = referenceSystem.mjdToEpoch(58849.5, "UTC"); d = struct('f107', 70, 'f107a', 71, 'ap', 3);
            tc.verifyErrorId(@() NonGravForces.env.density("msis00", ep, [0; NaN; 490e3], [], 'Drivers', d), 'NonGravForces:density:numericNaN');
            tc.verifyErrorId(@() NonGravForces.env.density("msis00", ep, [0; 0; 490e3], [], 'Drivers', struct('f107', NaN, 'f107a', 71, 'ap', 3)), 'NonGravForces:density:numericNaN');
            tc.verifyErrorId(@() NonGravForces.env.density("msis00", ep, [0; 0; 490e3], [], 'Drivers', struct('f107', [70 70], 'f107a', 71, 'ap', 3)), 'NonGravForces:density:inputDrivers');
            tc.verifyErrorId(@() NonGravForces.env.density("msis00", ep, [0; 0; 490e3], [], 'Drivers', d, 'LocalSolarTime_h', [1 2]), 'NonGravForces:density:inputDrivers');
            sj = NonGravForces.io.spaceWeather([string(fullfile(tc.root, 'data', 'spaceweather', 'SOLFSMY_2019-09_2020-03.TXT')) string(fullfile(tc.root, 'data', 'spaceweather', 'DTCFILE_2019-09_2020-03.TXT'))]);
            tc.verifyErrorId(@() NonGravForces.env.density("msis00", ep, [0; 0; 490e3], sj), 'NonGravForces:density:inputDrivers');
            % real drivers from the CelesTrak fixture, ground level and 1000 km both evaluate
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            [rho, ~, ~, info] = NonGravForces.env.density("msis00", referenceSystem.mjdToEpoch([58849.5; 58849.5], "UTC"), [0 1; 0 1; 0 1000e3], sw);
            tc.verifyTrue(rho(1) > 1 && rho(1) < 1.5, 'sea-level density about 1.2 kg/m^3'); tc.verifyTrue(rho(2) > 0 && rho(2) < 1e-13);
            tc.verifyEqual(info.drivers.source, "spaceWeather"); tc.verifySize(info.drivers.ap, [7 2]);
            % msis21 on the same real drivers: ground level has no atomic O, H or N (returned 0, not 9.999e-38);
            % 1000 km evaluates. The gateway re-initialises when the ap mode or the mass flag changes, and must
            % come back to the same state: a daily-Ap call after a history call and a gtd7-like call is bitwise the first.
            [rho21, sp21] = NonGravForces.env.density("msis21", referenceSystem.mjdToEpoch([58849.5; 58849.5], "UTC"), [0 1; 0 1; 0 1000e3], sw);
            tc.verifyTrue(rho21(1) > 1 && rho21(1) < 1.5, 'sea-level density about 1.2 kg/m^3'); tc.verifyTrue(rho21(2) > 0 && rho21(2) < 1e-13);
            tc.verifyEqual(sp21([2 6 7], 1), [0; 0; 0], 'no O, H, N at ground level: 0, not the model''s missing marker');
            tc.verifyErrorId(@() NonGravForces.env.density("msis21", ep, [0; NaN; 490e3], [], 'Drivers', d), 'NonGravForces:density:numericNaN');
            first = NonGravForces.env.density("msis21", ep, [0.3; 1.1; 490e3], [], 'Drivers', d);
            NonGravForces.env.density("msis21", ep, [0.3; 1.1; 490e3], [], 'Drivers', struct('f107', 70, 'f107a', 71, 'ap', 3 * ones(7, 1)));
            NonGravForces.env.density("msis21", ep, [0.3; 1.1; 490e3], [], 'Drivers', d, 'AnomalousOxygen', false);
            tc.verifyEqual(NonGravForces.env.density("msis21", ep, [0.3; 1.1; 490e3], [], 'Drivers', d), first, 're-initialisation returns the gateway to the same state');
            % the two climatologies on the real day at GRACE-FO height: tens of percent apart, never a factor
            e490 = referenceSystem.mjdToEpoch((58849:0.125:58849.875)', "UTC"); l490 = [linspace(-1.4, 1.4, 8); linspace(-3, 3, 8); 490e3 * ones(1, 8)];
            q = NonGravForces.env.density("msis21", e490, l490, sw) ./ NonGravForces.env.density("msis00", e490, l490, sw);
            tc.verifyBoundAndRecord("msis21 / msis00 at 490 km on 2020-01-01, largest departure of the ratio from 1", max(abs(log(q))), log(2), "-", "two empirical models of the same thermosphere: a plausibility bound on the wiring (species order, units), not a model comparison");
        end
        function magneticFieldRobustness(tc)
            % The geographic poles are ordinary points: the field there is the limit from a microradian away.
            rp = 6861.2e3;
            for sgn = [1 -1]
                Bpole = NonGravForces.env.magneticField(tc.igrf, [0; 0; sgn * rp]);
                Bnear = NonGravForces.env.magneticField(tc.igrf, rp * [1e-6 * cos(0.7); 1e-6 * sin(0.7); sgn * sqrt(1 - 1e-12)]);
                tc.verifyTrue(all(isfinite(Bpole)));
                tc.verifyBoundAndRecord("magneticField at a pole vs 1 microradian off it (" + string(sgn) + ")", norm(Bpole - Bnear) / norm(Bpole), 1e-5, "-", ...
                    "continuity: the horizontal gradient over 7 m of arc is below 1e-6 of the field");
            end
            % Far away the field falls as r^-3 and stays finite.
            Bfar = NonGravForces.env.magneticField(tc.igrf, [1e9; 0; 0]);
            tc.verifyTrue(all(isfinite(Bfar)) && norm(Bfar) < 1e-11);
        end
        function eclipseFactorRobustness(tc)
            s = [tc.au; 0; 0]; r = [-(6378137 + 490e3) * cosd(68.25); 0; (6378137 + 490e3) * sind(68.25)];
            nu = NonGravForces.env.eclipseFactor(r, s); tc.verifyTrue(nu > 0 && nu < 1, 'penumbra point');
            tc.verifyEqual(NonGravForces.env.eclipseFactor(r, s, 'Moon', [2 * tc.au; 0; 0]), nu, 'a body beyond the Sun casts no shadow');
            tc.verifyEqual(NonGravForces.env.eclipseFactor(r, s, 'Moon', [0; 0; -384400e3]), nu, 'a Moon off the line leaves nu');
            tc.verifyEqual(NonGravForces.env.eclipseFactor([7e6; 0; 0], s, 'Moon', [3e8; 0; 0]), 0, 'a Moon between satellite and Sun gives the lunar umbra');
            tc.verifyErrorId(@() NonGravForces.env.eclipseFactor([NaN; 0; 7e6], s), 'NonGravForces:eclipseFactor:numericNaN');
            tc.verifyErrorId(@() NonGravForces.env.eclipseFactor(r, [NaN; 0; 0]), 'NonGravForces:eclipseFactor:numericNaN');
            tc.verifyErrorId(@() NonGravForces.env.eclipseFactor(r, s, 'Moon', [1; 2]), 'NonGravForces:eclipseFactor:inputShape');
            tc.verifyEqual(NonGravForces.env.eclipseFactor(zeros(3, 0), s), zeros(1, 0), 'empty in, empty out');
        end
        function earthRadianceClimatologyRobustness(tc)
            clm = tc.syntheticClimatology(@(LAT, LON, m) 0.3 + 0 * LAT, @(LAT, LON, m) 220 + 0 * LAT, 20, 40, [58800 58830 58860]);
            ep = referenceSystem.mjdToEpoch(58810, "UTC");
            r = [6.87e6; 0; 0];
            f = @(varargin) NonGravForces.env.earthRadiance(varargin{:});
            % the contradictions: a grid to a backend that does not read it, and a backend without its grid
            tc.verifyErrorId(@() f("knocke", ep, [4 8], r, 'Climatology', clm), 'NonGravForces:earthRadiance:inputClimatology');
            tc.verifyErrorId(@() f("ceresClimatology", ep, [4 8], r, 'RotationItrsToGcrs', tc.eyes(1)), 'NonGravForces:earthRadiance:inputClimatology');
            tc.verifyErrorId(@() f("ceresClimatology", ep, [4 8], r, 'Climatology', clm), 'NonGravForces:earthRadiance:inputClimatology');
            % the two files are separate products and only this check pairs them
            bad = clm; bad.longwave.lat_deg = bad.longwave.lat_deg * 0.9;
            tc.verifyErrorId(@() f("ceresClimatology", ep, [4 8], r, 'Climatology', bad, 'RotationItrsToGcrs', tc.eyes(1)), 'NonGravForces:earthRadiance:inputClimatology');
            bad2 = clm; bad2.longwave.values = bad2.longwave.values(:, :, 1);
            tc.verifyErrorId(@() f("ceresClimatology", ep, [4 8], r, 'Climatology', bad2, 'RotationItrsToGcrs', tc.eyes(1)), 'NonGravForces:earthRadiance:inputClimatology');
            % a node count that does not bound the grids, and a non-uniform longitude axis
            bad3 = clm; bad3.albedo.mjd = bad3.albedo.mjd(1:2); bad3.longwave.mjd = bad3.longwave.mjd(1:2);
            tc.verifyErrorId(@() f("ceresClimatology", ep, [4 8], r, 'Climatology', bad3, 'RotationItrsToGcrs', tc.eyes(1)), 'NonGravForces:earthRadiance:inputClimatology');
            bad4 = clm; bad4.albedo.lon_deg(3) = bad4.albedo.lon_deg(3) + 1; bad4.longwave.lon_deg = bad4.albedo.lon_deg;
            tc.verifyErrorId(@() f("ceresClimatology", ep, [4 8], r, 'Climatology', bad4, 'RotationItrsToGcrs', tc.eyes(1)), 'NonGravForces:earthRadiance:inputClimatology');
            % radii that are not one height above the ellipsoid: the check that caught the convention
            bad5 = clm; bad5.albedo.radius_m(4) = bad5.albedo.radius_m(4) + 500;
            tc.verifyErrorId(@() f("ceresClimatology", ep, [4 8], r, 'Climatology', bad5, 'RotationItrsToGcrs', tc.eyes(1)), 'NonGravForces:earthRadiance:inputClimatology');
            % an epoch outside the nodes must refuse, not extrapolate silently
            far = referenceSystem.mjdToEpoch(59000, "UTC");
            tc.verifyErrorId(@() f("ceresClimatology", far, [4 8], r, 'Climatology', clm, 'RotationItrsToGcrs', tc.eyes(1)), 'NonGravForces:earthRadiance:dataRange');
            % and the closing node itself is inside, not outside
            edge = NonGravForces.env.earthRadiance("ceresClimatology", referenceSystem.mjdToEpoch(58860, "UTC"), [4 8], r, 'Climatology', clm, 'RotationItrsToGcrs', tc.eyes(1));
            tc.verifyEqual(edge.climatology.gridIndex, 2, 'the last node closes the last grid rather than falling past it');
        end
        function earthRadianceCeresRobustness(tc)
            ep = referenceSystem.mjdToEpoch([58849.5; 58849.5], "UTC");
            a = referenceSystem.ellipsoidParameters("WGS84").a;
            % over the pole, and over the Arctic in polar night where the product's SWdown is 0: finite, albedo 0 there
            r = [[0; 0; a + 490e3] [(a + 490e3) * cosd(80); 0; (a + 490e3) * sind(80)]];
            t = NonGravForces.env.earthRadiance("ceres", ep, [12 24], r, 'Ceres', tc.cer, 'RotationItrsToGcrs', tc.eyes(2));
            tc.verifyTrue(all(isfinite([t.albedo(:); t.exitanceLw_W_m2(:)])));
            dark = rad2deg(t.lat_rad) > 88;
            tc.verifyEqual(t.albedo(dark), zeros(nnz(dark), 1), 'no insolation in the month means albedo 0, not 0/0');
            % below the 20 km reference level the cap is refused, although it is above the ellipsoid
            tc.verifyErrorId(@() NonGravForces.env.earthRadiance("ceres", referenceSystem.mjdToEpoch(58849.5, "UTC"), [12 24], [a + 1e4; 0; 0], 'Ceres', tc.cer, 'RotationItrsToGcrs', eye(3)), 'NonGravForces:earthRadiance:inputRange');
            % the last node is inside the held months, one second past it is not
            eLast = referenceSystem.mjdToEpoch(tc.cer.mjd(end), "UTC");
            tc.verifyWarningFree(@() NonGravForces.env.earthRadiance("ceres", eLast, [2 4], [7e6; 0; 0], 'Ceres', tc.cer, 'RotationItrsToGcrs', eye(3)));
            ePast = referenceSystem.mjdToEpoch(tc.cer.mjd(end) + 1 / 86400, "UTC");
            tc.verifyErrorId(@() NonGravForces.env.earthRadiance("ceres", ePast, [2 4], [7e6; 0; 0], 'Ceres', tc.cer, 'RotationItrsToGcrs', eye(3)), 'NonGravForces:earthRadiance:dataRange');
        end
        function earthRadianceRobustness(tc)
            % Over the pole the tangent basis needs its fallback; a single tile; just above the surface.
            ep = referenceSystem.mjdToEpoch(58849.5, "UTC");
            tP = NonGravForces.env.earthRadiance("knocke", ep, [12 24], [0; 0; 7e6]);
            tc.verifyLessThan(abs(sum(tP.area_m2) / (2 * pi * tP.radius_m^2 * (1 - cos(tP.capAngle_rad))) - 1), 1e-14);
            tc.verifyTrue(all(isfinite(tP.n(:))));
            t1 = NonGravForces.env.earthRadiance("knocke", ep, [1 1], [7e6; 0; 0]); tc.verifyEqual(t1.M, 1);
            tc.verifyLessThan(abs(t1.area_m2 - 2 * pi * t1.radius_m^2 * (1 - cos(t1.capAngle_rad))), 1e-3);
            tLow = NonGravForces.env.earthRadiance("knocke", ep, [12 24], [6378137 + 1; 0; 0]); tc.verifyLessThan(tLow.capAngle_rad, 1e-3);
            tc.verifyErrorId(@() NonGravForces.env.earthRadiance("knocke", ep, [12 24], [NaN; 0; 7e6]), 'NonGravForces:earthRadiance:numericNaN');
            tc.verifyErrorId(@() NonGravForces.env.earthRadiance("knocke", struct('x', 1), [12 24], [7e6; 0; 0]), 'NonGravForces:earthRadiance:inputShape');
        end
        function radiationSourcesRobustness(tc)
            ep = referenceSystem.mjdToEpoch(58849.5, "UTC"); r = [-6.88e6; 0; 0]; q = [1; 0; 0; 0]; sun = [tc.au; 0; 0];
            sNone = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard", "erp.backend", "none"), 'Sun', sun);
            tc.verifyEqual(sNone.M, 1); tc.verifyEqual(sNone.P, 0, 'umbra: no solar pressure'); tc.verifyEqual(sNone.band, "vis");
            sPt = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard", "erp.backend", "pointSource"), 'Sun', sun);
            tc.verifyEqual(sPt.M, 3); tc.verifyEqual(sPt.band, ["vis","vis","ir"]); tc.verifyLessThan(sPt.P(2), 1e-20, 'anti-solar: albedo is the rounding of sin(pi) in the phase law'); tc.verifyGreaterThan(sPt.P(3), 0);
            sL = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard", "srp.lunarEclipse", true));
            tc.verifySize(sL.nu, [1 1]);
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard", "srp.lunarEclipse", true), 'Sun', sun), 'NonGravForces:radiationSources:inputShape');
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep, r, [NaN; 0; 0; 0], NonGravForces.Options("standard")), 'NonGravForces:radiationSources:numericNaN');
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("precise")), 'NonGravForces:radiationSources:depKernel');
            % solaarsCF through the ephemeris path (ecliptic normal measured from the Sun state) and the caller path
            sSol = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard", "srp.eclipse", "solaarsCF"));
            tc.verifyEqual(sSol.meta.eclipse.model, "solaarsCF"); tc.verifyTrue(sSol.nu >= 0 && sSol.nu <= 1);
            tc.verifyErrorId(@() NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard", "srp.eclipse", "solaarsCF"), 'Sun', sun), 'NonGravForces:radiationSources:inputShape');
            sSol2 = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard", "srp.eclipse", "solaarsCF"), 'Sun', sun, 'EclipticNormal', [0; 0; 1]);
            tc.verifyEqual(sSol2.nu, 0, 'the anti-solar point is deep umbra in SOLAARS-CF too');
            % a non-unit quaternion is normalised by referenceSystem: directions agree to 1e-15
            s1 = NonGravForces.env.radiationSources(ep, r, q, NonGravForces.Options("standard"), 'Sun', sun);
            s2 = NonGravForces.env.radiationSources(ep, r, 2 * q, NonGravForces.Options("standard"), 'Sun', sun);
            tc.verifyLessThan(max(abs(s1.dir(:) - s2.dir(:))), 1e-15); tc.verifyEqual(s1.P, s2.P);
        end
        function buildMexRobustness(tc)
            % A tampered digest makes the state STALE by content; restoring it makes it current again.
            rep = NonGravForces.buildMex('Check', true, 'Target', "msis");
            digestFile = fullfile(fileparts(rep.file), 'msisGateway.digest.txt');
            original = fileread(digestFile);
            tc.addTeardown(@() writelines(string(strtrim(original)), digestFile));
            writelines("0000", digestFile);
            tampered = NonGravForces.buildMex('Check', true, 'Target', "msis");
            tc.verifyTrue(tampered.present); tc.verifyFalse(tampered.upToDate); tc.verifyEqual(tampered.binaryDigest, "0000");
            writelines(string(strtrim(original)), digestFile);
            tc.verifyTrue(NonGravForces.buildMex('Check', true, 'Target', "msis").upToDate);
        end
    end

    % ===================================================================== metamorphic
    methods (Test, TestTags = {'metamorphic'})
        function windMetamorphic(tc)
            [ep, r, ~, eop, v] = tc.sampleState();
            [~, i1] = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
            [~, i2] = NonGravForces.env.wind("corotation", ep, 2 * r, v, 'Eop', eop);
            tc.verifyBoundAndRecord("co-rotation doubles with the position", max(abs(i2.corotation_m_s - 2 * i1.corotation_m_s), [], 'all') / max(abs(i1.corotation_m_s), [], 'all'), 1e-15, "-", "an exact power of two through two matrix products");
            [~, i3] = NonGravForces.env.wind("corotation", ep, r, 3 * v, 'Eop', eop);
            tc.verifyEqual(i3.corotation_m_s, i1.corotation_m_s, 'the co-rotation cannot depend on the satellite velocity');
            % Adding the co-rotation back as a wind field returns the inertial velocity, but (v - c) + c
            % is not v in floating point: the claim is one rounding of a 7.6 km/s quantity, not equality.
            vn = NonGravForces.env.wind("none", ep, r, v);
            vc = NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop, 'Wind_m_s', -i1.corotation_m_s);
            tc.verifyBoundAndRecord("co-rotation removed and added back", max(abs(vc - vn), [], 'all') / max(abs(vn), [], 'all'), 1e-15, "-", "two roundings of a 7.6 km/s velocity; measured 1.3e-16");
        end
        function windHwm14Metamorphic(tc)
            % The model wind is additive: hwm14 minus corotation is exactly the model wind; quiet plus DWM07 is the total
            % (to the model's single precision); and at ap = 0 the disturbance is not zero but DWM07's quiet level.
            [ep, r, ~, eop, v] = tc.sampleState();
            idx = 1:50:600;
            e = referenceSystem.daySecondToEpoch([ep.day(idx) ep.sec(idx)], ep.scale);
            vc = NonGravForces.env.wind("corotation", e, r(:, idx), v(:, idx), 'Eop', eop);
            [vh, ih] = NonGravForces.env.wind("hwm14", e, r(:, idx), v(:, idx), 'Eop', eop, 'Disturbed', false);
            tc.verifyBoundAndRecord("hwm14 minus corotation is the model wind", max(abs(vc - vh - ih.hwm14.gcrs_m_s), [], 'all'), 1e-9, "m/s", "(v - c) - (v - c - w) = w, three roundings of a 7.6 km/s velocity");
            [~, iq] = NonGravForces.env.wind("hwm14", e, r(:, idx), v(:, idx), 'Eop', eop, 'Ap', 30);
            [~, iq2] = NonGravForces.env.wind("hwm14", e, r(:, idx), v(:, idx), 'Eop', eop, 'Ap', 30);
            tc.verifyEqual(iq2.hwm14.northEast_m_s, iq.hwm14.northEast_m_s, 'repeated calls at the same inputs are bitwise (no internal cursor)');
            tc.verifyNotEqual(iq.hwm14.northEast_m_s, ih.hwm14.northEast_m_s, 'ap 30 must add a disturbance wind at 490 km');
        end
        function densityDtm2020Metamorphic(tc)
            % More solar flux, more density; higher, less. And beside NRLMSISE-00 on the reference day at GRACE-FO's
            % altitude -- the number O-19 asks about: CNES fitted DTM2020 to the TU Delft CHAMP, GOCE and Swarm
            % densities derived with high-fidelity satellite models.
            e1 = referenceSystem.mjdToEpoch(58849.5, "UTC");
            f = [70 120 180];
            rho = arrayfun(@(x) NonGravForces.env.density("dtm2020", e1, [0; 0; 490e3], [], 'Drivers', struct('f107', x, 'f107bar', x, 'kp3', 2, 'kp24', 2)), f);
            tc.verifyTrue(all(diff(rho) > 0), 'density rises with F10.7');
            rh = arrayfun(@(h) NonGravForces.env.density("dtm2020", e1, [0; 0; h], [], 'Drivers', struct('f107', 70, 'f107bar', 70, 'kp3', 2, 'kp24', 2)), [300 400 500] * 1e3);
            tc.verifyTrue(all(diff(rh) < 0), 'density falls with altitude');
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            [lat, lon] = meshgrid(deg2rad(-80:20:80), deg2rad(0:45:315));
            K = numel(lat);
            ep = referenceSystem.mjdToEpoch(58849 + (0:K-1)' / K, "UTC");
            lla = [lat(:)'; lon(:)'; 490e3 * ones(1, K)];
            rD = NonGravForces.env.density("dtm2020", ep, lla, sw);
            rM = NonGravForces.env.density("msis00", ep, lla, sw);
            rN = NonGravForces.env.density("msis21", ep, lla, sw);
            tc.verifyBoundAndRecord("dtm2020 / msis00 at 490 km, 2020-01-01: |log| of the median ratio", abs(log(median(rD ./ rM))), log(3), "-", ...
                "first run 2026-10-27: recorded; the solar-minimum density level O-19 asks about");
            tc.verifyBoundAndRecord("dtm2020 / msis21 at 490 km, 2020-01-01: |log| of the median ratio", abs(log(median(rD ./ rN))), log(3), "-", ...
                "first run 2026-10-27: recorded");
            fprintf('densityDtm2020Metamorphic: 490 km, 2020-01-01, %d points: median dtm2020/msis00 %.4f, dtm2020/msis21 %.4f, msis21/msis00 %.4f\n', ...
                K, median(rD ./ rM), median(rD ./ rN), median(rN ./ rM));
        end
        function densityMetamorphic(tc)
            % gtd7d - gtd7: the mass density gains exactly 16 * 1.66e-24 g times the anomalous-oxygen number density
            % (the reference's own arithmetic, nrlmsise-00.c gtd7d), in SI 16 * 1.66e-24 * 1e-3 * 1e6 kg per m^-3.
            [rhoD, spD] = tc.msisCases(1:15, true); [rho7, sp7] = tc.msisCases(1:15, false);
            tc.verifyEqual(spD, sp7, 'the number densities do not depend on the flag');
            % SI: 16 * 1.66e-24 g per particle = 16 * 1.66e-24 * 1e-3 kg, times the number density in m^-3
            tc.verifyBoundAndRecord("gtd7d - gtd7 vs 16 amu x anomalous O", max(abs((rhoD - rho7) - 16 * 1.66e-24 * 1e-3 * spD(8, :)) ./ rhoD), 1e-12, "rel", "two rounding operations, per case");
            % msis21: the model's own mass flag drops exactly the anomalous-oxygen term of its dot product,
            % whose weight is 31.9988/2 g/mol over 1e3 N_A (msis_constants.F90); number densities untouched.
            idx = find(tc.msis21In(:, 3) > 100)';
            [r21, s21] = tc.msis21Cases(idx); [r21n, s21n] = tc.msis21Cases(idx, 'AnomalousOxygen', false);
            tc.verifyEqual(s21n, s21, 'msis21 number densities do not depend on the mass flag');
            tc.verifyBoundAndRecord("msis21 rho - rho(no anomalous O) vs n(O*) m_O", max(abs((r21 - r21n) - s21(8, :) * (31.9988 / 2) / (1e3 * 6.02214076e23)) ./ r21), 1e-14, "rel", "a nine-term and an eight-term dot product; a few eps of the total");
            % exponential: one scale height down multiplies by e, whatever h0
            ep = referenceSystem.mjdToEpoch([58849.5; 58849.5], "UTC"); p = struct('rho0_kg_m3', 3e-12, 'h0_m', 450e3, 'H_m', 60e3);
            rho = NonGravForces.env.density("exponential", ep, [0 0; 0 0; 400e3 340e3], [], 'Exponential', p);
            tc.verifyBoundAndRecord("exponential scale-height ratio", abs(rho(2) / rho(1) - exp(1)) / exp(1), 1e-14, "rel", "");
        end
        function magneticFieldMetamorphic(tc)
            % A dipole is even in position and falls as r^-3; the rotation option is the rotation of the Earth-fixed field.
            r = [[6.8e6; 1e6; 2e6] [-3e6; 2e6; 5.5e6]];
            B1 = NonGravForces.env.magneticField(tc.igrf, r, 'MaxDegree', 1);
            tc.verifyBoundAndRecord("dipole field is even in position", max(vecnorm(NonGravForces.env.magneticField(tc.igrf, -r, 'MaxDegree', 1) - B1) ./ vecnorm(B1)), 1e-15, "-", "cos and sin of phi + pi and of pi - theta");
            tc.verifyBoundAndRecord("dipole field falls as r^-3", max(vecnorm(NonGravForces.env.magneticField(tc.igrf, 2 * r, 'MaxDegree', 1) - B1 / 8) ./ vecnorm(B1)), 1e-15, "-", "an exact power of two");
            M = referenceSystem.rotationFromQuaternion([cos([0.4 1.1]); zeros(1, 2); sin([0.4 1.1]); zeros(1, 2)]);
            [Bm, info] = NonGravForces.env.magneticField(tc.igrf, r, 'Rotation', M);
            tc.verifyEqual(Bm, reshape(pagemtimes(M, reshape(info.itrs_T, 3, 1, 2)), 3, 2), 'the rotation is applied to the Earth-fixed field and to nothing else');
            tc.verifyBoundAndRecord("rotation preserves the intensity", max(abs(vecnorm(Bm) - info.intensity_T) ./ info.intensity_T), 1e-15, "-", "an orthogonal matrix");
        end
        function eclipseFactorMetamorphic(tc)
            % Mirror symmetry across the Sun line (bitwise), rotation of the whole geometry (1e-15),
            % monotone through the penumbra along the oracle sweep.
            E = tc.eclipse; r = [E.rx E.ry E.rz]'; s = [E.sx E.sy E.sz]';
            nu = NonGravForces.env.eclipseFactor(r, s);
            tc.verifyEqual(NonGravForces.env.eclipseFactor(r .* [1; -1; -1], s), nu, 'mirror y,z');
            Rm = referenceSystem.rotationFromQuaternion([cos(0.4); sin(0.4) * [0.6; 0; 0.8]]);
            % HANDOVER accuracy: nu to 1e-12. Position rounding of a rotation (1e-9 m) through the penumbra
            % slope (~100 per rad) is 1e-14; the atan2 separation angle keeps it there (acos lost 8e-12).
            tc.verifyBoundAndRecord("eclipse rotation invariance", max(abs(NonGravForces.env.eclipseFactor(Rm * r, Rm * s) - nu)), 1e-12, "-", "HANDOVER accuracy figure; measured see last_run");
            sweep = nu(1:30);                                        % oracle cases 0..29 sweep psi from 0 to 180 deg
            tc.verifyTrue(all(diff(sweep) >= 0), 'nu is non-decreasing from the anti-Sun point outwards');
            % SOLAARS-CF: monotone and smooth through the penumbra; oblateness acts on the ecliptic-parallel component
            % only, so a displacement rE along the normal equals a displacement rE / sO in the ecliptic plane.
            sAu = [tc.au; 0; 0];
            rE = linspace(6.30e6, 6.50e6, 2001); rSw = [-2.56e6 * ones(1, 2001); zeros(1, 2001); rE];
            nuS = NonGravForces.env.eclipseFactor(rSw, sAu, 'Model', "solaarsCF", 'EclipticNormal', [0; 0; 1]);
            tc.verifyTrue(all(diff(nuS) >= 0), 'SOLAARS-CF monotone in rE');
            % tails of the tanh sum 100 km inside / outside the geometric limb (fit tails, not physics): measured 2026-09-07
            tc.verifyBoundAndRecord("SOLAARS-CF tail 100 km inside the limb", nuS(1), 1e-4, "-", "");
            tc.verifyBoundAndRecord("SOLAARS-CF tail 100 km outside the limb", 1 - nuS(end), 1e-4, "-", "");
            d1 = diff(nuS) / 100;
            tc.verifyBoundAndRecord("SOLAARS-CF derivative jump (C1 proxy, 100 m steps)", max(abs(diff(d1))) / max(abs(d1)), 1e-2, "-", "smooth tanh sum: second difference small against the slope");
            ell = referenceSystem.ellipsoidParameters("WGS84"); sO = ell.a / ell.b;
            rPar = [-2.56e6 * ones(1, 2001); rE / sO; zeros(1, 2001)];
            nuP = NonGravForces.env.eclipseFactor(rPar, sAu, 'Model', "solaarsCF", 'EclipticNormal', [0; 0; 1]);
            tc.verifyBoundAndRecord("SOLAARS-CF oblateness scaling relation", max(abs(nuP - nuS)), 1e-9, "-", "rE/sO*sO rounding through steep sigmoids (slope ~1e-4 per m)");
            nuCn = NonGravForces.env.eclipseFactor(rSw, sAu);
            [~, iS] = min(abs(nuS - 0.5)); [~, iC] = min(abs(nuCn - 0.5));
            tc.verifyBoundAndRecord("SOLAARS-CF vs conical half-shadow radius shift", abs(rE(iS) - rE(iC)), 5e4, "m", "refraction and oblateness; measured see last_run");
        end
        function earthRadianceClimatologyMetamorphic(tc)
            % PIECEWISE CONSTANT IN TIME is the property that separates this backend from ceres, and
            % it is asserted as an identity rather than described in a comment: two epochs inside one
            % month must give bitwise equal tiles, and the step must land exactly at the node.
            clm = tc.syntheticClimatology(@(LAT, LON, m) 0.2 + 0.1 * m, @(LAT, LON, m) 200 + 10 * m, 40, 80, [58800 58830 58860]);
            r = [6.87e6; 0; 0];
            R1 = tc.eyes(1);
            a = NonGravForces.env.earthRadiance("ceresClimatology", referenceSystem.mjdToEpoch(58801, "UTC"), [4 8], r, 'Climatology', clm, 'RotationItrsToGcrs', R1);
            b = NonGravForces.env.earthRadiance("ceresClimatology", referenceSystem.mjdToEpoch(58829, "UTC"), [4 8], r, 'Climatology', clm, 'RotationItrsToGcrs', R1);
            tc.verifyEqual(b.albedo, a.albedo, 'two epochs in one month are the same month, bitwise');
            tc.verifyEqual(b.exitanceLw_W_m2, a.exitanceLw_W_m2);
            c = NonGravForces.env.earthRadiance("ceresClimatology", referenceSystem.mjdToEpoch(58830, "UTC"), [4 8], r, 'Climatology', clm, 'RotationItrsToGcrs', R1);
            tc.verifyEqual(c.climatology.gridIndex, 2, 'the node opens the next month');
            tc.verifyBoundAndRecord("ceresClimatology steps by exactly one month's difference at the node", max(abs((c.albedo - a.albedo) - 0.1), [], 'all'), 4 * eps, "-", "0.1 per month by construction; a linear blend would give a fraction of it");
            % a field constant over the sphere comes back constant: the interpolation adds nothing
            flat = tc.syntheticClimatology(@(LAT, LON, m) 0.37 + 0 * LAT, @(LAT, LON, m) 231 + 0 * LAT, 30, 60, [58800 58830]);
            g = NonGravForces.env.earthRadiance("ceresClimatology", referenceSystem.mjdToEpoch(58810, "UTC"), [6 12], r, 'Climatology', flat, 'RotationItrsToGcrs', R1);
            tc.verifyBoundAndRecord("ceresClimatology of a constant field is that constant", max(abs(g.albedo - 0.37), [], 'all'), 4 * eps, "-", "bilinear weights sum to one");
            % longitude is periodic: the same point addressed from either side of the seam agrees
            seam = tc.syntheticClimatology(@(LAT, LON, m) 0.3 + 1e-3 * cosd(LON), @(LAT, LON, m) 220 + 0 * LAT, 30, 60, [58800 58830]);
            e1 = NonGravForces.env.earthRadiance("ceresClimatology", referenceSystem.mjdToEpoch(58810, "UTC"), [4 8], [-6.87e6; 0; 0], 'Climatology', seam, 'RotationItrsToGcrs', R1);
            tc.verifyTrue(all(isfinite(e1.albedo), 'all'), 'the anti-meridian is interpolated, not dropped');
        end
        function earthRadianceCeresMetamorphic(tc)
            % Turning the Earth by exactly one cell and shifting the fields by one column back leaves every tile value:
            % the ITRS longitudes all move by 1 deg and read the same numbers (rounding of atan2 and of the weights).
            ep = referenceSystem.mjdToEpoch([58845.3; 58858.7], "UTC"); r = [[6.87e6; 1e6; 5e5] [-2e6; 3e6; 5.8e6]];
            R = referenceSystem.rotationFromQuaternion([cos([0.3 1.2]); zeros(2, 2); sin([0.3 1.2])]);
            t1 = NonGravForces.env.earthRadiance("ceres", ep, [12 24], r, 'Ceres', tc.cer, 'RotationItrsToGcrs', R);
            c2 = tc.cer;
            for name = ["rsut" "rlut" "rsdt"], c2.(name) = circshift(tc.cer.(name), 1, 2); end
            Rz = [cosd(1) -sind(1) 0; sind(1) cosd(1) 0; 0 0 1];
            t2 = NonGravForces.env.earthRadiance("ceres", ep, [12 24], r, 'Ceres', c2, 'RotationItrsToGcrs', pagemtimes(R, Rz'));
            tc.verifyBoundAndRecord("ceres: Earth turned by one cell, fields shifted by one column", ...
                max(abs([t2.albedo(:) - t1.albedo(:); (t2.exitanceLw_W_m2(:) - t1.exitanceLw_W_m2(:)) / 300])), 1e-9, "-", ...
                "a longitude moved by 1 deg through atan2 and back into a weight: 1e-13 in the weight, times a field gradient of at most a few W/m^2 per cell");
            % The three fields are read independently: doubling SWup and SWdown together leaves the albedo, and the longwave.
            c3 = tc.cer; c3.rsut = 2 * c3.rsut; c3.rsdt = 2 * c3.rsdt;
            t3 = NonGravForces.env.earthRadiance("ceres", ep, [12 24], r, 'Ceres', c3, 'RotationItrsToGcrs', R);
            tc.verifyBoundAndRecord("ceres albedo is invariant to a common scale of SWup and SWdown", max(abs(t3.albedo(:) - t1.albedo(:))), 1e-14, "-", "an exact power of two in numerator and denominator");
            tc.verifyEqual(t3.exitanceLw_W_m2, t1.exitanceLw_W_m2, 'the longwave does not read the shortwave fields');
        end
        function earthRadianceMetamorphic(tc)
            % One period later the seasonal term repeats (1e-14); mirrored latitudes: the odd term cancels in the sum.
            r = [6.88e6; 0; 0];
            ep0 = referenceSystem.datetimeToEpoch(datetime(1981, 12, 22) + days([0 365.25]), "UTC");
            t = NonGravForces.env.earthRadiance("knocke", ep0, [12 24], [r r]);
            tc.verifyBoundAndRecord("Knocke annual periodicity", max(abs(t.albedo(:, 1) - t.albedo(:, 2))), 1e-14, "-", "cos of 2 pi");
            epN = referenceSystem.datetimeToEpoch(datetime(1981, 12, 22) + days(50), "UTC");
            tN = NonGravForces.env.earthRadiance("knocke", epN, [12 24], [0; 0; 7e6]); tS = NonGravForces.env.earthRadiance("knocke", epN, [12 24], [0; 0; -7e6]);
            tc.verifyLessThan(max(abs(tN.lat_rad + tS.lat_rad)), 1e-12, 'mirrored caps');
            co = tN.coefficients; even = 2 * (co.a0 + co.a2 * 0.5 * (3 * sin(tN.lat_rad).^2 - 1));
            tc.verifyBoundAndRecord("Knocke odd term cancels over mirrored caps", max(abs(tN.albedo + tS.albedo - even)), 1e-14, "-", "");
        end
        function radiationSourcesMetamorphic(tc)
            % Rotating the attitude rotates every direction and leaves every pressure (bitwise);
            % scaling tsi scales every pressure by the same factor.
            [ep, r, q] = tc.sampleState();
            ep5 = referenceSystem.daySecondToEpoch([ep.day(1:5) ep.sec(1:5)], "GPS"); r5 = r(:, 1:5); q5 = q(:, 1:5);
            qz = [cos(0.3); 0; 0; sin(0.3)];
            s1 = NonGravForces.env.radiationSources(ep5, r5, q5, NonGravForces.Options("standard"));
            s2 = NonGravForces.env.radiationSources(ep5, r5, referenceSystem.multiplyQuaternions(q5, qz), NonGravForces.Options("standard"));
            tc.verifyEqual(s2.P, s1.P);
            % convention-free: with R(q) BODY->GCRS, dir_gcrs = R1 dir1 and dir2 = R2' dir_gcrs
            R1 = referenceSystem.rotationFromQuaternion(q5); R2 = referenceSystem.rotationFromQuaternion(referenceSystem.multiplyQuaternions(q5, qz));
            expect = pagemtimes(R2, 'transpose', pagemtimes(R1, s1.dir), 'none');
            tc.verifyBoundAndRecord("source directions rotate with the attitude", max(abs(expect - s2.dir), [], 'all'), 1e-14, "-", "two rotations");
            % every pressure is linear in tsi: the ratio of two runs is the ratio of the two tsi values (1e-14: one division)
            s3 = NonGravForces.env.radiationSources(ep5, r5, q5, NonGravForces.Options("standard", "srp.tsi", 1400));
            nz = s1.P > 0;
            tc.verifyBoundAndRecord("pressures scale with tsi", max(abs(s3.P(nz) ./ s1.P(nz) - 1400 / 1361)), 1e-14, "-", "one division");
            % grid refinement moves the flux monotonically towards the oracle (first-order convergence)
            K = tc.knocke; k = find(~K.uniform, 1);
            err = zeros(1, 3); g = {[6 12], [12 24], [24 48]};
            for i = 1:3, [~, Ei] = tc.knockeCase(k, g{i}, "knocke"); err(i) = abs(Ei - K.E_ir_W_m2(k)); end
            tc.verifyTrue(all(diff(err) < 0), 'error decreases with the grid');
            tc.verifyBoundAndRecord("ERP convergence ratio (err 12x24 / err 24x48)", err(2) / err(3), 6, "-", "first order expects ~4; measured 4.1");
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function windSpeed(tc)
            % Its whole cost is the frame rotation it delegates, plus two page products.
            [ep, r, ~, eop, v] = tc.sampleState();
            base = @() referenceSystem.rotationItrsToGcrs(ep, 'Eop', eop);
            fcn = @() NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop);
            tc.verifyBoundAndRecord("wind vs the rotationItrsToGcrs it calls", tc.speedRatio(fcn, base), 4.5, "-", "budget 3x the expected 1.4: one rotation call plus two page products; 5 repeats, median");
        end
        function windHwm14Speed(tc)
            % HWM14 over the co-rotation for the 600-epoch arc: the geodetic conversion and 600 model calls on top of the
            % same rotation. Measured see last_run; budget 3x the 2026-09-15 figure.
            [ep, r, ~, eop, v] = tc.sampleState();
            q = tc.speedRatio(@() NonGravForces.env.wind("hwm14", ep, r, v, 'Eop', eop, 'Ap', 15), @() NonGravForces.env.wind("corotation", ep, r, v, 'Eop', eop));
            tc.verifyBoundAndRecord("wind hwm14 / corotation time ratio, 600 epochs", q, 6, "-", "measured 1.76 on 2026-09-15; budget 3x");
        end
        function densitySpeed(tc)
            % Linear scaling of the MEX loop: 86400 epochs cost at most 2x ten times 8640 epochs (per-call overhead amortised).
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            epBig = referenceSystem.daySecondToEpoch([repmat(58849, 86400, 1) (0:86399)'], "UTC"); llaBig = [zeros(2, 86400); 490e3 + 1e3 * sin((1:86400) / 1000)];
            epSmall = referenceSystem.daySecondToEpoch([repmat(58849, 8640, 1) (0:8639)'], "UTC");
            q = tc.speedRatio(@() NonGravForces.env.density("msis00", epBig, llaBig, sw), @() NonGravForces.env.density("msis00", epSmall, llaBig(:, 1:8640), sw));
            tc.verifyBoundAndRecord("msis00 86400 / 8640 epochs time ratio", q / 10, 2, "-", "linear scaling; measured see last_run");
            % NRLMSIS 2.1 against NRLMSISE-00 on the same 8640 epochs: measured 4.8 (19 us per epoch) on 2026-09-15; budget 3x.
            q21 = tc.speedRatio(@() NonGravForces.env.density("msis21", epSmall, llaBig(:, 1:8640), sw), @() NonGravForces.env.density("msis00", epSmall, llaBig(:, 1:8640), sw));
            tc.verifyBoundAndRecord("msis21 / msis00 time ratio, 8640 epochs", q21, 15, "-", "measured 4.8 on 2026-09-15; budget 3x");
        end
        function magneticFieldSpeed(tc)
            % Degree 13 over degree 1 for 86400 positions: 104 terms against 3 plus the shared Legendre setup.
            % Measured see last_run; budget 3x the 2026-09-15 figure.
            psi = linspace(0, 2 * pi, 86400); r = 6861.2e3 * [cos(psi) .* cos(0.3 * psi); sin(psi) .* cos(0.3 * psi); sin(0.3 * psi)];
            q = tc.speedRatio(@() NonGravForces.env.magneticField(tc.igrf, r), @() NonGravForces.env.magneticField(tc.igrf, r, 'MaxDegree', 1));
            tc.verifyBoundAndRecord("magneticField degree 13 / degree 1 time ratio, 86400 positions", q, 150, "-", "measured 50.8 on 2026-09-15; budget 3x");
        end
        function eclipseFactorSpeed(tc)
            % Conical over cylindrical for 86400 epochs: two asin, an acos and the lens area; measured see last_run. Bound 8.
            psi = linspace(0, pi, 86400); r = (6378137 + 490e3) * [-cos(psi); zeros(1, 86400); sin(psi)]; s = [tc.au; 0; 0];
            q = tc.speedRatio(@() NonGravForces.env.eclipseFactor(r, s), @() NonGravForces.env.eclipseFactor(r, s, 'Model', "cylindrical"));
            tc.verifyBoundAndRecord("eclipse conical / cylindrical time ratio", q, 8, "-", "");
        end
        function earthRadianceClimatologySpeed(tc)
            % Against knocke on the same cap and epochs: the rotation, a latitude SEARCH per tile and
            % 8 gathers. The search is what makes this dearer than ceres's arithmetic indexing.
            clm = tc.syntheticClimatology(@(LAT, LON, m) 0.3 + 1e-4 * LON, @(LAT, LON, m) 220 + 0.1 * LAT, 180, 360, [58800 58900 59000]);
            [ep, r] = tc.sampleState();
            R = tc.eyes(numel(ep.sec));
            q = tc.speedRatio(@() NonGravForces.env.earthRadiance("ceresClimatology", ep, [12 24], r, 'Climatology', clm, 'RotationItrsToGcrs', R), ...
                @() NonGravForces.env.earthRadiance("knocke", ep, [12 24], r));
            tc.verifyBoundAndRecord("earthRadiance ceresClimatology / knocke time ratio [12 24], 600 epochs", q, 30, "-", "see last_run; budget set above the measured figure, ratio not absolute");
        end
        function earthRadianceCeresSpeed(tc)
            % The CERES tiles over Knocke's on the same grid and 600 epochs: the rotation, the geodetic latitude and 8 x 3
            % gathers per tile. Measured see last_run; budget 3x the 2026-09-15 figure.
            [ep, r] = tc.sampleState();
            R = tc.eyes(numel(ep.sec));
            q = tc.speedRatio(@() NonGravForces.env.earthRadiance("ceres", ep, [12 24], r, 'Ceres', tc.cer, 'RotationItrsToGcrs', R), ...
                @() NonGravForces.env.earthRadiance("knocke", ep, [12 24], r));
            tc.verifyBoundAndRecord("earthRadiance ceres / knocke time ratio [12 24], 600 epochs", q, 15, "-", "measured 5.0 on 2026-09-15; budget 3x");
        end
        function earthRadianceSpeed(tc)
            % Tile count scales the cost: [24 48] over [12 24] for 600 epochs, four times the tiles, bound 8.
            [ep, r] = tc.sampleState();
            q = tc.speedRatio(@() NonGravForces.env.earthRadiance("knocke", ep, [24 48], r), @() NonGravForces.env.earthRadiance("knocke", ep, [12 24], r));
            tc.verifyBoundAndRecord("earthRadiance [24 48] / [12 24] time ratio", q, 8, "-", "4x the tiles");
        end
        function radiationSourcesSpeed(tc)
            % The full source list over the bare tiling for the 600-epoch sample: pressures, visibility, rotation. Bound 6.
            [ep, r, q] = tc.sampleState(); opt = NonGravForces.Options("standard");
            ratio = tc.speedRatio(@() NonGravForces.env.radiationSources(ep, r, q, opt), @() NonGravForces.env.earthRadiance("knocke", ep, [12 24], r));
            tc.verifyBoundAndRecord("radiationSources / earthRadiance time ratio (600 epochs)", ratio, 6, "-", "");
        end
    end
    methods
        function [S, files, sw] = swarmScale(tc, varargin)
            %SWARMSCALE env.densityScale over the committed Swarm A/B/C days (2020-01-01..03), msis00.
            files = string(fullfile(tc.root, 'data', 'samples', 'swarm', "SW_OPER_DNS" + ["A" "B" "C"] + "POD_2_2020-01-01_2020-01-03.csv"));
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            S = NonGravForces.env.densityScale(files, sw, varargin{:});
        end
    end
    methods (Test, TestTags = {'contract'})
        % densityScale (O-19, 2026-10-27)
        function densityScaleContract(tc)
            [S, files, sw] = tc.swarmScale();
            tc.verifyEqual(S.day, [58849 58850 58851]); tc.verifyEqual(S.satellite, ["A" "B" "C"]);
            tc.verifySize(S.factor, [3 3]); tc.verifySize(S.height_m, [3 3]); tc.verifyTrue(all(isfinite(S.factor), 'all'));
            tc.verifyEqual(S.model, "msis00"); tc.verifyEqual(S.stride, 10);
            e = referenceSystem.daySecondToEpoch([58849 0; 58850 43200; 58851 86399], "UTC");
            f = S.at(e, [490e3 490e3 490e3]);
            tc.verifySize(f, [1 3]); tc.verifyTrue(all(f > 0 & isfinite(f)));
            tc.verifyErrorId(@() NonGravForces.env.densityScale(files(2), sw), 'NonGravForces:densityScale:inputSatellites');
            tc.verifyErrorId(@() NonGravForces.env.densityScale(files([1 3]), sw), 'NonGravForces:densityScale:inputSatellites');
            tc.verifyErrorId(@() NonGravForces.env.densityScale(files, sw, 'Model', "exponential"), 'NonGravForces:densityScale:inputModel');
            tc.verifyErrorId(@() S.at(referenceSystem.mjdToEpoch(58852.5, "UTC"), 490e3), 'NonGravForces:densityScale:dataDay');
            tc.verifyErrorId(@() S.at(referenceSystem.mjdToEpoch(58849.5, "UTC"), 900e3), 'NonGravForces:densityScale:inputRange');
            tc.verifyErrorId(@() S.at(e, 490e3), 'NonGravForces:densityScale:inputRange');
        end
    end
    methods (Test, TestTags = {'reference'})
        function densityScaleReference(tc)
            % Oracle: an INDEPENDENT density model. CNES's DTM2020 (fitted to accelerometer densities, among them
            % Swarm A's) against NRLMSISE-00 at Swarm A's own positions must show the same ratio ESA's POD densities
            % show, to the 10-30 % both are uncertain by; the month gave 0.53 (DTM2020) against 0.49 (Swarm).
            [S, files, sw] = tc.swarmScale();
            [ep, lla, ~, aux] = NonGravForces.io.readSwarmDensity(files(1));
            k = find(aux.valid & ep.day(:)' == 58849); k = k(1:10:end);
            e = referenceSystem.daySecondToEpoch([ep.day(k) ep.sec(k)], "UTC");
            rD = NonGravForces.env.density("dtm2020", e, lla(:, k), sw);
            rM = NonGravForces.env.density("msis00", e, lla(:, k), sw);
            fD = sum(rD) / sum(rM);
            tc.verifyBoundAndRecord("Swarm A / msis00 vs DTM2020 / msis00 on 2020-01-01, |log|", abs(log(S.factor(1, 1) / fD)), log(1.3), "-", "two independent sources of the same ratio");
            % and the ratio itself is the solar-minimum excess of NRLMSISE-00 the month measured (0.41..0.61 per day)
            tc.verifyGreaterThan(min(S.factor, [], 'all'), 0.3); tc.verifyLessThan(max(S.factor, [], 'all'), 0.8);
        end
    end
    methods (Test, TestTags = {'precision'})
        function densityScalePrecision(tc)
            % At the bracket's own heights the interpolation returns the bracket's factors: an identity up to
            % exp(log(x)) rounding, a few ulp.
            S = tc.swarmScale();
            e = referenceSystem.daySecondToEpoch([S.day' repmat(43200, 3, 1)], "UTC");
            fl = S.at(e, S.hLow_m); fh = S.at(e, S.hHigh_m);
            err = max(abs([fl ./ S.low, fh ./ S.high] - 1));
            tc.verifyBoundAndRecord("densityScale.at at the bracket heights / bracket factors - 1", err, 8 * eps, "-", "exp(log(x)) round trip");
            tc.verifyEqual(S.low, sqrt(S.factor(1, :) .* S.factor(3, :)), 'RelTol', 8 * eps, 'the lower factor is the geometric mean of A and C');
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function densityScaleVectorization(tc)
            % S.at is element-wise: 300 epochs over three days at once equal 300 single calls, bitwise.
            S = tc.swarmScale();
            t = linspace(58849, 58851.999, 300);
            e = referenceSystem.daySecondToEpoch([floor(t') round((t' - floor(t')) * 86400)], "UTC");
            h = 470e3 + 30e3 * sin(1:300);
            f = S.at(e, h);
            g = zeros(1, 300);
            for k = 1:300
                g(k) = S.at(referenceSystem.daySecondToEpoch([e.day(k) e.sec(k)], "UTC"), h(k));
            end
            tc.verifyEqual(f, g);
        end
    end
    methods (Test, TestTags = {'robustness'})
        function densityScaleRobustness(tc)
            % A day with too few records is refused, not filled; one lower satellite suffices, and is then the bracket.
            [~, files, sw] = tc.swarmScale();
            S = NonGravForces.env.densityScale(files, sw, 'MinRecords', 2000);      % 288 records a day at Stride 10
            tc.verifyTrue(all(isnan(S.factor), 'all'));
            tc.verifyErrorId(@() S.at(referenceSystem.mjdToEpoch(58849.5, "UTC"), 490e3), 'NonGravForces:densityScale:dataDay');
            S2 = NonGravForces.env.densityScale(files([1 2]), sw);
            tc.verifyEqual(S2.satellite, ["A" "B"]); tc.verifyEqual(S2.low, S2.factor(1, :), 'RelTol', 4 * eps);
            % an epoch in GPS time is put on its UTC day: 2020-01-02T00:00:10 GPS is still 2020-01-01 in UTC
            Sall = tc.swarmScale();
            eG = referenceSystem.daySecondToEpoch([58850 10], "GPS");
            eU = referenceSystem.daySecondToEpoch([58849 86392], "UTC");
            tc.verifyEqual(Sall.at(eG, 490e3), Sall.at(eU, 490e3));
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function densityScaleMetamorphic(tc)
            % log f is linear in height: the midpoint of the bracket gets the geometric mean of its ends; the order of
            % the files changes nothing; a denser stride moves the daily factor by less than the product's own noise.
            [S, files, sw] = tc.swarmScale();
            e = referenceSystem.daySecondToEpoch([S.day' repmat(43200, 3, 1)], "UTC");
            fm = S.at(e, (S.hLow_m + S.hHigh_m) / 2);
            tc.verifyEqual(fm, sqrt(S.low .* S.high), 'RelTol', 1e-12);
            Sr = NonGravForces.env.densityScale(files([3 1 2]), sw);
            tc.verifyEqual(Sr.factor, S.factor);
            S1 = NonGravForces.env.densityScale(files, sw, 'Stride', 2);
            tc.verifyBoundAndRecord("densityScale Stride 2 vs 10, max |log ratio|", max(abs(log(S1.factor ./ S.factor)), [], 'all'), 0.05, "-", "5-min sampling of a ~20-min product");
        end
    end
    methods (Test, TestTags = {'speed'})
        function densityScaleSpeed(tc)
            % Three satellites x three days: three reads and 2592 msis00 points. Against the same number of msis00
            % points alone; the reads and the bookkeeping are what the ratio measures. Bound 20, first run.
            [~, files, sw] = tc.swarmScale();
            e = referenceSystem.daySecondToEpoch([repmat(58849, 2592, 1) (0:2591)' * 30], "UTC");
            lla = [zeros(2, 2592); 480e3 * ones(1, 2592)];
            q = tc.speedRatio(@() NonGravForces.env.densityScale(files, sw), @() NonGravForces.env.density("msis00", e, lla, sw));
            tc.verifyBoundAndRecord("densityScale / msis00 on as many points, time ratio", q, 20, "-", "first run");
        end
    end
    methods (Test, TestTags = {'contract'})
        % densityScale with 'Fallback': the combined scale (O-19, 2026-10-28)
        function densityScaleFallbackContract(tc)
            D = tc.driverScale();
            [C, files, sw] = tc.swarmScale('Fallback', D);
            tc.verifyEqual(C.kind, "combined"); tc.verifyEqual(C.fallback.kind, "driver");
            S = tc.swarmScale();
            tc.verifyEqual(S.kind, "swarm"); tc.verifyEmpty(S.fallback);
            % a fallback that is not a driver scale, or was formed for another model, is refused
            tc.verifyErrorId(@() NonGravForces.env.densityScale(files, sw, 'Fallback', S), 'NonGravForces:densityScale:inputFallback');
            tc.verifyErrorId(@() NonGravForces.env.densityScale(files, sw, 'Fallback', struct('x', 1)), 'NonGravForces:densityScale:inputFallback');
            tc.verifyErrorId(@() NonGravForces.env.densityScale(files, sw, 'Model', "msis21", 'Fallback', D), 'NonGravForces:densityScale:inputFallback');
            % the day the plain scale refuses (2020-01-05) the combined one takes from the fallback, and says so
            e = referenceSystem.daySecondToEpoch([58849 43200; 58853 43200], "UTC");
            tc.verifyErrorId(@() S.at(e, [490e3 490e3]), 'NonGravForces:densityScale:dataDay');
            f = C.at(e, [490e3 490e3]);
            tc.verifySize(f, [1 2]); tc.verifyTrue(all(f > 0 & isfinite(f)));
            tc.verifyEqual(C.source(e), ["swarm" "driver"]);
            tc.verifyEqual(S.source(e), ["swarm" "none"]);
        end
    end
    methods (Test, TestTags = {'reference'})
        function densityScaleFallbackReference(tc)
            % Oracle: the two scales it combines. On Swarm's days it IS the plain Swarm scale, on every other day it IS
            % the driver scale, bitwise -- epochs of both kinds interleaved in one call.
            D = tc.driverScale(); C = tc.swarmScale('Fallback', D); S = tc.swarmScale();
            e = referenceSystem.daySecondToEpoch([58853 3600; 58849 3600; 58854 70000; 58850 50000; 58851 100], "UTC");
            h = [480e3 485e3 490e3 495e3 500e3];
            f = C.at(e, h);
            iS = [2 4 5]; iD = [1 3];
            tc.verifyEqual(f(iS), S.at(referenceSystem.daySecondToEpoch([e.day(iS) e.sec(iS)], "UTC"), h(iS)));
            tc.verifyEqual(f(iD), D.at(referenceSystem.daySecondToEpoch([e.day(iD) e.sec(iD)], "UTC"), h(iD)));
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function densityScaleFallbackVectorization(tc)
            % 300 epochs over six days, three with Swarm and three without, at once equal 300 single calls, bitwise.
            C = tc.swarmScale('Fallback', tc.driverScale());
            t = linspace(58849, 58854.999, 300);
            e = referenceSystem.daySecondToEpoch([floor(t') round((t' - floor(t')) * 86400)], "UTC");
            h = 470e3 + 30e3 * sin(1:300);
            f = C.at(e, h);
            g = zeros(1, 300);
            for k = 1:300
                g(k) = C.at(referenceSystem.daySecondToEpoch([e.day(k) e.sec(k)], "UTC"), h(k));
            end
            tc.verifyEqual(f, g);
        end
    end
    methods (Test, TestTags = {'robustness'})
        function densityScaleFallbackRobustness(tc)
            % A Swarm day with too few records hands the whole day to the fallback, never half of it; a day neither
            % covers (2019-10-01: F30's 81 days before it are not in the excerpt) is refused by the fallback.
            D = tc.driverScale();
            [~, files, sw] = tc.swarmScale();
            C = NonGravForces.env.densityScale(files, sw, 'MinRecords', 2000, 'Fallback', D);
            e = referenceSystem.daySecondToEpoch([58849 600; 58849 80000], "UTC");
            tc.verifyEqual(C.source(e), ["driver" "driver"]);
            tc.verifyEqual(C.at(e, [490e3 495e3]), D.at(e, [490e3 495e3]));
            tc.verifyErrorId(@() C.at(referenceSystem.mjdToEpoch(58757.5, "UTC"), 490e3), 'NonGravForces:densityDriverScale:dataDay');
        end
    end
    methods
        function [S, sw, cls] = driverScale(tc, varargin)
            %DRIVERSCALE env.densityDriverScale over the committed space-weather excerpts (2019-09..2020-03).
            sw = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
            cls = NonGravForces.io.spaceWeather(fullfile(tc.root, 'data', 'spaceweather', 'radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt'));
            S = NonGravForces.env.densityDriverScale(sw, cls, varargin{:});
        end
    end
    methods (Test, TestTags = {'contract'})
        % densityDriverScale (O-19 step 4, 2026-10-28)
        function densityDriverScaleContract(tc)
            [S, sw, cls] = tc.driverScale();
            tc.verifyEqual([S.kind S.model S.proxy], ["driver" "msis00" "F30"]); tc.verifyEqual(S.offset, 0);
            e = referenceSystem.daySecondToEpoch([58849 0; 58850 43200; 58851 86399], "UTC");
            f = S.at(e, [490e3 490e3 490e3]);
            tc.verifySize(f, [1 3]); tc.verifyTrue(all(f > 0.2 & f < 2));
            tc.verifySize(S.drivers(e), [8 3]);
            tc.verifyErrorId(@() NonGravForces.env.densityDriverScale(sw, cls, 'Model', "msis21"), 'NonGravForces:densityDriverScale:inputModel');
            tc.verifyErrorId(@() NonGravForces.env.densityDriverScale(cls, sw), 'NonGravForces:densityDriverScale:inputFamily');
            % 2019-10-01 needs F30 from 2019-07-12; the excerpt starts 2019-09-01
            tc.verifyErrorId(@() S.at(referenceSystem.mjdToEpoch(58757.5, "UTC"), 490e3), 'NonGravForces:densityDriverScale:dataDay');
            tc.verifyErrorId(@() S.at(e, 490e3), 'NonGravForces:densityDriverScale:inputRange');
        end
    end
    methods (Test, TestTags = {'reference'})
        function densityDriverScaleReference(tc)
            % Oracle: tools/densityDriverFit.m's own drivers and prediction for 2020-01-01 on the runner's full
            % pool (run 40676249479, F30 coefficient record). Independent code -- the tool builds Fbar with a
            % trailing movmean over the file's rows, the unit with a calendar window -- over the same bytes.
            S = tc.driverScale();
            zRef = [44.592592592592595 1.2074074074074019 2 0 1 0 1 0.096843864597306581];
            logfRef = -0.54833027703885062;
            e = referenceSystem.mjdToEpoch(58849.5, "UTC");
            z = S.drivers(e);
            tc.verifyBoundAndRecord("densityDriverScale drivers vs densityDriverFit's on 2020-01-01, max rel", max(abs(z' - zRef) ./ max(abs(zRef), 1)), 1e-12, "-", "an 81-term mean summed in another order");
            tc.verifyBoundAndRecord("densityDriverScale log f vs densityDriverFit's on 2020-01-01", abs(log(S.at(e, 490e3)) - logfRef), 1e-12, "-", "transcribed to 17 digits");
        end
    end
    methods (Test, TestTags = {'precision'})
        function densityDriverScalePrecision(tc)
            % 'Offset' adds to log f and nothing else: an identity up to exp/log rounding of a number near 0.3.
            S = tc.driverScale(); So = tc.driverScale('Offset', -0.041);
            e = referenceSystem.daySecondToEpoch([(58849:58851)' repmat(43200, 3, 1)], "UTC");
            h = [490e3 490e3 490e3];
            err = max(abs((log(So.at(e, h)) - log(S.at(e, h))) + 0.041));
            tc.verifyBoundAndRecord("densityDriverScale Offset -0.041: log ratio + 0.041", err, 1e-14, "-", "exp then log of numbers below 2");
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function densityDriverScaleVectorization(tc)
            % S.at over 300 epochs on three days at once equals 300 single calls, bitwise.
            S = tc.driverScale();
            t = linspace(58849, 58851.999, 300);
            e = referenceSystem.daySecondToEpoch([floor(t') round((t' - floor(t')) * 86400)], "UTC");
            h = 470e3 + 30e3 * sin(1:300);
            f = S.at(e, h);
            g = zeros(1, 300);
            for k = 1:300
                g(k) = S.at(referenceSystem.daySecondToEpoch([e.day(k) e.sec(k)], "UTC"), h(k));
            end
            tc.verifyEqual(f, g);
        end
    end
    methods (Test, TestTags = {'robustness'})
        function densityDriverScaleRobustness(tc)
            % A missing F30 day inside the 81-day window or a missing Ap is refused, not averaged over; an epoch
            % in GPS time is put on its UTC day.
            [S, sw, cls] = tc.driverScale();
            c2 = cls; c2.index.f30(cls.mjd == 58800) = NaN;
            S2 = NonGravForces.env.densityDriverScale(sw, c2);
            tc.verifyErrorId(@() S2.at(referenceSystem.mjdToEpoch(58849.5, "UTC"), 490e3), 'NonGravForces:densityDriverScale:dataDay');
            s2 = sw; s2.index.apDaily(sw.mjd == 58849) = NaN;
            S3 = NonGravForces.env.densityDriverScale(s2, cls);
            tc.verifyErrorId(@() S3.at(referenceSystem.mjdToEpoch(58849.5, "UTC"), 490e3), 'NonGravForces:densityDriverScale:dataDay');
            tc.verifyEqual(S3.at(referenceSystem.mjdToEpoch(58850.5, "UTC"), 490e3), S.at(referenceSystem.mjdToEpoch(58850.5, "UTC"), 490e3));
            eG = referenceSystem.daySecondToEpoch([58850 10], "GPS");
            eU = referenceSystem.daySecondToEpoch([58849 86392], "UTC");
            tc.verifyEqual(S.at(eG, 490e3), S.at(eU, 490e3));
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function densityDriverScaleMetamorphic(tc)
            % The factor is the same at every height of a day, and raising the day's Ap by 10 multiplies it by
            % exp(10 c_ap) exactly; the day before's Ap does not enter.
            [S, sw, cls] = tc.driverScale();
            e = referenceSystem.mjdToEpoch(58849.5, "UTC");
            tc.verifyEqual(S.at(e, 300e3), S.at(e, 600e3));
            s2 = sw; j = sw.mjd == 58849; s2.index.apDaily(j) = sw.index.apDaily(j) + 10;
            S2 = NonGravForces.env.densityDriverScale(s2, cls);
            tc.verifyEqual(S2.at(e, 490e3) / S.at(e, 490e3), exp(10 * S.coef.ap), 'RelTol', 1e-13);
            s3 = sw; s3.index.apDaily(sw.mjd == 58848) = 400;
            tc.verifyEqual(NonGravForces.env.densityDriverScale(s3, cls).at(e, 490e3), S.at(e, 490e3));
        end
    end
    methods (Test, TestTags = {'speed'})
        function densityDriverScaleSpeed(tc)
            % A day of 10 s epochs: one driver evaluation per UTC day, then one exp per epoch. Against msis00 on
            % as many points, which the factor multiplies anyway. Bound 1, first run.
            S = tc.driverScale();
            e = referenceSystem.daySecondToEpoch([repmat(58849, 8640, 1) (0:8639)' * 10], "UTC");
            lla = [zeros(2, 8640); 490e3 * ones(1, 8640)];
            [~, sw] = tc.driverScale();
            q = tc.speedRatio(@() S.at(e, lla(3, :)), @() NonGravForces.env.density("msis00", e, lla, sw));
            tc.verifyBoundAndRecord("densityDriverScale.at / msis00 on 8640 points, time ratio", q, 1, "-", "first run");
        end
    end
end
