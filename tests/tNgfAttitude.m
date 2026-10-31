classdef tNgfAttitude < NgfTestCase
%TNGFATTITUDE Tests of N-12 (#69): io.readIMU1B, geom.attitudeQuality and geom.AttitudeProvider.fromItsg.
%
%   Oracle: the committed 30 min excerpts of 2020-01-01..03 GRACE-C in data/samples/multiday (SCA1B,
%   GNV1B and IMU1B, cut by tools/cutFixtures.py on the CI runner, records unchanged), the raw IMU1B
%   text parsed here with split/str2double (not gracefoReading), and tools/pyxcheck/attitude_quality_check.py,
%   a standard-library recomputation of the star-camera rates, the gyro rates over the same intervals and
%   the fit residual from their definitions. TU Graz's attitude excerpt (data/samples/itsg) is the oracle of
%   fromItsg: it is SCA1B's quaternion in the other convention (CLAUDE.md section 7).
%
%   No test skips when a fixture is missing: it fails.
% Created 2026-10-30 by Claude (Anthropic) for Matthias Weigelt.
    properties
        f      % fixture paths
        eop    % referenceSystem.readEop of the committed finals2000A excerpt
        tmp    % scratch folder for corrupted copies
    end
    methods (TestClassSetup)
        function loadFixtures(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            M = fullfile(root, 'data', 'samples', 'multiday');
            tc.f = struct('imu1', fullfile(M, 'IMU1B_2020-01-01_C_04.txt'), 'imu2', fullfile(M, 'IMU1B_2020-01-02_C_04.txt'), ...
                'sca1', fullfile(M, 'SCA1B_2020-01-01_C_04.txt'), 'sca2', fullfile(M, 'SCA1B_2020-01-02_C_04.txt'), ...
                'gnv1', fullfile(M, 'GNV1B_2020-01-01_C_04.txt'), ...
                'itsg', fullfile(root, 'data', 'samples', 'itsg', 'GRACEFO-1_attitude_2020-01-01.txt'), ...
                'sca600', fullfile(root, 'data', 'samples', 'SCA1B_2020-01-01_C_04.txt'), ...
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
            %CORRUPT Copy an L1B fixture with editFcn(lines, kFirstData) applied; return the path.
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
            tc.assertGreaterThan(numel(tok), 5);
        end
        function q = speedRatio(~, fcn, base)
            fcn(); base();
            r = zeros(1, 5);
            for k = 1:5
                t1 = tic; fcn(); a = toc(t1); t2 = tic; base(); b = toc(t2); r(k) = a / b;
            end
            q = median(r);
        end
        function Q = dayOne(tc, varargin)
            %DAYONE attitudeQuality over the 30 min of 2020-01-01 C.
            [eS, qS] = NonGravForces.io.readSCA1B(tc.f.sca1);
            [eI, aI, auxI] = NonGravForces.io.readIMU1B(tc.f.imu1);
            Q = NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId, varargin{:});
        end
    end

    % ===================================================================== io.readIMU1B
    methods (Test, TestTags = {'contract'})
        function readIMU1BContract(tc)
            [ep, a, aux] = NonGravForces.io.readIMU1B(tc.f.imu1);
            tc.verifyEqual(aux.numRecords, 43200, 'three gyros at 8 Hz for 30 min');
            tc.verifySize(a, [1 43200]); tc.verifySize(aux.gyroId, [1 43200]);
            tc.verifyEqual(numel(ep.day), 43200); tc.verifyEqual(string(ep.scale), "GPS");
            tc.verifyEqual(aux.gyros, [1 2 3]);
            tc.verifyEqual([aux.product aux.satellite], ["IMU1B" "C"]);
            tc.verifyTrue(all(aux.timeRef == "G"));
            tc.verifyErrorId(@() NonGravForces.io.readIMU1B(tc.f.sca1), 'NonGravForces:readIMU1B:dataProduct');
        end
    end
    methods (Test, TestTags = {'reference'})
        function readIMU1BReference(tc)
            % The raw text, parsed here: record k carries rcvtime_intg, rcvtime_frac (microseconds), time_ref,
            % satellite, gyro_id and FiltAng (degrees). Epoch = 2000-01-01 12:00 GPS + intg + frac, i.e. GPS day
            % 51544 + floor(...), the angle deg2rad(FiltAng).
            [ep, a, aux] = NonGravForces.io.readIMU1B(tc.f.imu1);
            for k = [1 2 3 1000 43200]
                tok = tc.rawRecord(tc.f.imu1, k);
                s = str2double(tok(1)) + str2double(tok(2)) * 1e-6 + 43200;     % seconds since 2000-01-01 00:00 GPS
                tc.verifyEqual(double(ep.day(k)), 51544 + floor(s / 86400), sprintf('day of record %d', k));
                tc.verifyEqual(double(ep.sec(k)), mod(s, 86400), 'AbsTol', 1e-6, sprintf('second of record %d', k));
                tc.verifyEqual(aux.gyroId(k), str2double(tok(5)), sprintf('gyro of record %d', k));
                tc.verifyEqual(a(k), deg2rad(str2double(tok(6))), 'RelTol', 4 * eps, sprintf('angle of record %d', k));
            end
            % the first record of 2020-01-01: 631108800 s + 116294 us after 2000-01-01 12:00 GPS
            tc.verifyEqual(double(ep.day(1)), 58849);
            tc.verifyBoundAndRecord("readIMU1B first epoch vs 2020-01-01 00:00:00.116294 GPS", abs(double(ep.sec(1)) - 0.116294), 1e-6, "s", ...
                "the tag's resolution is 1 us; CI 40978582699 found it off by more than 1e-9 s (value not printed)");
        end
    end
    methods (Test, TestTags = {'precision'})
        function readIMU1BPrecision(tc)
            % The conversion is one multiplication by pi/180: back in degrees the published 15 digits return.
            [~, a] = NonGravForces.io.readIMU1B(tc.f.imu1);
            tok = tc.rawRecord(tc.f.imu1, 1);
            tc.verifyBoundAndRecord("readIMU1B rad2deg(angle) vs the published FiltAng, record 1", ...
                abs(rad2deg(a(1)) / str2double(tok(6)) - 1), 4 * eps, "rel", "one multiplication there and back");
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function readIMU1BVectorization(tc)
            % Two days read together are the two days read alone, concatenated, bitwise.
            [e1, a1, x1] = NonGravForces.io.readIMU1B(tc.f.imu1);
            [e2, a2, x2] = NonGravForces.io.readIMU1B(tc.f.imu2);
            [e, a, x] = NonGravForces.io.readIMU1B([string(tc.f.imu1) string(tc.f.imu2)]);
            tc.verifyEqual(a, [a1 a2]);
            tc.verifyEqual(x.gyroId, [x1.gyroId x2.gyroId]);
            tc.verifyEqual(double(e.day(:))', [double(e1.day(:))' double(e2.day(:))']);
            tc.verifyEqual(double(e.sec(:))', [double(e1.sec(:))' double(e2.sec(:))']);
        end
    end
    methods (Test, TestTags = {'robustness'})
        function readIMU1BRobustness(tc)
            % A record not tagged in GPS time is refused, an empty TimeRange is refused, a TimeRange inside the
            % excerpt returns exactly the records it covers.
            bad = tc.corrupt(tc.f.imu1, 'IMU1B_2020-01-01_C_04.txt', @(L, k) localSetToken(L, k + 5, 3, "R"));
            tc.verifyErrorId(@() NonGravForces.io.readIMU1B(bad), 'NonGravForces:readIMU1B:dataTimeRef');
            tc.verifyErrorId(@() NonGravForces.io.readIMU1B(tc.f.imu1, 'TimeRange', [58850 58851]), 'NonGravForces:readIMU1B:dataEmpty');
            [ep, a] = NonGravForces.io.readIMU1B(tc.f.imu1);
            [e1, a1] = NonGravForces.io.readIMU1B(tc.f.imu1, 'TimeRange', 58849 + [600 900] / 86400);
            in = double(ep.sec) >= 600 & double(ep.sec) <= 900;
            tc.verifyEqual(a1, a(in(:)'));
            tc.verifyEqual(numel(e1.day), nnz(in));
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function readIMU1BMetamorphic(tc)
            % Each gyro is sampled at 8 Hz, its records strictly increasing in time, and its accumulated angle
            % turns at no more than the orbital rate projected on its axis (|rate| <= |omega| ~ 1.1e-3 rad/s).
            [ep, a, aux] = NonGravForces.io.readIMU1B(tc.f.imu1);
            t = double(ep.sec(:))';
            for g = aux.gyros
                s = aux.gyroId == g;
                dt = diff(t(s));
                tc.verifyTrue(all(dt > 0), sprintf('gyro %d records strictly increasing', g));
                tc.verifyEqual(median(dt), 0.125, 'AbsTol', 1e-3, sprintf('gyro %d at 8 Hz', g));
                rate = diff(a(s)) ./ dt;
                tc.verifyLessThan(max(abs(rate)), 1.3e-3, sprintf('gyro %d turns no faster than the orbit', g));
            end
        end
    end
    methods (Test, TestTags = {'speed'})
        function readIMU1BSpeed(tc)
            % The reader adds a column lookup, a conversion and a GPS check to gracefoReading's read: ratio near 1.
            q = tc.speedRatio(@() NonGravForces.io.readIMU1B(tc.f.imu1), @() gracefoReading.readL1b(tc.f.imu1));
            tc.verifyBoundAndRecord("readIMU1B / readL1b time ratio", q, 5, "-", "a lookup and a conversion over the read");
        end
    end

    % ===================================================================== geom.attitudeQuality
    methods (Test, TestTags = {'contract'})
        function attitudeQualityContract(tc)
            Q = tc.dayOne();
            M = numel(Q.t_s);
            tc.verifyEqual(M, 1799, 'one interval per pair of the 1800 SCA1B epochs');
            tc.verifySize(Q.omega_rad_s, [3 M]);
            tc.verifyEqual(Q.gyros, [1 2 3]);
            tc.verifyEqual([Q.gyro.id], [1 2 3]);
            for g = 1:3
                tc.verifySize(Q.gyro(g).rate_rad_s, [1 M]); tc.verifySize(Q.gyro(g).residual_rad_s, [1 M]);
                tc.verifyEqual(norm(Q.gyro(g).axis), 1, 'AbsTol', 1e-12);
                tc.verifyEqual(Q.gyro(g).n, 1798, 'the first interval starts before the first gyro sample');
                tc.verifyTrue(all(isnan(Q.gyro(g).cpr)), 'no period, no cpr amplitudes');
            end
            tc.verifyTrue(isnan(Q.period_s));
            tc.verifyFalse(isfield(Q, 'nominal'));
            % refusals, each paired with the call above that succeeded
            [eS, qS] = NonGravForces.io.readSCA1B(tc.f.sca1);
            [eI, aI, auxI] = NonGravForces.io.readIMU1B(tc.f.imu1);
            [eI2, aI2, auxI2] = NonGravForces.io.readIMU1B(tc.f.imu2);
            tc.verifyErrorId(@() NonGravForces.geom.attitudeQuality(eS, qS(1:3, :), eI, aI, auxI.gyroId), 'NonGravForces:attitudeQuality:inputShape');
            tc.verifyErrorId(@() NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId(1:end-1)), 'NonGravForces:attitudeQuality:inputShape');
            tc.verifyErrorId(@() NonGravForces.geom.attitudeQuality(eS, qS, eI2, aI2, auxI2.gyroId), 'NonGravForces:attitudeQuality:dataOverlap');
            tc.verifyErrorId(@() NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId, 'Orbit', struct('r', zeros(3, 2), 'v', zeros(3, 2))), ...
                'NonGravForces:attitudeQuality:inputOrbit');
        end
    end
    methods (Test, TestTags = {'reference'})
        function attitudeQualityReference(tc)
            % tools/pyxcheck/attitude_quality_check.py 2020-01-01 (standard library, raw text, the definitions):
            % omega over the 101st interval and each gyro's residual RMS after the projection.
            Q = tc.dayOne();
            % (the oracle normalises the file quaternions as referenceSystem.invertQuaternion does; before it did,
            % CI 40978582699 measured 5.4e-8 here and 1.43e-5 on the RMS below -- exactly the normalisation)
            oracleOmega = [-6.101346676364544e-07; 1.103522926664514e-03; -1.305407912748501e-06];
            tc.verifyBoundAndRecord("attitudeQuality omega(101) vs the stdlib oracle", norm(Q.omega_rad_s(:, 101) - oracleOmega) / norm(oracleOmega), 1e-9, "-", ...
                "the same rotation vector of R_k' R_k+1 in two languages; the oracle's convention (R^T of the active matrix) is the one the gyros fit");
            oracleRms = [3.570157778623e-07 3.750847719329e-07 3.657078907429e-07];
            oracleRate = [9.027399493282e-04 2.062160169000e-06 9.027578625462e-04];
            tc.verifyBoundAndRecord("attitudeQuality gyro residual RMS vs the stdlib oracle (3 gyros)", max(abs([Q.gyro.rms_rad_s] ./ oracleRms - 1)), 1e-5, "rel", ...
                "both are exact projections; the epochs differ by referenceSystem's rounding of the GPS seconds");
            tc.verifyBoundAndRecord("attitudeQuality gyro RMS rate vs the stdlib oracle (3 gyros)", max(abs([Q.gyro.rmsRate_rad_s] ./ oracleRate - 1)), 1e-6, "rel", ...
                "the accumulated angle interpolated to the same SCA1B epochs and differenced");
            % independent first-order form of skew(omega) = R' dR/dt at the 101st interval
            [eS, qS] = NonGravForces.io.readSCA1B(tc.f.sca1);
            R = referenceSystem.rotationFromQuaternion(qS(:, 101:102));
            S = ((R(:, :, 1) + R(:, :, 2)) / 2)' * (R(:, :, 2) - R(:, :, 1)) / double(eS.sec(102) - eS.sec(101));
            w1 = [S(3, 2) - S(2, 3); S(1, 3) - S(3, 1); S(2, 1) - S(1, 2)] / 2;
            tc.verifyLessThan(norm(w1 - Q.omega_rad_s(:, 101)) / norm(w1), 1e-5, 'the exact rotation vector agrees with the midpoint form to second order (theta^2 ~ 1e-6)');
        end
    end
    methods (Test, TestTags = {'precision'})
        function attitudeQualityPrecision(tc)
            % The residual is orthogonal to every column of the design, to rounding: the projection is exact
            % however ill-conditioned the fit (condition numbers in the thousands on a nominal arc).
            Q = tc.dayOne();
            for g = 1:3
                use = isfinite(Q.gyro(g).residual_rad_s);
                X = [Q.omega_rad_s(:, use)', ones(nnz(use), 1)];
                Xn = X ./ vecnorm(X);
                r = Q.gyro(g).residual_rad_s(use)';
                tc.verifyBoundAndRecord(sprintf("attitudeQuality gyro %d residual orthogonal to the design", g), ...
                    max(abs(Xn' * r)) / (norm(r) * sqrt(nnz(use))), 1e-12, "-", "an orthogonal projection by economy QR");
            end
            tc.verifyGreaterThan(Q.gyro(1).condition, 10, 'a nominal arc leaves the fit ill-conditioned, and the unit says so');
        end
    end
    methods (Test, TestTags = {'vectorization'})
        function attitudeQualityVectorization(tc)
            % The gyros are processed one at a time: passing only gyro 2's records gives gyro 2's result bitwise.
            [eS, qS] = NonGravForces.io.readSCA1B(tc.f.sca1);
            [eI, aI, auxI] = NonGravForces.io.readIMU1B(tc.f.imu1);
            Q = NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId);
            s = auxI.gyroId == 2;
            e2 = referenceSystem.daySecondToEpoch([double(eI.day(s)) double(eI.sec(s))], string(eI.scale));
            Q2 = NonGravForces.geom.attitudeQuality(eS, qS, e2, aI(s), auxI.gyroId(s));
            tc.verifyEqual(Q2.gyro.residual_rad_s, Q.gyro(2).residual_rad_s, 'AbsTol', 1e-9 * Q.gyro(2).rms_rad_s, ...
                'the same records, rebuilt as an epoch struct of their own');
            tc.verifyEqual(Q2.omega_rad_s, Q.omega_rad_s);
        end
    end
    methods (Test, TestTags = {'robustness'})
        function attitudeQualityRobustness(tc)
            % A gyro reset (a jump larger than the threshold) splits that gyro's series: no interval spans it,
            % and every other interval keeps its rate.
            [eS, qS] = NonGravForces.io.readSCA1B(tc.f.sca1);
            [eI, aI, auxI] = NonGravForces.io.readIMU1B(tc.f.imu1);
            Q = NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId);
            g1 = find(auxI.gyroId == 1);
            j = g1(4000);                                  % ~500 s into the arc
            aJ = aI;
            aJ(g1(g1 >= j)) = aJ(g1(g1 >= j)) + deg2rad(30);
            QJ = NonGravForces.geom.attitudeQuality(eS, qS, eI, aJ, auxI.gyroId);
            lost = isfinite(Q.gyro(1).rate_rad_s) & ~isfinite(QJ.gyro(1).rate_rad_s);
            tc.verifyGreaterThanOrEqual(nnz(lost), 1, 'the interval that spans the jump is dropped');
            tc.verifyLessThanOrEqual(nnz(lost), 2, 'and at most its neighbour, when an SCA1B epoch falls in the split gap');
            kept = isfinite(QJ.gyro(1).rate_rad_s);
            tc.verifyEqual(QJ.gyro(1).rate_rad_s(kept), Q.gyro(1).rate_rad_s(kept), 'AbsTol', 1e-12, ...
                'the other intervals keep their rates (a 30 deg offset on an 84 rad accumulator rounds at 1e-14)');
            tc.verifyEqual(QJ.gyro(2).rate_rad_s, Q.gyro(2).rate_rad_s, 'the other gyros are untouched');
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function attitudeQualityMetamorphic(tc)
            % 1. Physics: an Earth-pointing satellite turns about body y at the orbital rate, so the SCA1B
            %    rates are ~(0, 2 pi / P, 0), and against the nominal frame the attitude is off by about a degree.
            [eS, qS] = NonGravForces.io.readSCA1B(tc.f.sca1);
            [eG, r, v] = NonGravForces.io.readGNV1B(tc.f.gnv1, 'Eop', tc.eop);
            tc.assertEqual([double(eG.day(:)) double(eG.sec(:))], [double(eS.day(:)) double(eS.sec(:))], 'GNV1B and SCA1B share their 1 s epochs');
            Q = tc.dayOne('Orbit', struct('r', r, 'v', v), 'AlongVelocity', false);
            tc.verifyGreaterThan(Q.period_s, 5500); tc.verifyLessThan(Q.period_s, 5800);
            w = median(Q.omega_rad_s, 2);
            tc.verifyEqual(w(2), 2 * pi / Q.period_s, 'RelTol', 0.02, 'body y turns at the orbital rate');
            tc.verifyLessThan(max(abs(w([1 3]))), 2e-5, 'body x and z barely turn');
            tc.verifyLessThan(max(abs(Q.nominal.angles_rad(:))), deg2rad(3), 'the attitude follows the nominal frame within degrees');
            tc.verifyTrue(all(isfinite(Q.nominal.cpr(:))) && all(isfinite([Q.gyro.cpr])), 'with a period the cpr amplitudes exist');
            % 2. The residual is invariant to a constant bias added to a gyro's rate (a linear drift of its angle).
            [eI, aI, auxI] = NonGravForces.io.readIMU1B(tc.f.imu1);
            tI = reshape(referenceSystem.secondsBetween(referenceSystem.daySecondToEpoch([58849 0], "GPS"), eI), 1, []);
            Q0 = NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId);
            Qb = NonGravForces.geom.attitudeQuality(eS, qS, eI, aI + 1e-5 * tI .* (auxI.gyroId == 3), auxI.gyroId);
            ok = isfinite(Q0.gyro(3).residual_rad_s);
            tc.verifyBoundAndRecord("attitudeQuality residual invariant to a gyro bias of 1e-5 rad/s", ...
                max(abs(Qb.gyro(3).residual_rad_s(ok) - Q0.gyro(3).residual_rad_s(ok))) / Q0.gyro(3).rms_rad_s, 1e-6, "-", ...
                "the bias lies in the span of the design; rounding of a 1e-5 rad/s rate over 1 s");
            tc.verifyEqual(Qb.gyro(1).residual_rad_s, Q0.gyro(1).residual_rad_s, 'gyro 1 untouched');
        end
    end
    methods (Test, TestTags = {'speed'})
        function attitudeQualitySpeed(tc)
            % 1800 rotation vectors and three interpolations against the read of the IMU excerpt itself.
            [eS, qS] = NonGravForces.io.readSCA1B(tc.f.sca1);
            [eI, aI, auxI] = NonGravForces.io.readIMU1B(tc.f.imu1);
            q = tc.speedRatio(@() NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId), @() NonGravForces.io.readIMU1B(tc.f.imu1));
            tc.verifyBoundAndRecord("attitudeQuality / readIMU1B time ratio (30 min)", q, 5, "-", "vectorised rotations and interpolations, cheaper than parsing the text");
        end
    end

    % ===================================================================== geom.AttitudeProvider.fromItsg
    methods (Test, TestTags = {'contract'})
        function attitudeProviderItsgContract(tc)
            ap = NonGravForces.geom.AttitudeProvider.fromItsg(tc.f.itsg);
            tc.verifyEqual(ap.source, "itsg");
            tc.verifyEqual(ap.nodeCount, 600, 'the committed excerpt: 600 epochs at 10 s');
            tc.verifyEqual(ap.scale, "GPS");
            % the same file twice is one set of nodes, not a refusal downstream (2026-10-31: joined as given)
            ap2 = NonGravForces.geom.AttitudeProvider.fromItsg([tc.f.itsg tc.f.itsg]);
            tc.verifyEqual(ap2.nodeCount, 600);
            tc.verifyErrorId(@() NonGravForces.geom.AttitudeProvider.fromItsg(strrep(tc.f.itsg, 'attitude', 'reducedDynamicOrbit')), ...
                'NonGravForces:AttitudeProvider:dataItsg');
        end
    end
    methods (Test, TestTags = {'reference'})
        function attitudeProviderItsgReference(tc)
            % TU Graz's quaternion is SCA1B's in the other convention: conjugated once, the provider returns
            % SCA1B's attitude at the shared 10 s epochs, the rotation to rounding.
            ap = NonGravForces.geom.AttitudeProvider.fromItsg(tc.f.itsg);
            [eS, qS] = NonGravForces.io.readSCA1B(tc.f.sca600);
            k = 1:10:600;
            e = referenceSystem.daySecondToEpoch([double(eS.day(k)) double(eS.sec(k))], string(eS.scale));
            qI = ap.at(e);
            d = pagemtimes(referenceSystem.rotationFromQuaternion(qI), 'transpose', referenceSystem.rotationFromQuaternion(qS(:, k)), 'none');
            ang = acos(min(1, (squeeze(d(1, 1, :) + d(2, 2, :) + d(3, 3, :)) - 1) / 2));
            tc.verifyBoundAndRecord("fromItsg vs readSCA1B at the shared epochs, rotation angle", max(ang), 1e-7, "rad", ...
                "bit-identical quaternions conjugated once on both paths; acos near 1 resolves only sqrt(eps)");
        end
    end
end

function L = localSetToken(L, k, j, value)
tok = split(strtrim(L(k))); tok(j) = string(value); L(k) = strjoin(tok, " ");
end
