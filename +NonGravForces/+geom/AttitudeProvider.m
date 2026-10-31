classdef AttitudeProvider
%ATTITUDEPROVIDER BODY->GCRS quaternions at query epochs: SCA1B nodes interpolated by referenceSystem, or a nominal attitude from the state.
%
%   ap = NonGravForces.geom.AttitudeProvider(epochNodes, qNodes)                 % from io.readSCA1B
%   ap = NonGravForces.geom.AttitudeProvider(epochNodes, qNodes, 'Method', "squad")
%   q  = ap.at(epoch)                                                            % 4xK BODY->GCRS
%   q  = NonGravForces.geom.AttitudeProvider.nominal(r, v, xAlongVelocity)       % 4xK from the state
%   ap = NonGravForces.geom.AttitudeProvider.fromItsg(files)                     % TU Graz's attitude product
%
%   The attitude source of NonGravForces.Model (HANDOVER D.4: attitude is
%   stored BODY->GCRS, scalar first, and referenceSystem operates on it).
%   Interpolation is referenceSystem.interpolateQuaternions (slerp or squad),
%   which makes the node series sign-continuous first and refuses to
%   extrapolate: a query outside the node range raises its
%   queryOutOfRange error unwrapped, because nothing in the series
%   constrains the attitude there (HANDOVER E.2 rule 5).
%
%   NOMINAL ATTITUDE. nominal(r, v, xAlongVelocity) builds the ideal
%   Earth-pointing frame from the GCRS state: body +Z to nadir (-r/|r|), body
%   +X along the in-plane horizontal direction, with the sign the caller
%   states (GRACE-C flew with +X against the velocity on 2020-01-01, measured
%   by tests/tNgfIo.m; GRACE-D the other way), body +Y completing the
%   right-handed triad (= -(orbit normal) for +X along the velocity). The
%   quaternion comes from referenceSystem.quaternionFromRotation of the
%   matrix whose columns are the body axes in the GCRS, the same convention
%   io.readSCA1B delivers (checked by the round trip in the tests).
%
%   TU GRAZ'S ATTITUDE (N-12, 2026-10-30). fromItsg(files) reads the ITSG
%   operational attitude (GROOPS instrument files, GPS time as MJD, 10 s) through
%   gracefoReading.readGroopsInstrument and conjugates it EXACTLY ONCE: TU Graz's
%   "SRF -> CRF" quaternion is this toolbox's "ICRF -> SRF" in the other
%   quaternion convention, and for GRACE-FO it is BIT-IDENTICAL to SCA1B's on the
%   same epochs (PR #108; CLAUDE.md section 7). For GRACE-FO it is therefore not a
%   second attitude source, only a second carrier of the first; for classic GRACE
%   (GRACE-1/-2) it is TU Graz's own product. The independent sensor of the
%   GRACE-FO attitude change is the IMU (io.readIMU1B, geom.attitudeQuality).
%
% CONSTRUCTOR INPUTS
%   Name         Size  Type    Units  Description
%   ----         ----  ----    -----  -----------
%   epochNodes   1x1   struct  -      referenceSystem epoch struct, M entries, strictly increasing
%   qNodes       4xM   double  -      BODY->GCRS quaternions at the nodes (io.readSCA1B output)
%
% NAME-VALUE
%   Method       1x1   string  -      "slerp" (default) | "squad"
%
% PROPERTIES (read-only)
%   source       string   "nodes", or "itsg" from fromItsg
%   method       string
%   epochNodes   struct
%   qNodes       4xM double, as given (sign continuity is the interpolator's)
%   nodeCount    double
%   scale        string   time scale of the nodes; queries are converted to it
%   span         1x2 double decimal day numbers of the first and last node, in .scale
%   convention   struct   referenceSystem.attitudeConvention()
%   provenance   struct
%
% METHODS
%   q = at(epoch)              4xK unit quaternions BODY->GCRS at the query epochs
%   [q, info] = at(epoch)      info from interpolateQuaternions plus .converted (logical)
%   q = nominal(r, v, xAlongVelocity)   static, 4xK from a 3xK state and a logical
%   ap = fromItsg(files, 'Method', m)   static, a provider over TU Graz's attitude files (one satellite, in
%                                       any order: sorted, an epoch held twice with the same rotation kept
%                                       once), conjugated once; source "itsg"
%
% ACCURACY
%   Interpolation: at the nodes the result reproduces the input to 1e-16;
%   between 1 Hz SCA1B nodes the slerp error against a smoothly turning
%   attitude is that of the sampling (referenceSystem measured 1.2e-2 rad for
%   1 s nodes over a fast synthetic turn; a GRACE-FO attitude turning at the
%   orbital rate 1.1e-3 rad/s is linear to 1e-7 rad over 1 s). squad is C^1
%   in the rate. Nominal: exact geometry; the real attitude differs from it by
%   the pointing offsets (~1 deg on the sample, tests/tNgfGeom.m).
%
% ERRORS
%   NonGravForces:AttitudeProvider:inputNodes      qNodes not 4xM, epochNodes not an epoch struct with M entries, M < 2
%   NonGravForces:AttitudeProvider:inputEpoch      at(): not an epoch struct
%   NonGravForces:AttitudeProvider:inputShape      nominal(): r or v not 3xK, or degenerate (r x v = 0)
%   NonGravForces:AttitudeProvider:dataItsg        fromItsg(): a file is not a 4-column quaternion product, or
%                                                  two files hold one epoch with different attitudes
%   gracefoReading:*                               fromItsg(): propagated unwrapped
%   referenceSystem:interpolateQuaternions:*       propagated unwrapped (queryOutOfRange in particular)
%   referenceSystem:*                              propagated unwrapped
%
% REFERENCES
%   HANDOVER D.4; referenceSystem.interpolateQuaternions, .quaternionFromRotation, .attitudeConvention.
%
% EXAMPLE
%   [ep, q] = NonGravForces.io.readSCA1B(file);
%   ap = NonGravForces.geom.AttitudeProvider(ep, q);
%   qk = ap.at(referenceSystem.daySecondToEpoch([58849 100.5], "GPS"));
%
% See also NonGravForces.io.readSCA1B, NonGravForces.geom.checkAttitudeSense, referenceSystem.interpolateQuaternions.
%
% KEYWORDS attitude provider quaternion interpolation slerp squad nominal
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
    properties (SetAccess = private)
        source (1,1) string = "nodes"
        method (1,1) string = "slerp"
        epochNodes struct
        qNodes double
        nodeCount (1,1) double = 0
        scale (1,1) string = ""
        span (1,2) double = [NaN NaN]
        convention struct
        provenance struct
    end
    methods
        function obj = AttitudeProvider(epochNodes, qNodes, options)
            arguments
                epochNodes (1,1) struct
                qNodes double
                options.Method (1,1) string {mustBeMember(options.Method, ["slerp","squad"])} = "slerp"
            end
            unit = "AttitudeProvider";
            if ~all(isfield(epochNodes, ["day","sec","scale"]))
                NonGravForces.fail(unit, "inputNodes", "epochNodes must be a referenceSystem epoch struct (.day .sec .scale) -> io.readSCA1B returns one");
            end
            M = numel(epochNodes.day);
            if ~(ismatrix(qNodes) && isequal(size(qNodes), [4 M]) && all(isfinite(qNodes(:)))) || M < 2
                NonGravForces.fail(unit, "inputNodes", "qNodes must be a finite 4x%d array for the %d node epochs (at least 2) -> io.readSCA1B output, one column per node", M, M);
            end
            obj.method = options.Method;
            obj.epochNodes = epochNodes;
            obj.qNodes = qNodes;
            obj.nodeCount = M;
            obj.scale = string(epochNodes.scale);
            t = double(epochNodes.day(:)) + double(epochNodes.sec(:)) / 86400;
            obj.span = [t(1) t(end)];
            obj.convention = referenceSystem.attitudeConvention();
            obj.provenance = NonGravForces.provenance("geom.AttitudeProvider", 'Options', struct('method', options.Method, 'nodes', M, 'span', obj.span));
        end
        function [q, info] = at(obj, epoch)
            %AT BODY->GCRS quaternions at the query epochs (converted to the node scale first).
            if ~(isstruct(epoch) && all(isfield(epoch, ["day","sec","scale"])))
                NonGravForces.fail("AttitudeProvider", "inputEpoch", "at() takes a referenceSystem epoch struct -> build one with daySecondToEpoch or mjdToEpoch");
            end
            converted = string(epoch.scale) ~= obj.scale;
            if converted
                epoch = referenceSystem.convertTimeScale(epoch, obj.scale);
            end
            [q, ~, info] = referenceSystem.interpolateQuaternions(obj.qNodes, obj.epochNodes, epoch, 'Method', obj.method);
            info.converted = converted;
        end
    end
    methods (Static)
        function obj = fromItsg(files, options)
            %FROMITSG A provider over TU Graz's operational attitude files, the quaternion conjugated exactly once.
            arguments
                files (1,:) string {mustBeNonempty}
                options.Method (1,1) string {mustBeMember(options.Method, ["slerp","squad"])} = "slerp"
            end
            t = zeros(0, 2); qf = zeros(0, 4);
            for k = 1:numel(files)
                [tk, qk] = gracefoReading.readGroopsInstrument(files(k));
                if size(qk, 2) ~= 4
                    NonGravForces.fail("AttitudeProvider", "dataItsg", "%s carries %d data columns, not the 4 of a quaternion -> pass TU Graz's attitude product", files(k), size(qk, 2));
                end
                t = [t; tk]; qf = [qf; qk]; %#ok<AGROW>
            end
            % Files given out of order, or overlapping at a day boundary, are put in order and an epoch held twice is
            % kept once -- when both copies are the same rotation (q or -q). Until 2026-10-31 the files were joined as
            % given and the interpolator downstream refused or, worse, interpolated across the overlap.
            [t, ord] = sortrows(t); qf = qf(ord, :);
            dup = [false; all(diff(t, 1, 1) == 0, 2)];
            if any(dup)
                k = find(dup);
                same = abs(sum(qf(k, :) .* qf(k - 1, :), 2)) >= 1 - 1e-12;     % |q1 . q2| = 1: the same rotation
                if ~all(same)
                    NonGravForces.fail("AttitudeProvider", "dataItsg", "the files hold epoch (day %d, %.3f s) twice with different attitudes -> pass files that do not overlap, or one version of each day", ...
                        t(k(find(~same, 1)), 1), t(k(find(~same, 1)), 2));
                end
                t = t(~dup, :); qf = qf(~dup, :);
            end
            ep = referenceSystem.daySecondToEpoch(t, "GPS");
            obj = NonGravForces.geom.AttitudeProvider(ep, referenceSystem.invertQuaternion(qf'), 'Method', options.Method);
            obj.source = "itsg";
            obj.provenance = NonGravForces.provenance("geom.AttitudeProvider.fromItsg", 'Files', files, 'Options', options);
        end
        function q = nominal(r, v, xAlongVelocity)
            %NOMINAL Ideal Earth-pointing attitude: +Z nadir, +X along (true) or against (false) the velocity.
            arguments
                r double
                v double
                xAlongVelocity (1,1) logical
            end
            unit = "AttitudeProvider";
            K = size(r, 2);
            if ~(ismatrix(r) && size(r, 1) == 3 && isequal(size(v), [3 K]) && all(isfinite([r(:); v(:)])))
                NonGravForces.fail(unit, "inputShape", "r and v must be finite 3xK arrays of equal K -> one GCRS state per column");
            end
            zb = -r ./ vecnorm(r);                                   % nadir
            h = cross(r, v, 1);
            if any(vecnorm(h) < 1e-6 * vecnorm(r) .* vecnorm(v))
                NonGravForces.fail(unit, "inputShape", "r x v vanishes at column %d -> a radial velocity has no orbit frame", find(vecnorm(h) < 1e-6 * vecnorm(r) .* vecnorm(v), 1));
            end
            along = cross(h ./ vecnorm(h), r ./ vecnorm(r), 1);     % in-plane horizontal, sense of motion
            if xAlongVelocity, xb = along; else, xb = -along; end
            yb = cross(zb, xb, 1);
            R = zeros(3, 3, K);
            R(:, 1, :) = reshape(xb, 3, 1, K); R(:, 2, :) = reshape(yb, 3, 1, K); R(:, 3, :) = reshape(zb, 3, 1, K);   % columns: body axes in GCRS
            q = referenceSystem.quaternionFromRotation(R);
        end
    end
end
