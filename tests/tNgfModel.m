classdef tNgfModel < NgfTestCase
%TNGFMODEL Tests of NonGravForces.Model construction and its stubbed methods.
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    methods (Test, TestTags = {'contract'})
        function modelContract(tc)
            mm = NonGravForces.geom.graceFO("C");
            M = NonGravForces.Model(mm, NonGravForces.Options("fast"));
            tc.verifyEqual(M.options.preset, "fast");
            tc.verifyEqual(M.provenance.unit, "Model");
            tc.verifyErrorId(@() NonGravForces.Model(struct()), 'NonGravForces:Model:inputMacro');
            % Raised by Model itself since 2026-09-07, so the identifier its ERRORS block documents
            % is the one a caller actually sees (it used to be MATLAB:validation:UnableToConvert).
            tc.verifyErrorId(@() NonGravForces.Model(mm, struct()), 'NonGravForces:Model:inputOptions');
            for m = ["accel","torque","series","partials","rhs","budget"]
                tc.verifyErrorId(@() M.(m)(), 'NonGravForces:Model:notImplemented');
            end
        end
    end
    methods (Test, TestTags = {'robustness'})
        function modelRobustness(tc)
            % Default options and absent attitude are legal at construction.
            M = NonGravForces.Model(NonGravForces.geom.graceFO("D"));
            tc.verifyEqual(M.options.preset, "standard");
            tc.verifyEmpty(M.attitude);
        end
    end
end
