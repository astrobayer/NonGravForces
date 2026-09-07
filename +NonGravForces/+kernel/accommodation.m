function [alpha, info] = accommodation(model, species, T_K, options)
%ACCOMMODATION Energy accommodation coefficient per epoch: fixed (sesam and langmuir refuse until their constants are in hand).
%
%   alpha = NonGravForces.kernel.accommodation("fixed", species, T_K, 'Alpha', 0.85)
%   [alpha, info] = NonGravForces.kernel.accommodation(___)
%
%   Returns the coefficient the gas-surface interaction kernel uses per epoch,
%   expanded to 1xK so that the caller never branches on the model.
%
%   MODELS
%     fixed     alpha = Alpha at every epoch (Alpha scalar or 1xK). The value is
%               the caller's statement; nothing is assumed when it is omitted
%               (rule R5: the option is empty by default and its omission errors).
%     sesam     Pilinski, Argrow & Palo (2010; 2013): a Langmuir isotherm in the
%               atomic-oxygen partial pressure p_O = n_O k T interpolating
%               between a clean and an adsorbate-covered surface. REFUSED
%               (notImplemented) until the constants are read from the paper
%               (HANDOVER O-16); none is guessed.
%     langmuir  Same isotherm with caller-supplied constants: refused with sesam.
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   model     1x1   string  -      "fixed" | "sesam" | "langmuir"
%   species   8xK   double  1/m^3  Number densities in const.speciesName order (from env.density);
%                                  read by sesam/langmuir, validated for shape only by fixed
%   T_K       1xK   double  K      Gas kinetic temperature at the satellite (density's T(2,:))
%
% NAME-VALUE ARGUMENTS
%   Name      Size      Type    Units  Description
%   ----      ----      ----    -----  -----------
%   Alpha     1x1|1xK   double  -      Coefficient in [0, 1] for "fixed"; required by that model
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   alpha     1xK   double  -      Accommodation coefficient per epoch
%   info      1x1   struct  -      .model, .K, .provenance
%
% ACCURACY
%   fixed: exact pass-through. The published spread of the coefficient itself
%   is the model error: 0.82-0.93 at 270-500 km (March et al. 2019; Pilinski &
%   Argrow 2013), i.e. 5-10 % in the drag coefficient.
%
% ERRORS
%   NonGravForces:accommodation:notImplemented   sesam or langmuir (O-16)
%   NonGravForces:accommodation:inputShape       species not 8xK or T_K not 1xK with the same K
%   NonGravForces:accommodation:inputAlpha       fixed without 'Alpha', or Alpha not in [0, 1] / not 1x1 or 1xK
%   NonGravForces:accommodation:numericNaN       NaN/Inf in Alpha, species or T_K
%
% REFERENCES
%   Pilinski M.D., Argrow B.M., Palo S.E. (2010) Semiempirical model for satellite
%   energy-accommodation coefficients, J. Spacecraft Rockets 47(6). Pilinski et al.
%   (2013) JSR 50(3). March G. et al. (2019) J. Space Weather Space Clim. 9, A36.
%
% EXAMPLE
%   [~, sp, T] = NonGravForces.env.density("msis00", ep, lla, sw);
%   alpha = NonGravForces.kernel.accommodation("fixed", sp, T(2, :), 'Alpha', 0.85);
%
% See also NonGravForces.kernel.gsiCoefficients, NonGravForces.kernel.dragPanelForce, NonGravForces.env.density.
%
% KEYWORDS accommodation coefficient SESAM Langmuir gas-surface interaction
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    model (1,1) string {mustBeMember(model, ["fixed","sesam","langmuir"])}
    species double
    T_K double
    options.Alpha double = []
end
unit = "accommodation";
if any(model == ["sesam","langmuir"])
    NonGravForces.fail(unit, "notImplemented", ...
        "model %s needs the Langmuir-isotherm constants of Pilinski et al. (2010, 2013), which are not in hand (HANDOVER O-16) -> use fixed with an explicit 'Alpha', or supply the paper", model);
end
K = size(T_K, 2);
if ~(ismatrix(species) && size(species, 1) == 8 && size(species, 2) == K && isrow(T_K) || (K == 0 && isempty(species)))
    NonGravForces.fail(unit, "inputShape", "species must be 8xK and T_K 1xK with the same K -> got species %s, T_K %s", mat2str(size(species)), mat2str(size(T_K)));
end
if any(~isfinite(species(:))) || any(~isfinite(T_K(:)))
    NonGravForces.fail(unit, "numericNaN", "species or T_K contain NaN/Inf -> screen the density output");
end
a = options.Alpha;
if isempty(a) || ~(isscalar(a) || isequal(size(a), [1 K]))
    NonGravForces.fail(unit, "inputAlpha", "fixed needs 'Alpha' as a scalar or 1xK -> pass 'Alpha', 0.85 (the value is your statement, not a default)");
end
if any(~isfinite(a(:)))
    NonGravForces.fail(unit, "numericNaN", "Alpha contains NaN/Inf -> check the input");
end
if any(a(:) < 0 | a(:) > 1)
    NonGravForces.fail(unit, "inputAlpha", "Alpha must lie in [0, 1] -> got %g", a(find(a < 0 | a > 1, 1)));
end
alpha = a + zeros(1, K);
info.model = model;
info.K = K;
info.provenance = NonGravForces.provenance("kernel.accommodation", 'Options', struct('model', model, 'Alpha', options.Alpha));
end
