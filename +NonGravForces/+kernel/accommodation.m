function [acc, info] = accommodation(model, species, T_K, options)
%ACCOMMODATION Energy accommodation of the gas-surface interaction per epoch: fixed, either Pilinski fit, the Walker/Mehta adsorbate mixture, or SESAM.
%
%   acc = NonGravForces.kernel.accommodation("fixed", species, T_K, 'Alpha', 0.85)
%   acc = NonGravForces.kernel.accommodation("pilinski2010", species, T_K)
%   acc = NonGravForces.kernel.accommodation("pilinski2011", species, T_K)
%   acc = NonGravForces.kernel.accommodation("walker2014", species, T_K, 'SurfaceMass_amu', 263.3223)
%   acc = NonGravForces.kernel.accommodation("sesam", species, T_K, 'RelativeSpeed_ms', 7600)
%   [acc, info] = NonGravForces.kernel.accommodation(___)
%
%   Returns what kernel.dragPanelForce needs to form the re-emitted momentum:
%   either ONE accommodation coefficient per epoch, or a surface that is a
%   mixture of two states. The kernel never branches on the model name; it reads
%   acc.coverage, which is empty for a single state.
%
%   MODELS
%     fixed         alpha = Alpha at every epoch (Alpha scalar or 1xK). The value is
%                   the caller's statement; nothing is assumed when it is omitted.
%     pilinski2010  alpha = K n_O T / (1 + K n_O T) with the atomic-oxygen number
%     pilinski2011  density n_O and the gas temperature T. ONE MODEL, TWO PUBLISHED
%                   FITS, and each option is named after the source that fits it:
%                     pilinski2010  K = 7.50e-17 m^3 K^-1  Pilinski, Argrow & Palo
%                                   (2010) eq. 8, fitted to the Bowman & Moe sphere
%                                   accommodations, mean error -0.9 %, sd 1.7 %
%                     pilinski2011  K = 4.98e-17 m^3 K^-1  Pilinski (2011) sect. 7.3,
%                                   eqs 7.1-7.2, the author's later refit
%                   Both assume a CLEAN SUBSTRATE of alpha 0 and are fitted between
%                   150 and 500 km, where the fitted alpha spans 0.85-1.00; outside
%                   that range it is an extrapolation and info says so. Until
%                   2026-09-16 a single option pilinski2010 carried the 2011 value,
%                   so its name and its constant cited different papers; the owner
%                   ruled both should ship (HANDOVER LP-D, WP-D1).
%                   K IS TIED TO THE DENSITY MODEL IT WAS FITTED WITH. Pilinski
%                   (2010) is explicit that n_O and T came from NRLMSISE-00 and that
%                   "the application of this value should only be performed with
%                   NRLMSISE-00"; a different model needs K refitted. This toolbox
%                   does not refuse another backend -- the precise preset evaluates
%                   msis21, which runs ~10 % below NRLMSISE-00 at GRACE altitudes --
%                   but that combination is outside the fit and is stated here
%                   because nothing in the numbers would reveal it.
%     walker2014    The adsorbate model as Paul, Sheridan, Licata & Mehta (2022, eqs
%                   10-13) state it: a fraction theta of the surface carries adsorbed
%                   atomic oxygen and accommodates fully, the rest is clean substrate,
%                     theta   = K_L P_O / (1 + K_L P_O),  P_O = n_O k T  [Pa],  K_L = 1.44e6 Pa^-1
%                     alpha_a = 1
%                     alpha_s = 3 mu / (1 + mu)^2,  mu = (sum_j x_j m_j) / m_s
%                   with x_j the mole fractions of H, He, N, N2, O and O2 (normalised
%                   over those six) and m_s the mass of a surface lattice particle; the
%                   drag of the two states is mixed with weight theta. ATTRIBUTION, as
%                   that paper gives it: the coverage law and K_L to Walker, Mehta &
%                   Koller (2014, JSR 51(5) *1544*, the CLL drag-coefficient model), the
%                   substrate coefficient to Mehta et al. (2022); the name of this option
%                   follows the first. THE PAGE MATTERS: that issue carries two papers by
%                   the same three authors, and the other one (51(5) 1522, DOI
%                   10.2514/1.A32668, read 2026-09-16) does NOT contain this law -- it
%                   mentions the Langmuir isotherm once and credits it to Pilinski. 1544
%                   is not in hand, so the attribution is still Paul et al. 2022's and
%                   not this project's reading (HANDOVER O-17). K_L's unit is not
%                   printed; Pa^-1 follows from P_O formed in SI in the WVU RSM Suite
%                   that implements it (read as Torr^-1, theta would be 0.13 instead of
%                   0.95 at n_O = 1e15 m^-3, 1000 K). m_s has no universally
%                   agreed value: 263.3223 amu (Mehta et al. 2022, as cited by Paul
%                   2022), 65 amu (Pilinski 2011, Table 11.1); the owner chose the first
%                   for the precise preset (HANDOVER Part I, 2026-09-15), and it is
%                   required here, never defaulted.
%     sesam         SESAM (Pilinski, Argrow, Palo & Bowman 2013). The same Langmuir
%                   mixture as walker2014, with two differences that matter more than
%                   the shared algebra:
%                     1. K_L IS A FUNCTION OF INCIDENT KINETIC ENERGY, not a constant:
%                          K_L(E) = s_o(E, Eb, Ta) K_L,o + K_L,f          (eq. 22)
%                        where s_o is the initial sticking coefficient, the fraction of
%                        the incident Maxwellian whose ram energy falls below the binding
%                        energy Eb (eq. 26b, closed form). A faster satellite therefore
%                        adsorbs less in the same gas, which is the effect SESAM exists
%                        to model and the reason it needs the FLOW and not only the gas.
%                     2. P_O IS A RAM PRESSURE, NOT n_O k T (eq. 15, Appendix B):
%                          P_O = 1/2 rho_O V_r^2 C_D,in,O(s)
%                        with C_D,in,O the free-molecular sphere drag coefficient with
%                        the re-emission term set to zero. At 520 km this is 140x the
%                        n_O k T that walker2014's constant was fitted against, so the
%                        two pressures are NOT interchangeable and neither are their
%                        Langmuir constants (K_L,o = 5e6 torr^-1 = 3.75e4 Pa^-1 against
%                        Walker's 1.44e6 Pa^-1). Feeding this model the other pressure
%                        would under-predict its coverage by that factor, silently.
%                   Table C1 parameters, all fixed and none of them options: Eb = 5.7 eV,
%                   K_L,o = 5e6 torr^-1, K_L,f = 3e4 torr^-1, m_s = 65 amu. The substrate
%                   coefficient is eq. (19), alpha_s = 2.4 mu/(1+mu)^2 -- SESAM'S OWN 2.4,
%                   where walker2014 ships 3.0 from Paul et al. (2022); the two sources
%                   disagree and each option follows its own (HANDOVER O-42).
%                   NOT IMPLEMENTED FROM THAT PAPER: the DSMC lookup tables and the shape
%                   factor for tumbling spheres, which map alpha to a drag coefficient.
%                   This toolbox does not need them -- it has its own panel geometry and
%                   Sentman kernel, and takes only alpha from SESAM (HANDOVER O-16).
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   model     1x1   string  -      "fixed" | "pilinski2010" | "pilinski2011" | "walker2014" | "sesam"
%   species   8xK   double  1/m^3  Number densities in const.speciesName order (env.density)
%   T_K       1xK   double  K      Gas kinetic temperature at the satellite (density's T(2,:))
%
% NAME-VALUE ARGUMENTS
%   Name              Size      Type    Units  Description
%   ----              ----      ----    -----  -----------
%   Alpha             1x1|1xK   double  -      Coefficient in [0, 1]; required by "fixed", refused by the others
%   SurfaceMass_amu   1x1       double  amu    Surface lattice particle mass, > 0; required by "walker2014", refused by "sesam" (which fixes 65)
%   RelativeSpeed_ms  1x1|1xK   double  m/s    Satellite speed relative to the gas, > 0; required by "sesam", refused by the others
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   acc       1x1   struct  -      .model; .alpha 1xK (fixed, pilinski20xx; for walker2014 and sesam
%                                  the adsorbate state's, 1); .coverage 1xK (walker2014, sesam: theta;
%                                  [] for a single state); .alphaSubstrate 1xK (walker2014, sesam; []
%                                  otherwise). The EFFECTIVE coefficient of a two-state surface is
%                                  theta + (1 - theta) alphaSubstrate; the kernel mixes the two
%                                  states' forces rather than their coefficients, so this unit does
%                                  not form it.
%   info      1x1   struct  -      .K, .pressureO_Pa (walker2014, sesam), .meanMass_amu (walker2014,
%                                  sesam), .insideFit (1xK logical, EVERY model: whether each epoch
%                                  lies inside the range the model was fitted over -- alpha in
%                                  0.85-1.00 for pilinski20xx, P_O >= 1e-5 Pa for walker2014, true
%                                  everywhere for fixed, which extrapolates nothing, and for sesam,
%                                  whose source states no domain this project can transcribe),
%                                  .stickingCoefficient .langmuir_invPa .speedRatioO (sesam), .provenance
%
% ACCURACY
%   Exact evaluation of the published formulas, and the model error is the
%   gas-surface interaction itself. March et al. (2021) OPTIMISE A CONSTANT
%   alpha per mission rather than computing one -- 0.85 (CHAMP) and 0.82 (GOCE),
%   recommending 0.85 and superseding Doornbos (2011)'s 0.93 -- because the
%   in-situ composition a computed model needs is not measured on those
%   satellites. That is what all three presets use.
%   THE COMPUTED MODELS ARE STRONGLY ALTITUDE- AND SOLAR-CYCLE DEPENDENT, and
%   the spread between them is not small. Measured 2026-09-17 at the equator at
%   noon, F10.7 = 70 (low) and 180 (high), Ap 4 and 9:
%     300 km, low   pilinski2010 0.948   walker2014 0.857  (March 2021 reports
%                   SESAM oscillating between 0.85 and 0.93 at these altitudes)
%     520 km, low   pilinski2010 0.162   walker2014 0.167
%     520 km, high  pilinski2010 0.830   walker2014 0.634
%   So these models reproduce the published range where it was published, and
%   predict a factor-of-five swing over the solar cycle at GRACE-FO's altitude.
%   On the reference day (2020-01-01, deep minimum) walker2014 gives 17.5 % more
%   drag than the fitted 0.85. Choosing between them is a statement about the
%   surface, not a fidelity tier (HANDOVER O-42).
%   The substrate mass of walker2014 alone moves its effective coefficient by
%   ~0.1 at GRACE-FO altitude near solar minimum (measured with the three values
%   in circulation, 2026-09-15).
%
% WARNINGS
%   NonGravForces:accommodation:outsideFit  a computed model was evaluated outside the range it was
%       fitted in. Raised ONCE PER CALL, naming how many of the K epochs and which range; silence it
%       with warning('off', 'NonGravForces:accommodation:outsideFit') when the extrapolation is
%       intended. Owner ruling of 2026-09-17 (HANDOVER O-42): warn rather than refuse, because
%       extrapolating a fitted model knowingly is legitimate and refusing would forbid it. THIS AND
%       dependencyReport ARE THE ONLY UNITS THAT MAY WARN; a warning is not printing, since it
%       carries an identifier and can be silenced.
%       IT WILL FIRE ON GRACE-FO. At 520 km in solar minimum P_O is ~1.4e-8 Pa, three orders below
%       walker2014's fitted domain, and the adsorption models return 0.11-0.17 where the fitted
%       literature reports 0.82-0.93 (measured 2026-09-17). That is the extrapolation this warning
%       exists to name, and it is why all three presets take a fixed alpha instead.
%
% ERRORS
%   NonGravForces:accommodation:inputShape       species not 8xK or T_K not 1xK with the same K
%   NonGravForces:accommodation:inputAlpha       fixed without 'Alpha', Alpha outside [0, 1] or not 1x1/1xK, or Alpha given to another model
%   NonGravForces:accommodation:inputSurfaceMass walker2014 without a positive 'SurfaceMass_amu', or a mass given to sesam
%   NonGravForces:accommodation:inputSpeed       sesam without a positive 'RelativeSpeed_ms', or a speed given to another model
%   NonGravForces:accommodation:numericNaN       NaN/Inf in Alpha, species or T_K
%   NonGravForces:accommodation:inputRange       a negative density or a non-positive temperature (models that read them)
%
% REFERENCES
%   Pilinski M.D., Argrow B.M., Palo S.E. (2010) J. Spacecraft Rockets 47(6) 951-956,
%   DOI 10.2514/1.49330, eq. 8 and the K of the Results section;
%   Pilinski M.D. (2011) Dynamic gas-surface interaction modeling for satellite
%   aerodynamic computations, PhD thesis, University of Colorado, sect. 7.3 and
%   Table 11.1. Walker A., Mehta P., Koller J. (2014) J. Spacecraft Rockets 51(5)
%   1544. Paul S.N., Sheridan P.L., Licata R.J., Mehta P.M. (2022) arXiv:2210.08364v2,
%   eqs 10-13. March G., van den IJssel J., Siemes C., Visser P.N.A.M., Doornbos E.N.,
%   Pilinski M. (2021) Gas-surface interactions modelling influence on satellite
%   aerodynamics and thermosphere mass density, J. Space Weather Space Clim. 11, 54,
%   DOI 10.1051/swsc/2021035 -- the fitted alpha = 0.85 the presets use (corrected
%   2026-09-17; this block cited "9, A36", which resolves to no paper, HANDOVER O-42).
%   Pilinski M.D., Argrow B.M., Palo S.E., Bowman B.R. (2013) Semi-empirical satellite
%   accommodation model for spherical and randomly tumbling objects, J. Spacecraft
%   Rockets 50(3) 556-571, DOI 10.2514/1.A32348: eqs (15)-(22) and (26b), Table C1.
%   Eq. (15)'s second bracket is transcribed as the classical Schaaf-Chambre sphere
%   coefficient (4s^4 + 4s^2 - 1)/(2s^4), which is what Appendix B derives; the scan in
%   hand reads +1 there. At LEO s ~ 8 the two differ by 2.4e-04 relative.
%
% EXAMPLE
%   [~, sp, T] = NonGravForces.env.density("msis00", ep, lla, sw);
%   acc = NonGravForces.kernel.accommodation("walker2014", sp, T(2, :), 'SurfaceMass_amu', 263.3223);
%   [acc.coverage; acc.alphaSubstrate]
%
% See also NonGravForces.kernel.dragPanelForce, NonGravForces.kernel.gsiCoefficients, NonGravForces.env.density.
%
% KEYWORDS accommodation coefficient adsorption Langmuir atomic oxygen gas-surface interaction
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic); pilinski2010 and walker2014 2026-09-15,
% pilinski2011 and the 2010 constant from the paper itself 2026-09-16, by Claude Opus 5 (Anthropic).
arguments
    model (1,1) string {mustBeMember(model, ["fixed","pilinski2010","pilinski2011","walker2014","sesam"])}
    species double
    T_K double
    options.Alpha double = []
    options.SurfaceMass_amu double = []
    options.RelativeSpeed_ms double = []
end
unit = "accommodation";
K = size(T_K, 2);
if ~(ismatrix(species) && size(species, 1) == 8 && size(species, 2) == K && isrow(T_K) || (K == 0 && isempty(species)))
    NonGravForces.fail(unit, "inputShape", "species must be 8xK and T_K 1xK with the same K -> got species %s, T_K %s", mat2str(size(species)), mat2str(size(T_K)));
end
if any(~isfinite(species(:))) || any(~isfinite(T_K(:)))
    NonGravForces.fail(unit, "numericNaN", "species or T_K contain NaN/Inf -> screen the density output");
end
if model ~= "sesam" && ~isempty(options.RelativeSpeed_ms)
    NonGravForces.fail(unit, "inputSpeed", "model %s does not read a relative speed -> remove 'RelativeSpeed_ms', or use model sesam whose pressure is a ram pressure", model);
end
c = NonGravForces.const();
acc = struct('model', model, 'alpha', [], 'coverage', [], 'alphaSubstrate', []);
info.K = K;
switch model
    case "fixed"
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
        acc.alpha = a + zeros(1, K);
        info.insideFit = true(1, K);            % a stated value extrapolates nothing
    case {"pilinski2010", "pilinski2011"}
        % ONE MODEL, TWO PUBLISHED FITS, and each option is named after the source that fits it.
        % Until 2026-09-16 the single option pilinski2010 carried the 2011 thesis's constant, so its
        % name and its value cited different papers; the owner ruled both should ship (LP-D WP-D1).
        localNoAlpha(unit, model, options.Alpha);
        localRange(unit, species, T_K);
        if model == "pilinski2010"
            KL = 7.50e-17;                                 % m^3 K^-1, Pilinski, Argrow & Palo (2010) eq. 8
        else
            KL = 4.98e-17;                                 % m^3 K^-1, Pilinski (2011) sect. 7.3
        end
        x = KL * species(2, :) .* T_K;                      % n_O T in m^-3 K
        acc.alpha = x ./ (1 + x);
        info.insideFit = acc.alpha >= 0.85 & acc.alpha <= 1;
        localWarnOutsideFit(model, info.insideFit, "the fitted alpha of 0.85 to 1.00, reached between 150 and 500 km");
    case "walker2014"
        localNoAlpha(unit, model, options.Alpha);
        ms = options.SurfaceMass_amu;
        if ~(isscalar(ms) && isfinite(ms) && ms > 0)
            NonGravForces.fail(unit, "inputSurfaceMass", "walker2014 needs a positive 'SurfaceMass_amu' -> the sources disagree (263.3223 amu Mehta et al. 2022, 65 amu Pilinski 2011); state which one you use");
        end
        localRange(unit, species, T_K);
        KL = 1.44e6;                                        % Pa^-1, Walker et al. 2014 via Paul et al. 2022 eq. 11
        PO = species(2, :) .* c.kB .* T_K;                  % Pa
        acc.coverage = KL * PO ./ (1 + KL * PO);
        acc.alpha = ones(1, K);                             % the adsorbate state accommodates fully
        six = [6 1 7 3 2 4];                                % H He N N2 O O2, the species of Paul et al. 2022 eq. 13
        m = c.speciesMass(six)' / c.NA;                     % kg per molecule, 6x1
        n6 = species(six, :);
        meanMass = sum(n6 .* m, 1) ./ max(sum(n6, 1), realmin);   % mole-fraction weighted mean molecular mass
        mu = meanMass / (ms * c.amu);
        acc.alphaSubstrate = 3 * mu ./ (1 + mu).^2;
        info.pressureO_Pa = PO;
        info.meanMass_amu = meanMass / c.amu;
        % [W14a] fig. 17 fits theta well for P_O above ~1e-5 Pa (below ~300 km) and OVERESTIMATES
        % it between 1e-6 and 1e-5; below that the isotherm is extrapolated. The domain is a
        % pressure, not an alpha, so this model's insideFit is a statement about the environment it
        % was handed and not about the number it returned.
        info.insideFit = PO >= 1e-5;
        localWarnOutsideFit(model, info.insideFit, "the atomic-oxygen partial pressure of 1e-5 Pa and above, which is below about 300 km");
    case "sesam"
        % Pilinski et al. (2013), the model O-16 has named since LP-5. It is the Langmuir isotherm
        % again, but with the adsorption constant itself a function of the INCIDENT KINETIC ENERGY,
        % so a fast satellite adsorbs less than a slow one in the same gas.
        localNoAlpha(unit, model, options.Alpha);
        if ~isempty(options.SurfaceMass_amu)
            NonGravForces.fail(unit, "inputSurfaceMass", "sesam fixes its substrate mass at 65 amu (Pilinski et al. 2013 Table C1) -> remove 'SurfaceMass_amu', or use model walker2014 which takes one");
        end
        vr = options.RelativeSpeed_ms;
        if isempty(vr) || ~(isscalar(vr) || isequal(size(vr), [1 K]))
            NonGravForces.fail(unit, "inputSpeed", "sesam needs 'RelativeSpeed_ms' as a scalar or 1xK -> its pressure is a RAM pressure, so it needs the flow and not only the gas state");
        end
        if any(~isfinite(vr(:))) || any(vr(:) <= 0)
            NonGravForces.fail(unit, "inputSpeed", "RelativeSpeed_ms must be finite and > 0 -> got %g", vr(find(~isfinite(vr) | vr <= 0, 1)));
        end
        localRange(unit, species, T_K);
        vr = vr + zeros(1, K);
        Eb = 5.7 * c.eV;                                    % J, adsorption energy, Table C1
        KLo = 5e6 / c.torr;                                 % Pa^-1, initial Langmuir parameter, Table C1
        KLf = 3e4 / c.torr;                                 % Pa^-1, final Langmuir parameter, Table C1
        ms = 65;                                            % amu, average substrate mass, Table C1
        mO = c.speciesMass(2) / c.NA;                       % kg per atom of atomic oxygen
        kT = c.kB * T_K;                                    % J
        Er = 0.5 * mO * vr.^2;                              % J, ram kinetic energy, eq. (15) text
        % --- eq. (26b), the initial sticking coefficient, in the only form that survives LEO.
        % As printed it evaluates exp(Eb/kT) ~ 1e36 and exp(2 sqrt(Eb Er)/kT) ~ 1e65 against an
        % exp(-(Eb+Er)/kT) ~ 1e-66. Combining the exponents first gives two terms in [0, 1]:
        %   exp(-(Eb+Er)/kT) exp(Eb/kT)            = exp(-Er/kT)
        %   exp(-(Eb+Er)/kT) exp(2 sqrt(Eb Er)/kT) = exp(-(sqrt(Eb) - sqrt(Er))^2/kT)
        % and the printed erf(a) + erf(b) cancels to nothing once Er > Eb, so it is written with the
        % identity erf(a) + erf(b) = erfc(-a) - erfc(b), which carries no +-1 at all.
        a = (sqrt(Eb) - sqrt(Er)) ./ sqrt(kT);
        b = sqrt(Er ./ kT);
        g = sqrt(pi * kT .* Er);
        num = g .* (erfc(-a) - erfc(b)) + kT .* (exp(-Er ./ kT) - exp(-(sqrt(Eb) - sqrt(Er)).^2 ./ kT));
        den = g .* (1 + erf(b)) + kT .* exp(-Er ./ kT);
        s0 = min(max(num ./ den, 0), 1);                    % a normalised fraction by construction
        KL = s0 * KLo + KLf;                                % eq. (22), Pa^-1
        % --- eq. (15), and this is the fact that makes sesam a different model and not a constant
        % swap: its pressure is a RAM pressure, 140x the n_O kB T that walker2014's constant was
        % fitted against (HANDOVER Part G, LP-E, WP-E3 scope).
        sO = vr ./ sqrt(2 * c.kB * T_K / mO);               % atomic-oxygen speed ratio
        cdIn = (2 * sO.^2 + 1) ./ (sqrt(pi) * sO.^3) .* exp(-sO.^2) ...
             + (4 * sO.^4 + 4 * sO.^2 - 1) ./ (2 * sO.^4) .* erf(sO);   % sphere, re-emission removed
        rhoO = species(2, :) * mO;                          % kg m^-3
        PO = 0.5 * rhoO .* vr.^2 .* cdIn;                   % Pa
        acc.coverage = KL .* PO ./ (1 + KL .* PO);          % eq. (17)
        acc.alpha = ones(1, K);                             % eq. (18): the adsorbate state accommodates fully
        six = [6 1 7 3 2 4];                                % H He N N2 O O2, eq. (16)
        m6 = c.speciesMass(six)' / c.NA;
        n6 = species(six, :);
        meanMass = sum(n6 .* m6, 1) ./ max(sum(n6, 1), realmin);
        mu = meanMass / (ms * c.amu);                       % eq. (20)
        acc.alphaSubstrate = 2.4 * mu ./ (1 + mu).^2;       % eq. (19) -- SESAM's own 2.4, not Paul's 3.0
        info.pressureO_Pa = PO;
        info.meanMass_amu = meanMass / c.amu;
        info.stickingCoefficient = s0;
        info.langmuir_invPa = KL;
        info.speedRatioO = sO;
        % SESAM states no domain of validity this project can transcribe: [P13] fits its parameters
        % over a set of objects and orbits rather than over a stated range of a variable, so there
        % is nothing to compare an epoch against. insideFit is therefore TRUE everywhere and says
        % so in the help, rather than being quietly absent or quietly invented.
        info.insideFit = true(1, K);
end
if nargout > 1
    info.provenance = NonGravForces.provenance("kernel.accommodation", 'Options', struct('model', model, 'Alpha', options.Alpha, 'SurfaceMass_amu', options.SurfaceMass_amu, 'RelativeSpeed_ms', options.RelativeSpeed_ms));
end
end

function localWarnOutsideFit(model, inside, domain)
%LOCALWARNOUTSIDEFIT Warn once per call when a computed accommodation leaves the range it was fitted in.
%   OWNER RULING 2026-09-17 (HANDOVER O-42): warn, rather than refuse or report only. The flag has
%   been computed since 2026-09-07 and NO CALLER HAS EVER READ IT, so a user who selected an
%   adsorption model at GRACE-FO altitude received a coefficient five times below every published
%   fit with nothing to say so -- a silent degradation under R4, with the instrument to end it
%   already present.
%
%   THIS IS THE SECOND UNIT IN THE TOOLBOX ALLOWED TO WARN, after dependencyReport, and it is a
%   deliberate exception to "no printing from library code" (CLAUDE.md sect. 3): a warning is not
%   printing -- it carries an identifier, the caller can silence it with warning('off', id), and it
%   cannot be mistaken for a result. Refusing was the alternative and was rejected: extrapolating a
%   fitted model is a legitimate thing to do knowingly, which is exactly what R5 asks for.
%
%   It fires ONCE PER CALL and not once per epoch: this unit is vectorised over K, so an arc of
%   86400 seconds outside the fit is one warning naming how many epochs, not 86400 of them.
if all(inside)
    return
end
warning("NonGravForces:accommodation:outsideFit", ...
    "accommodation %s: %d of %d epochs lie outside %s. The returned coefficient is an extrapolation of a fitted model -> read info.insideFit, choose a fixed alpha (March et al. 2021 recommend 0.85), or silence this with warning('off', 'NonGravForces:accommodation:outsideFit').", ...
    model, nnz(~inside), numel(inside), domain);
end

function localNoAlpha(unit, model, a)
if ~isempty(a)
    NonGravForces.fail(unit, "inputAlpha", "model %s computes its coefficient from the environment -> remove 'Alpha', or use model fixed", model);
end
end

function localRange(unit, species, T_K)
if any(species(:) < 0) || any(T_K <= 0)
    NonGravForces.fail(unit, "inputRange", "number densities must be >= 0 and T_K > 0 -> check the density output");
end
end
