classdef Options
%OPTIONS Validated option set for NonGravForces.Model, with fidelity presets.
%
%   opt = NonGravForces.Options                       % preset "standard"
%   opt = NonGravForces.Options("precise")
%   opt = NonGravForces.Options("fast", "drag.density", "msis00")
%   opt = opt.set("srp.tsi", 1367)                    % legacy TSI, explicit
%
%   Configure once, evaluate many (principle P5). Every field is validated
%   at construction; an unknown field name errors, and EVERY ADMISSIBLE VALUE
%   EVALUATES: a value this toolbox cannot compute is not offered (rule R4,
%   HANDOVER LP-C finding F5). Presets fix the DEFAULT tier of every model
%   family; extended models (HWM14 winds, transient thermal, body torques) are
%   opt-ins no preset enables, so preset results never move when one lands.
%
%   FIDELITY LADDER
%     field                   fast          standard        precise
%     drag.density            exponential   msis00          msis21
%     drag.wind               corotation    corotation      corotation (+hwm14 opt-in)
%     drag.gsi                constantCd    sentman         sentman
%     drag.accommodation      fixed 0.85    fixed 0.85      fixed 0.85 (March et al. 2021)
%     srp.eclipse             conical       conical         solaarsCF
%     erp.backend             none          knocke          ceres
%     erp.grid                -             12x24           24x48
%     thermal.model           none          steadyState     steadyState (+transient opt-in)
%     geom.selfShadow         none          lut 16          lut 32
%     sun.source              rs.analytic   rs.analytic     rs.de440
%
% PROPERTIES (read-only after construction; change through set())
%   preset      string    "fast" | "standard" | "precise"
%   drag        struct    .density  "exponential" | "msis00" | "msis21"
%                         .wind     "none" | "corotation" | "hwm14"
%                         .gsi      "constantCd" | "sentman" | "maxwell"
%                         .accommodation "fixed" | "pilinski2010" | "pilinski2011" | "walker2014"
%                         .alpha    accommodation coefficient, read by "fixed" only ([] otherwise)
%                         .surfaceMass_amu  substrate particle mass, read by "walker2014" only
%                         .specularFraction specular fraction, read by "maxwell" only
%                         .Tw_K     wall temperature [K]
%                         .cdScale  estimated scale on the whole drag force
%                         .cd       the coefficient constantCd uses (2.2 for a convex LEO body;
%                                   Montenbruck & Gill 2005)
%   srp         struct    .tsi [W/m^2] .eclipse "cylindrical"|"conical"|"solaarsCF"
%                         .lunarEclipse (logical) .crScale (estimated scale on the whole
%                         radiation force; the optics are per panel, there is no lumped Cr)
%   erp         struct    .backend "none"|"pointSource"|"knocke"|"ceres" .grid [nRings nSectors].
%                         Earth radiance is Lambertian in every backend (no angular
%                         distribution models: the field that pretended otherwise is gone)
%   thermal     struct    .model "none"|"steadyState"|"transient" .internalPower_W
%   torques     struct    .gravityGradient .magnetic (logical)
%   geom        struct    .selfShadow "none"|"lut"|"raycast" .lutResolution (cube-map cells per face edge)
%   sun         struct    .source "rs.analytic"|"rs.de440"
%   series      struct    .useParallel (logical) .chunk (double, epochs)
%   partials    struct    .method "analytic" | "finiteDiff"
%
%   A field a model does not read is NOT silently ignored when it contradicts the
%   choice: Model refuses a non-empty drag.alpha beside a non-fixed accommodation,
%   a maxwell gas without drag.specularFraction, walker2014 without
%   drag.surfaceMass_amu, and internal power on a macro model with no shares.
%
% ERRORS
%   NonGravForces:Options:inputPreset   unknown preset
%   NonGravForces:Options:inputField    unknown dotted field name, or name-value overrides not in pairs
%   NonGravForces:Options:inputValue    value outside its admissible set
%
% EXAMPLE
%   opt = NonGravForces.Options("standard", "drag.wind", "hwm14");
%   opt.drag.wind
%
% See also NonGravForces.Model.
%
% KEYWORDS options presets fidelity validation
% Created 2026-08-27 by Claude Fable 5 (Anthropic); admissible sets reduced to what evaluates 2026-09-15 by Claude Opus 5 (Anthropic).
    properties (SetAccess = private)
        preset   (1,1) string = "standard"
        drag     struct
        srp      struct
        erp      struct
        thermal  struct
        torques  struct
        geom     struct
        sun      struct
        series   struct
        partials struct
    end
    methods
        function obj = Options(preset, varargin)
            arguments
                preset (1,1) string = "standard"
            end
            arguments (Repeating)
                varargin
            end
            if ~any(preset == ["fast","standard","precise"])
                NonGravForces.fail("Options", "inputPreset", ...
                    "unknown preset '%s' -> use fast, standard or precise", preset);
            end
            obj.preset = preset;
            obj = obj.applyPreset(preset);
            if mod(numel(varargin), 2) ~= 0
                NonGravForces.fail("Options", "inputField", ...
                    "name-value overrides must come in pairs -> check the argument list");
            end
            for k = 1:2:numel(varargin)
                obj = obj.set(varargin{k}, varargin{k+1});
            end
        end

        function obj = set(obj, field, value)
            %SET Override one dotted field, e.g. set("drag.gsi","maxwell"). Validated.
            arguments
                obj
                field (1,1) string
                value
            end
            parts = split(field, ".");
            if numel(parts) ~= 2 || ~isprop(obj, parts(1)) || parts(1) == "preset"
                NonGravForces.fail("Options", "inputField", ...
                    "unknown option '%s' -> use <group>.<field>, groups: drag srp erp thermal torques geom sun series partials", field);
            end
            grp = obj.(parts(1));
            if ~isfield(grp, parts(2))
                NonGravForces.fail("Options", "inputField", ...
                    "unknown option '%s' -> fields of %s are: %s", field, parts(1), strjoin(fieldnames(grp), ", "));
            end
            NonGravForces.Options.validateValue(field, value);
            if isstring(value) || ischar(value), value = string(value); end
            grp.(parts(2)) = value;
            obj.(parts(1)) = grp;
        end

        function s = toStruct(obj)
            %TOSTRUCT Plain struct (for provenance hashing and display).
            s = struct('preset', obj.preset, 'drag', obj.drag, 'srp', obj.srp, ...
                'erp', obj.erp, 'thermal', obj.thermal, 'torques', obj.torques, ...
                'geom', obj.geom, 'sun', obj.sun, 'series', obj.series, 'partials', obj.partials);
        end
    end
    methods (Access = private)
        function obj = applyPreset(obj, preset)
            switch preset
                case "fast"
                    dens = "exponential"; gsi = "constantCd"; acc = "fixed"; alpha = 0.85; ms = [];
                    ecl = "conical"; erpB = "none"; grid = [0 0]; th = "none"; ss = "none"; lutN = 16; sunS = "rs.analytic";
                case "standard"
                    dens = "msis00"; gsi = "sentman"; acc = "fixed"; alpha = 0.85; ms = [];
                    ecl = "conical"; erpB = "knocke"; grid = [12 24]; th = "steadyState"; ss = "lut"; lutN = 16; sunS = "rs.analytic";
                case "precise"
                    % The FITTED constant of March et al. (2021), not a computed accommodation, by owner
                    % ruling of 2026-09-17 (HANDOVER O-42). That paper is the state of the art for this
                    % exact application and fits a constant per mission BECAUSE the in-situ composition a
                    % computed model needs is not measured on these satellites; 0.85 is its recommendation
                    % and supersedes Doornbos (2011)'s 0.93. The adsorption models remain available by
                    % name -- and disagree with this by 17.5 % of the drag at 520 km in solar minimum,
                    % which is measured, not feared (HANDOVER Part G, LP-E).
                    dens = "msis21"; gsi = "sentman"; acc = "fixed"; alpha = 0.85; ms = [];
                    ecl = "solaarsCF"; erpB = "ceres"; grid = [24 48]; th = "steadyState"; ss = "lut"; lutN = 32; sunS = "rs.de440";
            end
            obj.drag     = struct('density', dens, 'wind', "corotation", 'gsi', gsi, ...
                                  'accommodation', acc, 'alpha', alpha, 'surfaceMass_amu', ms, ...
                                  'specularFraction', [], 'Tw_K', 300, 'cdScale', 1.0, 'cd', 2.2);
            obj.srp      = struct('tsi', NonGravForces.const().S0, 'eclipse', ecl, 'lunarEclipse', false, 'crScale', 1.0);
            obj.erp      = struct('backend', erpB, 'grid', grid);
            obj.thermal  = struct('model', th, 'internalPower_W', 0);
            obj.torques  = struct('gravityGradient', false, 'magnetic', false);
            obj.geom     = struct('selfShadow', ss, 'lutResolution', lutN);
            obj.sun      = struct('source', sunS);
            obj.series   = struct('useParallel', false, 'chunk', 8640);
            obj.partials = struct('method', "analytic");
        end
    end
    methods (Static, Access = private)
        function validateValue(field, value)
            sets = struct( ...
                'drag_density',       ["exponential","msis00","msis21"], ...
                'drag_wind',          ["none","corotation","hwm14"], ...
                'drag_gsi',           ["constantCd","sentman","maxwell"], ...
                'drag_accommodation', ["fixed","pilinski2010","pilinski2011","walker2014"], ...
                'srp_eclipse',        ["cylindrical","conical","solaarsCF"], ...
                'erp_backend',        ["none","pointSource","knocke","ceres"], ...
                'thermal_model',      ["none","steadyState","transient"], ...
                'geom_selfShadow',    ["none","lut","raycast"], ...
                'sun_source',         ["rs.analytic","rs.de440"], ...
                'partials_method',    ["analytic","finiteDiff"]);
            key = replace(field, ".", "_");
            if isfield(sets, key)
                if ~(isstring(value) || ischar(value)) || ~isscalar(string(value)) || ~any(string(value) == sets.(key))
                    NonGravForces.fail("Options", "inputValue", ...
                        "%s = '%s' not admissible -> one of %s", field, strjoin(string(value), " "), strjoin(sets.(key), ", "));
                end
                return
            end
            % numeric fields; an empty range row admits [] (the "not stated" of R5)
            ranges = struct('drag_alpha', [0 1 1], 'drag_surfaceMass_amu', [eps 1e6 1], 'drag_specularFraction', [0 1 1], ...
                'drag_Tw_K', [1 2000 0], 'drag_cdScale', [0.1 10 0], 'drag_cd', [0.5 10 0], 'srp_crScale', [0.1 10 0], ...
                'srp_tsi', [1300 1400 0], 'thermal_internalPower_W', [0 1e4 0], 'series_chunk', [1 1e7 0], ...
                'geom_lutResolution', [1 256 0]);
            if isfield(ranges, key)
                r = ranges.(key);
                if r(3) && isempty(value) && isnumeric(value)
                    return
                end
                if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= r(1) && value <= r(2))
                    NonGravForces.fail("Options", "inputValue", ...
                        "%s must be a finite scalar in [%g, %g]%s -> got %s", field, r(1), r(2), string(repmat(' or []', 1, r(3))), mat2str(value));
                end
                if key == "geom_lutResolution" && value ~= round(value)
                    NonGravForces.fail("Options", "inputValue", "%s must be an integer number of cells per cube-face edge -> got %g", field, value);
                end
                return
            end
            switch key
                case {'srp_lunarEclipse','torques_gravityGradient','torques_magnetic','series_useParallel'}
                    if ~(islogical(value) && isscalar(value))
                        NonGravForces.fail("Options", "inputValue", "%s must be a logical scalar -> use true/false", field);
                    end
                case 'erp_grid'
                    if ~(isnumeric(value) && isequal(size(value), [1 2]) && all(value >= 0) && all(value == round(value)))
                        NonGravForces.fail("Options", "inputValue", "%s must be [nRings nSectors] non-negative integers over the visible cap -> e.g. [12 24]", field);
                    end
            end
        end
    end
end
