classdef Options
%OPTIONS Validated option set for NonGravForces.Model, with fidelity presets.
%
%   opt = NonGravForces.Options                       % preset "standard"
%   opt = NonGravForces.Options("precise")
%   opt = NonGravForces.Options("fast", "drag.density", "exponential")
%   opt = opt.set("srp.tsi", 1367)                    % legacy TSI, explicit
%
%   Configure once, evaluate many (principle P5). Every field is validated
%   at construction; an unknown field name errors (rule R4: never degrade
%   silently). Presets fix the DEFAULT tier of every model family; extended
%   models (winds, CLL, IGRF torques, transient thermal) are opt-ins that no
%   preset enables, so preset results never move when a new model lands.
%
%   FIDELITY LADDER
%     field                 fast          standard        precise
%     drag.density          exponential   msis00          msis21
%     drag.wind             corotation    corotation      corotation (+hwm14 opt-in)
%     drag.gsi              constantCd    sentman         sentman (+cll opt-in)
%     drag.accommodation    0.85 fixed    0.85 fixed      sesam
%     srp.eclipse           conical       conical         solaarsCF
%     erp.backend           none          knocke          ceres
%     erp.grid              -             12x24           24x48
%     thermal.model         none          steadyState     steadyState (+transient)
%     geom.selfShadow       none          lut             raycast
%     sun.source            rs.analytic   rs.analytic     rs.de440
%
% PROPERTIES (read-only after construction; change through set())
%   preset      string    "fast" | "standard" | "precise"
%   drag        struct    .density .wind .gsi .accommodation .alpha .Tw_K
%                         .cdScale
%   srp         struct    .tsi [W/m^2] .eclipse .lunarEclipse (logical)
%   erp         struct    .backend .grid [1x2] .anisotropy
%   thermal     struct    .model .internalPower_W
%   torques     struct    .gravityGradient .magnetic (logical)
%   geom        struct    .selfShadow .lutResolution
%   sun         struct    .source
%   series      struct    .useParallel (logical) .chunk (double, epochs)
%   partials    struct    .method "analytic" | "complexStep" | "finiteDiff"
%
% ERRORS
%   NonGravForces:Options:inputPreset   unknown preset
%   NonGravForces:Options:inputField    unknown dotted field name
%   NonGravForces:Options:inputValue    value outside its admissible set
%
% EXAMPLE
%   opt = NonGravForces.Options("standard", "drag.wind", "hwm14");
%   opt.drag.wind
%
% See also NonGravForces.Model.
%
% KEYWORDS options presets fidelity validation
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
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
            %SET Override one dotted field, e.g. set("drag.gsi","cll"). Validated.
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
                case "fast",     dens = "exponential"; gsi = "constantCd"; acc = "fixed"; ecl = "conical";   erpB = "none";   grid = [0 0];   th = "none";        ss = "none";    sunS = "rs.analytic";
                case "standard", dens = "msis00";      gsi = "sentman";    acc = "fixed"; ecl = "conical";   erpB = "knocke"; grid = [12 24]; th = "steadyState"; ss = "lut";     sunS = "rs.analytic";
                case "precise",  dens = "msis21";      gsi = "sentman";    acc = "sesam"; ecl = "solaarsCF"; erpB = "ceres";  grid = [24 48]; th = "steadyState"; ss = "raycast"; sunS = "rs.de440";
            end
            obj.drag     = struct('density', dens, 'wind', "corotation", 'gsi', gsi, ...
                                  'accommodation', acc, 'alpha', 0.85, 'Tw_K', 300, 'cdScale', 1.0);
            obj.srp      = struct('tsi', NonGravForces.const().S0, 'eclipse', ecl, 'lunarEclipse', false);
            obj.erp      = struct('backend', erpB, 'grid', grid, 'anisotropy', "lambert");
            obj.thermal  = struct('model', th, 'internalPower_W', 0);
            obj.torques  = struct('gravityGradient', false, 'magnetic', false);
            obj.geom     = struct('selfShadow', ss, 'lutResolution', 642);
            obj.sun      = struct('source', sunS);
            obj.series   = struct('useParallel', false, 'chunk', 8640);
            obj.partials = struct('method', "analytic");
        end
    end
    methods (Static, Access = private)
        function validateValue(field, value)
            sets = struct( ...
                'drag_density',       ["exponential","msis00","msis21","dtm2020","jb2008"], ...
                'drag_wind',          ["none","corotation","hwm14"], ...
                'drag_gsi',           ["constantCd","sentman","cll","maxwell","schamberg"], ...
                'drag_accommodation', ["fixed","sesam","langmuir"], ...
                'srp_eclipse',        ["cylindrical","conical","solaarsCF"], ...
                'erp_backend',        ["none","pointSource","knocke","ceres"], ...
                'erp_anisotropy',     ["lambert","adm"], ...
                'thermal_model',      ["none","foldedIntoSrp","steadyState","transient"], ...
                'geom_selfShadow',    ["none","lut","raycast"], ...
                'sun_source',         ["rs.analytic","rs.de440"], ...
                'partials_method',    ["analytic","complexStep","finiteDiff"]);
            key = replace(field, ".", "_");
            if isfield(sets, key)
                if ~(isstring(value) || ischar(value)) || ~any(string(value) == sets.(key))
                    NonGravForces.fail("Options", "inputValue", ...
                        "%s = '%s' not admissible -> one of %s", field, string(value), strjoin(sets.(key), ", "));
                end
                return
            end
            ranges = struct('drag_alpha', [0 1], 'drag_Tw_K', [1 2000], 'drag_cdScale', [0.1 10], ...
                'srp_tsi', [1300 1400], 'thermal_internalPower_W', [0 1e4], 'series_chunk', [1 1e7], ...
                'geom_lutResolution', [12 1e5]);
            if isfield(ranges, key)
                r = ranges.(key);
                if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= r(1) && value <= r(2))
                    NonGravForces.fail("Options", "inputValue", ...
                        "%s must be a finite scalar in [%g, %g] -> got %s", field, r(1), r(2), mat2str(value));
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
                        NonGravForces.fail("Options", "inputValue", "%s must be [nLat nLon] non-negative integers -> e.g. [12 24]", field);
                    end
            end
        end
    end
end
