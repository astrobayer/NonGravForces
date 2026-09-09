function out = packageToolbox(options)
%PACKAGETOOLBOX Build NonGravForces.mltbx from a programmatic options object.
%
%   out = packageToolbox
%   out = packageToolbox('Output', "/tmp/NonGravForces.mltbx", 'Build', false)
%
%   NO .prj FILE. MATLAB's toolbox packager normally keeps its settings in a
%   generated XML project file, and that file is a second authority: the version
%   would live there as well as in toolboxVersion.m, the file list there as well
%   as in this repository's layout, and the two would part company on the first
%   release nobody double-checked. So the options are built here, from the
%   sources that already own the facts, and the .prj is not committed.
%
%   WHAT IS SHIPPED, AND WHAT IS DELIBERATELY NOT.
%
%     shipped      +NonGravForces (the toolbox), demos, doc/html, mex (this
%                  project's own C), setupNonGravForces.m, LICENSE,
%                  THIRD-PARTY-NOTICES.md, README.md
%     NOT shipped  data/            third-party excerpts (IERS, CelesTrak,
%                                   NASA/JPL). Redistributing them inside a
%                                   package is the same question NRLMSISE-00
%                                   raised, and it is not this function's to
%                                   answer -- see obligation O-35
%                  mex/nrlmsise00   not ours to redistribute at all (O-34);
%                                   setupNonGravForces fetches it
%                  built MEX        platform-specific; the installer runs
%                                   NonGravForces.buildMex once per machine
%                  tests/, tools/   a released toolbox is not a development tree
%
%   THE CONSEQUENCE OF LEAVING data/ OUT IS THAT THE DEMOS DO NOT RUN FROM AN
%   INSTALLED COPY, and that is stated rather than discovered: they need the
%   2020-01-01 GRACE-C arc. Either the data ships or the demos carry a pointer
%   to where a reader downloads it. O-35 is that decision.
%
% NAME-VALUE ARGUMENTS
%   Name       Size  Type     Units  Description
%   ----       ----  ----     -----  -----------
%   Output     1x1   string   -      Destination .mltbx; default <root>/NonGravForces.mltbx
%   Build      1x1   logical  -      Actually package (default true). false returns the
%                                    options for inspection, which is what a dry run is for
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ------  -----  -----------
%   out   1x1   struct  -      .options (matlab.addons.toolbox.ToolboxOptions), .file (string),
%                              .built (logical), .version (string), .excluded (1xN string)
%
% ERRORS
%   NonGravForces:packageToolbox:depMissing   this MATLAB has no programmatic toolbox packager
%   NonGravForces:packageToolbox:dataMissing  a folder that must be shipped is not there
%
% EXAMPLE
%   out = packageToolbox('Build', false);
%   out.options.ToolboxFiles                  % what would go in, before anything is written
%
% See also NonGravForces.toolboxVersion, NonGravForces.buildMex, setupNonGravForces.
%
% KEYWORDS package mltbx release toolbox installer
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
arguments
    options.Output (1,1) string = ""
    options.Build (1,1) logical = true
end
unit = "packageToolbox";
if exist('matlab.addons.toolbox.ToolboxOptions', 'class') ~= 8
    NonGravForces.fail(unit, "depMissing", "this MATLAB has no matlab.addons.toolbox.ToolboxOptions -> it arrived in R2023a; package from the Add-On Manager instead, and do not commit the .prj it writes");
end
root = fileparts(fileparts(mfilename('fullpath')));                % tools/ -> repository root
ship = ["+NonGravForces", "demos", fullfile("doc", "html"), "mex", ...
        "setupNonGravForces.m", "LICENSE", "THIRD-PARTY-NOTICES.md", "README.md"];
paths = fullfile(root, ship);
absent = paths(~(isfolder(paths) | isfile(paths)));
if ~isempty(absent)
    NonGravForces.fail(unit, "dataMissing", "these must be shipped and are not present: %s -> run from a complete checkout", strjoin(absent, ", "));
end
v = string(NonGravForces.toolboxVersion('string'));
dest = options.Output;
if strlength(dest) == 0
    dest = string(fullfile(root, "NonGravForces.mltbx"));
end

opts = matlab.addons.toolbox.ToolboxOptions(root, "b5f0a1c2-9d3e-4a71-8c62-1f0d5e7a4b93", ...
    'ToolboxName', "NonGravForces", ...
    'ToolboxVersion', v, ...
    'OutputFile', dest);
opts.ToolboxFiles = paths;
opts.Description = "Non-gravitational forces and torques on GRACE-FO satellites: gas-surface interaction, solar and Earth radiation pressure, thermal recoil, and the sensors that measure them.";
opts.Summary = "Non-gravitational force modelling for low Earth orbit";
opts.AuthorName = "Matthias Weigelt";
opts.ToolboxGettingStartedGuide = fullfile(root, "doc", "html", "GettingStarted.html");
% Everything the toolbox needs on the path, and nothing that is not: tests and tools stay in the
% development tree, and mex/nrlmsise00 cannot be here because it was never ours (O-34).
opts.ToolboxMatlabPath = fullfile(root, ["", "demos"]);
excluded = ["data (third-party excerpts, O-35)", "mex/nrlmsise00 (not redistributable, O-34)", ...
            "built MEX binaries (per-machine; run NonGravForces.buildMex)", "tests, tools"];

out = struct('options', opts, 'file', dest, 'built', false, 'version', v, 'excluded', excluded);
if ~options.Build
    return
end
matlab.addons.toolbox.packageToolbox(opts);
out.built = isfile(dest);
if ~out.built
    NonGravForces.fail(unit, "dataMissing", "the packager reported success but %s does not exist -> check the output path is writable", dest);
end
end
