function pool = dataPool(explicit, options)
%DATAPOOL The root of the shared data pool: the explicit argument, else DATAPOOL, else READING_DATA_POOL.
%
%   pool = NonGravForces.dataPool
%   pool = NonGravForces.dataPool(options.Pool, 'Option', "Pool")
%   pool = NonGravForces.dataPool(options.DataPool, 'Required', false)
%
%   EVERYTHING DOWNLOADED THAT IS DATA LIVES IN THE POOL, never in the toolbox (owner,
%   2026-11-08: "All data go always to the data pool. Only model code or libraries go to the
%   toolbox."). The pool is laid out by product and shared by every astrobayer toolbox, so one
%   copy of each product serves all of them:
%     <pool>/SPICE                              de440s.bsp (FetchKernels)
%     <pool>/EOP                                finals2000A.all
%     <pool>/SpaceWeather                       SW-All.csv, Dst/, F30/
%     <pool>/IGRF                               igrf14coeffs.txt
%     <pool>/CERES/EBAF                         CERES EBAF 4.2.1 rsut / rlut / rsdt
%     <pool>/EarthRadiationPressure             the GROOPS grids
%     <pool>/GRACEFO/L1B/RL04/<yyyy>/<yyyy-mm-dd>   Level-1B days, the reference day included
%     <pool>/ITSG, <pool>/Swarm, <pool>/TOLEOS  the observations of FetchPool
%   Model sources with the coefficient files they ship with, compiled MEX and the files committed
%   to the repository stay in the toolbox; they are not downloads of data.
%
%   THE ROOT IS RESOLVED IN ONE ORDER, here and nowhere else (R1):
%     1. explicit, when it is not empty (a caller's 'DataPool' or 'Pool');
%     2. getenv("DATAPOOL") (the toolbox installer writes it into startup.m);
%     3. getenv("READING_DATA_POOL") (the older name; still honoured, never required);
%     4. otherwise REFUSE. Nothing is ever written into the toolbox folder as a fallback.
%   NonGravForces keeps no stored preference of its own, so the contract's step "the
%   repository's own stored preference" is empty here. A value of blanks only counts as unset.
%   The folder is not required to exist: setup creates the product folders it fills, and a
%   reader checks the file it needs.
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   explicit  1x1   string  -      The caller's own pool option; "" (default) defers to the
%                                  environment.
%
% NAME-VALUE ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   Required  1x1   logical  -      true (default) refuses when no root is named; false returns ""
%                                   instead (setupNonGravForces, which REPORTS a skipped fetch and
%                                   never raises for absence).
%   Option    1x1   string   -      The caller's option name, quoted in the refusal (default
%                                   "DataPool").
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   pool  1x1   string  -      The root, surrounding blanks removed; "" only when Required is false
%                              and nothing is named.
%
% ACCURACY
%   Not applicable: no numerical work.
%
% ERRORS
%   NonGravForces:dataPool:dataPoolUnset   no explicit root, DATAPOOL and READING_DATA_POOL empty,
%                                          and Required true.
%
% REFERENCES
%   toolboxInstaller docs/DATA_POOL.md (the data-pool contract, 2026-11-08), sections 2 and 3;
%   HANDOVER Part I (2026-11-08).
%
% EXAMPLE
%   pool = NonGravForces.dataPool("", 'Required', false);
%   kernel = fullfile(NonGravForces.dataPool(), "SPICE", "de440s.bsp");
%
% See also setupNonGravForces, NonGravForces.fetchCatalog.
%
% KEYWORDS data pool DATAPOOL READING_DATA_POOL environment
% Created 2026-11-08 by Claude Opus 5.5 (Anthropic).
arguments
    explicit (1,1) string = ""
    options.Required (1,1) logical = true
    options.Option (1,1) string = "DataPool"
end
pool = strtrim(explicit);
if strlength(pool) == 0
    pool = strtrim(string(getenv('DATAPOOL')));
end
if strlength(pool) == 0
    pool = strtrim(string(getenv('READING_DATA_POOL')));
end
if strlength(pool) == 0 && options.Required
    NonGravForces.fail("dataPool", "dataPoolUnset", ...
        "no data pool is named: '%s' is empty and neither DATAPOOL nor READING_DATA_POOL is set -> " + ...
        "set the DATAPOOL environment variable (the toolbox installer's Data pool does) or pass '%s'", ...
        options.Option, options.Option);
end
end
