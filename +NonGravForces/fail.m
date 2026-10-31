function ME = fail(unit, mnemonic, message, varargin)
%FAIL Build (and by default throw) a NonGravForces error with a uniform identifier.
%
%   NonGravForces.fail(unit, mnemonic, message, args...)          throws
%   ME = NonGravForces.fail(unit, mnemonic, message, args...)     returns
%
%   Every error raised by the toolbox goes through this function, so that
%   (1) the identifier is always 'NonGravForces:<unit>:<mnemonic>', (2) the
%   message always ends with a CORRECTION HINT telling the caller what to
%   change, and (3) a test can verify the identifier without parsing text.
%
%   MNEMONIC TAXONOMY (the first token of the mnemonic names the class):
%     input*       contract violation: shape, class, range, unknown option
%     dep*         a dependency is absent, too old, or answered unexpectedly
%     data*        a data file is missing, malformed, or outside its range
%     numeric*     NaN/Inf produced, matrix singular, tolerance not reachable
%     unsupported* a valid request the toolbox refuses by design (no fallback)
%     notImplemented   a specified unit that is not yet delivered (stubs)
%
%   Use: if nargout == 0 the exception is thrown from the CALLER's frame
%   (throwAsCaller), so the stack points at the offending function, not at
%   fail itself.
%
% INPUT ARGUMENTS
%   Name      Size   Type    Units  Description
%   ----      ----   ----    -----  -----------
%   unit      1x1    string  -      Function or class name, e.g. "eclipseFactor"
%   mnemonic  1x1    string  -      camelCase reason, first token from the
%                                   taxonomy above, e.g. "inputShape"
%   message   1x1    string  -      sprintf-style message. Must contain the
%                                   correction hint after the marker ' -> '.
%   args      any    -       -      sprintf arguments for message
%
% OUTPUT ARGUMENTS
%   Name      Size   Type        Units  Description
%   ----      ----   ----        -----  -----------
%   ME        1x1    MException  -      Only when an output is requested; then
%                                       nothing is thrown.
%
% ERRORS
%   NonGravForces:fail:inputMnemonic   mnemonic does not start with a
%                                      taxonomy token
%   NonGravForces:fail:inputHint       message lacks the ' -> ' hint marker
%   MATLAB:validators:mustBeNonzeroLengthText   unit, mnemonic or message empty
%   These two identifiers are raised with error(), not through fail itself.
%
% ACCURACY
%   Not applicable.
%
% REFERENCES
%   HANDOVER principle P10 (one raise point); CLAUDE.md section 3 (error conventions).
%
% EXAMPLE
%   NonGravForces.fail("eclipseFactor", "inputShape", ...
%       "r must be 3xK, got %dx%d -> transpose the input", 4, 10)
%
% See also MException, throwAsCaller.
%
% KEYWORDS error handling identifier taxonomy
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
arguments
    unit (1,1) string {mustBeNonzeroLengthText}
    mnemonic (1,1) string {mustBeNonzeroLengthText}
    message (1,1) string {mustBeNonzeroLengthText}
end
arguments (Repeating)
    varargin
end
classes = ["input","dep","data","numeric","unsupported","notImplemented"];
if ~any(startsWith(mnemonic, classes))
    error("NonGravForces:fail:inputMnemonic", ...
        "mnemonic '%s' must start with one of %s -> pick the class that names the failure", ...
        mnemonic, strjoin(classes, ", "));
end
if ~contains(message, " -> ")
    error("NonGravForces:fail:inputHint", ...
        "message for %s:%s lacks the ' -> ' correction hint -> append what the caller should change", ...
        unit, mnemonic);
end
id  = "NonGravForces:" + unit + ":" + mnemonic;
txt = sprintf(message, varargin{:});
ME  = MException(id, "%s", txt);
if nargout == 0
    throwAsCaller(ME);
end
end
