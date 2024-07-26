function f = friction_factor(re, roughness, lQd, opts)
% friction_factor, calculate friction factor or k factor for pipe flow
%
% f = friction_factor(re) returns darcy friction factor for a smooth
% pipe at the given reynolds number.
%
% f = friction_factor(re, roughness) returns friction factor with the given
% relative surface rougness (e/D).
%
% k = friction_factor(re, roughness, lQd) returns the k factor for a given 
% L/D of pipe.
%
% f = friction_factor(..., name=value, ...) accepts optional name-value
% pairs.
%
% Options:
%  - implicit = true: Set to true (default) for the implicit colebrook
%  equation or set to false for the explicit Swamee-Jain approximation.
%  The explicit option is substantially faster.
%  - solverOpts = optimset(): Provide modified options structure for the
%  implicit fzero solution.  Not used for explicit.
%  - reTransition = [2000, 4000]: Reynolds number range for laminar to 
%  turbulent transition.  Can alter to force laminar or turbulent answer.
%  - smallLaminar = 0.1: Smoothing factor used between laminar and
%  transition region.
%  - smallTurb = 0.4: Smoothing factor used between transition and
%  turbulent region.

arguments
    re double {mustBePositiveOrNan};
    roughness double {mustBePositiveOrNan} = 0.0;
    lQd double {mustBePositive} = 1.0; 
    opts.implicit (1,1) = true;
    opts.solverOpts (1, 1) = optimset();
    opts.reTransition (2, 1) double = [2000, 4000];
    opts.smallLaminar (1, 1) double {mustBePositive} = 0.1;
    opts.smallTurb (1, 1) double {mustBePositive} = 0.4;
end

%% Calculate turbulent friction factor
% Use the Swamee-Jain equation for an approximation or guess
fturb = 0.25./log10(roughness./3.7 + 5.74./re.^0.9).^2;

% solve for the implicit values if needed
if opts.implicit
    f0 = fturb;
    % force expansion of all inputs to the same size
    re = re + zeros(size(roughness));
    roughness = roughness + zeros(size(re));
    for i = 1:numel(re)
        if isnan(re(i)) || isnan(roughness(i))
            % caon't run implicit solver for nan
            continue;
        end
        fun = @(f) colebrook(re(i), roughness(i), f);
        fturb(i) = fzero(fun, f0(i), opts.solverOpts);
    end
end

%% Calculate laminar to turbulent transition
 flaminar = 64./re;
 transFactor = (re - opts.reTransition(1))./(opts.reTransition(2)-opts.reTransition(1));
 transFactor = max_s(0.0, min_s(1.0, transFactor, opts.smallTurb), opts.smallLaminar);
 f = (1.0-transFactor).*flaminar + transFactor.*fturb;

% scale for lQd if needed, this will let the user get a k-factor instead of a friction factor
f = f.*lQd;

end

function err = colebrook(re, roughness, f)
% implementation of the implicit colebrook equation to pass to solver
rf = sqrt(f);
err = 1.0./rf + 2.0.*log10(roughness./3.7+2.51./re./rf);
end

function y = step_s(x, e)
% smooth approximation of step function
y = 0.5.*x./sqrt(x.^2 + e.^2) + 0.5;
end

function y = min_s(x1, x2, e)
% smooth approximation of min function
y = step_s(x1-x2, e).*x2 + step_s(x2-x1, e).*x1;
end

function y = max_s(x1, x2, e)
% smooth approximation of max function
y = step_s(x1-x2, e).*x1 + step_s(x2-x1, e).*x2;
end
