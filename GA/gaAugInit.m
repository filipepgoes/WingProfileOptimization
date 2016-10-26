function [Iterate,run_outer,step,penalty,minTolFun,currentTolFun,currentTolCon, ...
    currentOmega,lambda,shift,lambdabar,alphaL,penaltyFactor,betaconstr, ...
    betafun,alphafun,alphaconstr,startTolCon,startOmega,startTolFun, ...
    numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr,reduceFunCount,infMessage] = ... 
    gaAugInit(Iterate,options,verbosity);
%GAAUGINIT Initial values of parameters for augmented Lagragian GA 
%   Private to GA

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/21 19:21:34 $
%    Rakesh Kumar

% How many nonlinear constraints?
numNonlinIneqcstr = length(Iterate.cineq);
numNonlinEqcstr   = length(Iterate.ceq);
numNonlinCstr = numNonlinIneqcstr + numNonlinEqcstr;
% Initialize shift 
shift = zeros(numNonlinCstr,1);
% Initialize Lagrange estimates
lambda = zeros(numNonlinCstr,1);
% Initialize Lagrange estimate updates (to be used in AugUpdate)
lambdabar = zeros(numNonlinCstr,1);
% Penalty factor to increase penalty
penaltyFactor = gaoptimget(options,'PenaltyFactor',gaoptimset,'fast');
% Parameter used to update currentTolFun (successful subproblem)
betafun  = 1.01; % 1.0 0.9 0.25
% Parameter used to update currentTolCon (successful subproblem)
betaconstr  = 0.25; %0.9, 0.25;
% Parameter used to update currentTolFun (unsuccessful subproblem)
alphafun = 1.0; % 1.0 0.75 0.1
% Parameter used to update currentTolCon (unsuccessful subproblem)
alphaconstr = 0.75; %0.1, 0.75;
% Used to check constraint satisfaction after a sub-problem is solved
alphaL = min(1,alphaconstr/(1-alphaconstr)); % <= 1
% Initial value of tolerances
startTolCon  = 1.0;  % etas
startOmega   = 1.0;  % omega
% Initial tolerance on function value should not be too small (multiply by this factor)
funtolConstant = 1e3;
% Positive constant added to shift
shiftConstant = 1e-4;

% Initial values for changing parameters
% Recommended values are 1, 5, and 10  < 1 (rho)
penalty = gaoptimget(options,'InitialPenalty',gaoptimset,'fast');
% Error checking on rho and tau
if penalty < 1
    if verbosity > 0
        warning('gads:GAAUGINIT:smallPenalty','InitialPenalty must be greater than or equal to one; using default.\n');
    end
    penalty = gaoptimget(gaoptimset,'InitialPenalty');
end
if penaltyFactor <= 1
    if verbosity > 0
        warning('gads:GAAUGINIT:smallPenaltyFactor','PenaltyFactor must be greater than one; using default.\n');
    end
    penalty = gaoptimget(gaoptimset,'PenaltyFactor');
end

TolCon = gaoptimget(options,'TolCon',gaoptimset,'fast');
minTolFun = gaoptimget(options,'TolFun',gaoptimset,'fast');

% Inverse of penalty (to be used in equations)
invPenalty = 1/penalty;
% Shifts s.t. c_i(x) - s_i <0; s_i > 0
shift(1:numNonlinIneqcstr) = max(0,Iterate.cineq) + shiftConstant;
shift(numNonlinIneqcstr+1:numNonlinCstr) = max(shiftConstant,Iterate.ceq); 
lambda = (shift./invPenalty).^(1/alphaL);
% Parameters to be changed in outer loop (main problem)
currentTolCon  = startTolCon*(invPenalty^alphaconstr);
currentOmega   = startOmega*(invPenalty^alphafun);
% Lower bound on tolerance for first generation
currentTolFun =  max(minTolFun,currentOmega*funtolConstant);
% Upper bound on tolerance for first generation
currentTolFun = min(1,currentTolFun); 
startTolFun = currentTolFun;
run_outer = true;
step = '';

% Temporary parameters to track function evaluation
reduceFunCount = 0;
% Parameter to track if nonlcon returns non-real value to fmincon
infMessage = [];