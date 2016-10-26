function [Iterate,run_outer,step,penalty,currentTolMesh,currentTolCon, ...
    currentOmega,lambda,shift,lambdabar,alphaL,penaltyFactor,betaconstr, ...
    betamesh,alphamesh,alphaconstr,startTolCon,startOmega,startTolMesh, ...
    numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr,OuterIter,reduceFunCount] = ...
    psAugInit(Iterate,options,verbosity);
%PSAUGINIT Initial values of parameters for augmented Lagragian PS 
%   Private to PATTERNSEARCH

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/21 19:21:50 $


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
penaltyFactor = psoptimget(options,'PenaltyFactor',psoptimset,'fast');
% Parameter used to update currentTolMesh (successful subproblem)
betamesh  = 1.01; % 1.0 0.9 0.25
% Parameter used to update currentTolCon (successful subproblem)
betaconstr  = 0.25; %0.9, 0.25;
% Parameter used to update currentTolMesh (unsuccessful subproblem)
alphamesh = 1.0; % 1.0, 0.75, 0.1
% Parameter used to update currentTolCon (unsuccessful subproblem)
alphaconstr = 0.75; %1.0, 0.75, 0.1;
% Parameter used for constraint satisfaction after a sub-problem is solved
alphaL = min(1,alphaconstr/(1-alphaconstr)); % <= 1
% Constants used for tolerances
startTolCon  = 1.0;  % etas
startOmega   = 1.0;  % omega
% Initial mesh size should not be too small (multiply by this factor)
meshConstant = 1e3;
% Positive constant added to shift
shiftConstant = 1e-4;
% Initial values for penalty parameter
penalty = psoptimget(options,'InitialPenalty',psoptimset,'fast');

% Error checking on InitialPenalty and PenaltyFactor
if penalty < 1
    if verbosity > 0
        warning('gads:PSAUGINIT:smallPenalty','InitialPenalty must be greater than or equal to one; using default.\n');
    end
    penalty = psoptimget(psoptimset,'InitialPenalty');
end
if penaltyFactor <= 1
    if verbosity > 0
        warning('gads:PSAUGINIT:smallPenaltyFactor','PenaltyFactor must be greater than one; using default.\n');
    end
    penalty = psoptimget(psoptimset,'PenaltyFactor');
end

TolCon  = psoptimget(options,'TolCon',psoptimset,'fast');
minMesh = psoptimget(options,'TolMesh',psoptimset,'fast');
InitialMeshSize = psoptimget(options,'InitialMeshSize',psoptimset,'fast');
% Inverse of penalty (to be used in equations)
invPenalty = 1/penalty;
% Shifts s.t.  s_i - c_i(x) > 0; s_i > 0
shift(1:numNonlinIneqcstr) = max(0,Iterate.cineq) + shiftConstant;
shift(numNonlinIneqcstr+1:numNonlinCstr) = max(shiftConstant,Iterate.ceq); 
lambda = (shift./invPenalty).^(1/alphaL);
% Parameters to be changed in outer loop (main problem)
currentTolCon  = min(1e-1,startTolCon*(invPenalty^alphaconstr));
currentOmega   = min(1e-1,startOmega*(invPenalty^alphamesh));
currentTolMesh = max(minMesh,meshConstant*currentOmega/(1 + norm(lambda) + penalty));
currentTolMesh = min(currentTolMesh,InitialMeshSize);
startTolMesh = currentTolMesh;
OuterIter = 0;
run_outer = true;
step = '';

% Parameter to track function evaluation
reduceFunCount = 0;