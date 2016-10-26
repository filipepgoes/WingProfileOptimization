function [X,Fval,maxConstr,reasonToStop,run_outer] = gaAugConverged(options,state,currentTolFun, ...
    Iterate,numNonlinCstr,numNonlinIneqcstr,lambdabar,X,FVAL,run_outer,step,type,infMessage,verbosity);
%GAAUGCONVERGED Augmented lagrangian convergence test.
%   Private to GA

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/21 19:21:33 $
%    Rakesh Kumar

X(:) = Iterate.x;
Fval = Iterate.f;
reasonToStop = '';
comp_slack = 0;
maxConstr = 0;
Gen = state.Generation;
stallTol = min(options.TolFun,eps);
if strcmpi(type,'nonlinearconstr')
    % Calculate complementary condition and constraint violation
    if numNonlinIneqcstr
        comp_slack = norm(Iterate.cineq.*lambdabar(1:numNonlinIneqcstr));
        maxConstr = max([maxConstr;Iterate.cineq(:)]);
    end
    if numNonlinCstr > numNonlinIneqcstr
        maxConstr = max([maxConstr;abs(Iterate.ceq(:))]);
    end
else % Linearly constrained problems do not enter in this function
    return;
end
% Print iterative information
if  verbosity > 1 && Gen > 0
    FunEval  = state.FunEval;
    MeanFval = mean(state.Score);
    StallGen = Gen - state.LastImprovement;
    fprintf('%5.0f       %5.0f  %12.6g %12.4g    %3.0f\n', ...
        Gen, FunEval, Fval, maxConstr, StallGen);
end
% Converged at an optimum
if currentTolFun < options.TolFun && maxConstr <= options.TolCon && ...
        comp_slack <= sqrt(options.TolCon)
    run_outer = false;
    reasonToStop = sprintf('%s','Optimization terminated: ');
    reasonToStop = [reasonToStop,sprintf('%s %6.5g%s', 'current tolerance on f(x)', currentTolFun, ' is less than options.TolFun')];
    % Check if linear constraints are satisfied
    linCon = options.LinearConstr;
    if linCon.linconCheck && ~feasibleLinearConstraints
        reasonToStop = [reasonToStop,sprintf('%s\n','but linear constraints are not satisfied.')];
    else
        reasonToStop = [reasonToStop,sprintf('\n%s', ' and constraint violation is less than options.TolCon.')];
    end
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Stall but constraints are satisfied
if currentTolFun < stallTol && maxConstr <= options.TolCon
    run_outer = false;
    reasonToStop = sprintf('%s %6.5g\n','Optimization terminated: norm of the step is less than ',eps);
    % Check if linear constraints are satisfied
    linCon = options.LinearConstr;
    if linCon.linconCheck && ~feasibleLinearConstraints
        reasonToStop = [reasonToStop,sprintf('%s\n','but linear constraints are not satisfied.')];
    else
        reasonToStop = [reasonToStop, sprintf('%s', ' and constraints violation is less than options.TolCon.')];
    end
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% fmincon enocntered NaN or Inf and could not continue; error here
if ~isempty(infMessage) & strmatch('optimlib:optimfcnchk',infMessage)
    msg = sprintf('%s\n','Constraint function returned non-real value;can not continue.');
    error('gads:GAAUGCONVERGED:NaNFval',msg);
end
% Stalls and constraints are not satisfied
if strcmpi(step,'Infeasible')
    run_outer  = false;
    reasonToStop = sprintf('%s\n','Optimization terminated: no feasible point found.');
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Maximum generations
if(Gen >= options.Generations)
    reasonToStop = sprintf(['Optimization terminated: ','maximum number of generations exceeded.']);
    run_outer  = false;
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Maximum time limit
if((cputime-state.StartTime) > options.TimeLimit)
    reasonToStop = sprintf(['Optimization terminated: ','time limit exceeded.']);
    run_outer  = false;
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Stall generation limit
if((cputime-state.LastImprovementTime) > options.StallTimeLimit)
    reasonToStop = sprintf(['Optimization terminated: ','stall time limit exceeded.']);
    run_outer  = false;
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Stall time limit
if((state.Generation  - state.LastImprovement) > options.StallGenLimit)
    reasonToStop = sprintf(['Optimization terminated: ','stall generations limit exceeded.']);
    run_outer  = false;
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Minimum fitness limit
if(min(min(state.Score)) <= options.FitnessLimit )
    reasonToStop = sprintf(['Optimization terminated: ','minimum fitness limit reached.']);
    run_outer  = false;
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end
% Stop requested from user
if(~isempty(state.StopFlag))
    reasonToStop = sprintf(['Optimization terminated: ',state.StopFlag]);
    run_outer  = false;
    if verbosity > 0
        fprintf('%s\n',reasonToStop);
    end
    return;
end

% Print header again
if verbosity > 1 && rem(Gen,20)==0 && Gen > 0
    fprintf('\n                           Best       max        Stall\n');
    fprintf('Generation  f-count        f(x)     constraint  Generations\n');
end
%---------------------------------------------------------
    function feasible  = feasibleLinearConstraints
        % Function to check if linear constraints are satisfied at final point
        % If it is a constrained problem and we want to check linear constraints
        A = linCon.A;
        L = linCon.L;
        U = linCon.U;
        IndEqcstr = linCon.IndEqcstr;
        tol = sqrt(options.TolCon);
        [FVAL,best] = min(state.Score);
        X = state.Population(best,:);
        feasible = isfeasible(Iterate.x,A,L,U,tol,IndEqcstr);
    end % End of feasibleLinearConstraints
end