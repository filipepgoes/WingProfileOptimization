function [X,FVAL,REASON,OUTPUT,POPULATION,SCORES] = ...
    gacon(FitnessFcn,GenomeLength,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options)
%GACON Genetic algorithm generalized constrained solver.
%   Private function to GA
%
%  GACON solves problems of the form:
%           min F(X)    subject to:    A*X <= b
%            X                         Aeq*X = beq
%                                      LB <= X <= UB
%                                      C(X) <= 0
%                                      Ceq(X) = 0
%

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2.2.2 $  $Date: 2005/07/24 20:56:56 $
%    Rakesh Kumar

defaultopt = gaoptimset;

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(FitnessFcn,'defaults')
    X = defaultopt;
    return
end
% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
    error('gads:GACON:tenthInputNotStruct','Invalid input argument to GACON.');
elseif isempty(options)
    options = gaoptimset;
end
% Initialize output args
X = []; FVAL = []; REASON = [];POPULATION=[];SCORES=[];state = [];
user_options = options; EXITFLAG = [];
% Determine the verbosity
switch  gaoptimget(options,'Display',defaultopt,'fast')
    case {'off','none'}
        verbosity = 0;
    case 'final'
        verbosity = 1;
    case 'iter'
        verbosity = 2;
    case 'diagnose'
        verbosity = 3;
    otherwise
        verbosity = 1;
end
% Validate options and fitness function
[GenomeLength,FitnessFcn,options] = validate(GenomeLength,FitnessFcn,options);

% If nonlcon is a cell array with additional arguments, handle them
if iscell(nonlcon)
    conFcnArg = nonlcon(2:end);
    constr = nonlcon{1};
else
    conFcnArg = {};
    constr = nonlcon;
end
% Determine a start point
if ~isempty(options.InitialPopulation)
   initialX = options.InitialPopulation(1,:); 
else
   initialX = randn(1,GenomeLength);    
end
% Type of optimization sub-problem
X = initialX;
Iterate.x = initialX(:);
GenomeLength = length(Iterate.x);
if ~isempty(Aineq) || ~isempty(Aeq)
    subtype = 'linearconstraints';
elseif ~isempty(LB) || ~isempty(UB)
    subtype = 'boundconstraints';
else
    subtype = 'unconstrained';
end
% If we do not have any nonlinear constraint then there is no
% sub-problem; main problem type
if isempty(nonlcon)
    type = subtype; % Should not reach here
else
    type = 'nonlinearconstr';
end

% Remember the random number states used
OUTPUT.randstate  = rand('state');
OUTPUT.randnstate = randn('state');
OUTPUT.generations = 0;
OUTPUT.funccount   = 0;
OUTPUT.message   = '';

% Bound correction
[LB,UB,msg,EXITFLAG] = checkbound(LB,UB,GenomeLength,verbosity);
if EXITFLAG < 0
    FVAL     =    [];
    X(:)     =  Iterate.x;
    REASON = msg;
    OUTPUT.message = msg;
    if verbosity > 0
        fprintf('%s\n',msg)
    end
    return;
end
% Make A, L, U as L <= A*X <= U
[Iterate,A,L,U,nineqcstr,neqcstr,ncstr,IndIneqcstr,IndEqcstr,msg,EXITFLAG] = ...
    aluform(Iterate.x,Aineq,Bineq,Aeq,Beq,LB,UB,GenomeLength,subtype,verbosity, ...
    constr,X,conFcnArg);
if EXITFLAG < 0
    FVAL     =    [];
    X(:)     =  Iterate.x;
    REASON = msg;
    OUTPUT.message = msg;
    if verbosity > 0
        fprintf('%s\n',msg)
    end
    return;
else
   options.InitialPopulation(1,:) = Iterate.x'; 
end
% Validate constraints and add it to options structure
options = constrValidate(A,L,U,IndEqcstr,IndIneqcstr,nonlcon,options,subtype,type);
how = '';
% Initialize augmented lagrangian parameters
[Iterate,run_outer,step,penalty,minTolFun,currentTolFun,currentTolCon, ...
    currentOmega,lambda,shift,lambdabar,alphaL,penaltyFactor,betaconstr, ...
    betafun,alphafun,alphaconstr,startTolCon,startOmega,startTolFun, ...
    numNonlinIneqcstr,numNonlinEqcstr,numNonlinCstr,reduceFunCount,infMessage] = ...
    gaAugInit(Iterate,options,verbosity);
% Fitness funtion formulation as augmented lagrangian
if ~strcmpi(options.Vectorized,'on')
    SubFitness = @augLagFun;
else
    SubFitness = @augLagVecFun;
end
% Create initial state: population, scores, status data
state = makeState(GenomeLength,FitnessFcn,options);
Iterate.f = state.Score(1);
% Constraint information in state
state.NonlinIneq = Iterate.cineq;
state.NonlinEq =   Iterate.ceq;
state.Best(1) = Iterate.f;

% Give the plot/output Fcns a chance to do any initialization they need.
state = gaplot(FitnessFcn,options,state,'init');
[state,options] = gaoutput(FitnessFcn,options,state,'init');
% Print some diagnostic information if asked for
if verbosity > 2
    gadiagnose(FitnessFcn,nonlcon,GenomeLength,type,nineqcstr,neqcstr,ncstr,user_options);
end
% Setup display header
if  verbosity > 1
    fprintf('\n                           Best       max        Stall\n');
    fprintf('Generation  f-count        f(x)     constraint  Generations\n');
end

% Outer loop setup the augmeted Lagrangian formulation and inner loop
% solves the sub-problem
while run_outer
    % Check for augmented Lagrangian convergence
    [X,FVAL,maxConstr,REASON,run_outer] = gaAugConverged(options,state,currentTolFun, ...
        Iterate,numNonlinCstr,numNonlinIneqcstr,lambdabar,X,FVAL,run_outer, ...
        step,type,infMessage,verbosity);
    if ~run_outer, continue; end
    state.Generation = state.Generation + 1;
    % Update Iterate with value from sub-problem
    Iterate.f = SubFitness(Iterate.x');
   % Reset some parameters before solving sub-problem
    [exitFlag,innerState,options,step] = gaAugReset(Iterate,state, ...
        options,step,currentTolFun);
    % Solve the sub-problem
    while isempty(exitFlag)
        innerState.Generation = innerState.Generation + 1;
        offset = 0; totalPop = options.PopulationSize;
        % Adjust function evaluation counter
        innerState.FunEval = innerState.FunEval-reduceFunCount;
        reduceFunCount = 0;
        % Each sub-population loop
        for pop = 1:length(totalPop)
            populationSize =  totalPop(pop);
            thisPopulation = 1 + (offset:(offset + populationSize - 1));
            population = innerState.Population(thisPopulation,:);
            score = innerState.Score( thisPopulation );
            % Empty population is also possible
            if isempty(thisPopulation), continue;  end
            [score,population,innerState] = stepGA(score,population,options, ...
                innerState,GenomeLength,SubFitness);
            % store the results for this sub-population
            innerState.Population(thisPopulation,:) = population;
            innerState.Score(thisPopulation) = score;
            offset = offset + populationSize;
        end
        scores = innerState.Score;
        % Remember the best score
        best = min(innerState.Score);
        generation = innerState.Generation;
        innerState.Best(generation) = best;
        % Keep track of improvement in the best
        if((generation > 1) && finite(best))
            if(innerState.Best(generation-1) > best)
                innerState.LastImprovement = generation;
                innerState.LastImprovementTime = cputime;
            end
        end
        % Do any migration
        innerState = migrate(FitnessFcn,GenomeLength,options,innerState);
        % Check to see if any stopping criteria have been met for
        % sub-propblem
        exitFlag = isItTimeToStop(options,innerState);
    end %End sub-problem
    % Find the best solution found for sub-problem
    [unused,best] = min(innerState.Score);
    Iterate.x = innerState.Population(best,:)';
    % Update lagrange parameters and state; make sure that intermediate
    % point is feasible for next sub-problem
    updateIntermediateSolution;
    % Update the state (similar to innerState)
    state.Population = repmat(Iterate.x',options.PopulationSize,1);
    state.Score = repmat(Iterate.f,options.PopulationSize,1);
    scores = state.Score;
    best = min(state.Score);
    generation = state.Generation;
    state.Best(generation) = best;
    state.NonlinIneq = Iterate.cineq;
    state.NonlinEq =   Iterate.ceq;
    state.FunEval = state.FunEval + innerState.FunEval;
    % Keep track of improvement in the best
    if((generation > 1) && finite(best))
        if(state.Best(generation-1) ~= best)
            state.LastImprovement = generation;
            state.LastImprovementTime = cputime;
        end
    end
    % Update plot and output functions
    state = gaplot(FitnessFcn,options,state,'iter');
    [state,options,optchanged] = gaoutput(FitnessFcn,options,state,'iter');
    if optchanged
        % Determine the verbosity
        switch  gaoptimget(options,'Display',defaultopt,'fast')
            case {'off','none'}
                verbosity = 0;
            case 'final'
                verbosity = 1;
            case 'iter'
                verbosity = 2;
            case 'diagnose'
                verbosity = 3;
            otherwise
                verbosity = 1;
        end
        % Add constraint related options
        options = constrValidate(A,L,U,IndEqcstr,IndIneqcstr,nonlcon,options,subtype,type);
    end
end % End outer while loop

% Update output structure
OUTPUT.generations = state.Generation;
OUTPUT.message     = REASON;
OUTPUT.funccount   = state.FunEval;
OUTPUT.maxconstraint = maxConstr;
% Load up outputs
if(nargout > 4)
    POPULATION = state.Population;
    if(nargout > 5)
        SCORES = state.Score;
    end
end
% Call hybrid function
if(strcmpi(options.PopulationType,'doubleVector') && ~isempty(options.HybridFcn))
    [X,FVAL] = callHybridFunction;
end
% Give the Output functions a chance to finish up
gaplot(FitnessFcn,options,state,'done');
gaoutput(FitnessFcn,options,state,'done');

%----------------------------------------------------------------
% Function to update Lagrange multipliers and penalty. Also create next
% sub-problem and make sure that intermediate solution is feasible for next
% sub-problem
    function updateIntermediateSolution
        % Make sure that the current point is feasible for the sub-problem
        [Iterate.cineq(:),Iterate.ceq(:)] = feval(constr,Iterate.x',conFcnArg{:});
        % Reset TolFun in options which was changed for previous sub-problem
        options.TolFun = minTolFun;
        % Update lagrange and penalty parameters
        [lambda,lambdabar,penalty,currentOmega,currentTolFun,currentTolCon,shift,how] = ...
            gaAugUpdate(Iterate,lambda,penalty,startOmega,startTolCon,startTolFun, ...
            currentOmega,currentTolFun,currentTolCon,options.TolFun,options.TolCon, ...
            shift,penaltyFactor,betafun,betaconstr,alphafun,alphaconstr,alphaL, ...
            numNonlinIneqcstr,numNonlinCstr);
        if numNonlinIneqcstr
            shiftedConstr = -Iterate.cineq + shift(1:numNonlinIneqcstr);
            % Infeasible next point?
            if any(shiftedConstr <= 0) %&& strcmpi(how,'Increase penalty') 
                e = -2; f = Inf; infMessage = [];
                % Solve a problem to find a feasible solution
                warnstate = warning;
                warning off;
                try
                    [x,f,e] = fmincon(@(x) x(end),[Iterate.x;penaltyFactor], ...
                        [Aineq zeros(nineqcstr,1)],Bineq,[Aeq zeros(neqcstr,1)],Beq,[LB;0],[UB;penaltyFactor], ...
                        @infeasConstr,optimset('LargeScale','off','FunValCheck','on','TolFun',0,'Display','off'));
                catch
                    [unused,infMessage] = lasterr;
                end
                warning(warnstate);
                % Accept this condition as feasible iterate for next
                % iteration if f is close to 0
                if e > 0 && abs(f) >= 1e-3
                    f = 0;
                end
                % If exitflag is negative or 'f' is not close to zero,
                % then the problem does not have a feasible solution
                if e < 0 || abs(f) >= 1e-3
                    step = 'Infeasible';
                    Iterate.f = feval(FitnessFcn,Iterate.x',options.FitnessFcnArgs{:});
                    return;
                else
                    Iterate.x(:) = x(1:end-1);
                end
            end
        end
        % Compute the constraints at the intermediate point and function value
        [Iterate.cineq(:),Iterate.ceq(:)] = feval(constr,Iterate.x',options.NonconFcnArgs{:});
        % This should not happen
        if numNonlinIneqcstr
            shiftedConstr = -Iterate.cineq + shift(1:numNonlinIneqcstr);
            if any (shiftedConstr <= 0)
                % Adjust shift
                shift(1:numNonlinIneqcstr) = max(0,Iterate.cineq) + 1e-4;
                how = 'Infeasible point';
            end
        end
        Iterate.f = feval(FitnessFcn,Iterate.x',options.FitnessFcnArgs{:});
    end % End of updateIntermediateSolution
%---------------------------------------------------------------
% Nested function to model augmented Lagrangian
    function lagrangeFval = augLagFun(input,varargin)
        ceq = zeros(numNonlinEqcstr,1);
        cin = zeros(numNonlinIneqcstr,1);
        [cin(:),ceq(:)] = feval(options.NonconFcn,input,options.NonconFcnArgs{:});
        % Inequality constraint must satisfy this condition (see log term
        % in the augmented lagrangian formulation)
        if numNonlinIneqcstr
            shiftedConstr = -cin + shift(1:numNonlinIneqcstr); % must be > 0
            if any(shiftedConstr <= 0)
                lagrangeFval = NaN;
                reduceFunCount = reduceFunCount + 1;
                return;
            end
            fval = feval(FitnessFcn,input,options.FitnessFcnArgs{:});
            % lagrangeFval = f(x) - sum(lambda_i*shift_i*log(shift_i - c(x)_i))
            % i: nonlinear inequality constraints
            lagrangeFval = fval - sum(lambda(1:numNonlinIneqcstr).*shift(1:numNonlinIneqcstr).*log(shiftedConstr));
            if numNonlinEqcstr
                % lagrangeFval = lagrangeFval + sum(lambda_i*c(x)_i) +
                % (penalty/2)*sum(c(x)_i^2); i: equality constr
                lagrangeFval = lagrangeFval + sum(lambda(numNonlinIneqcstr+1:numNonlinCstr).*ceq) ...
                    + sum(ceq.^2 )*(penalty/2);
            end
        else
            % lagrangeFval = sum(lambda_i*c(x)_i) +
            % (penalty/2)*sum(c(x)_i^2); i: equality constr
            fval = feval(FitnessFcn,input,options.FitnessFcnArgs{:});
            lagrangeFval = fval + sum(lambda.*ceq) ...
                + sum(ceq.^2 )*(penalty/2);
        end
    end % End of augLagFcn
%---------------------------------------------------------------
% Vectorized version of augmented lagrangian formulation
    function phi = augLagVecFun(input,varargin)
        % How many points are to be evaluated?
        [nPoints,unused] = size(input);
        phi = nan(nPoints,1);
        % We are not doing any vectorized evaluation at this time for
        % nonlinear consraints
        for i = 1:nPoints
            phi(i) =  augLagFun(input(i,:));
        end
    end % End of augLagVecfun

%---------------------------------------------------------------
% Constraint formulation to handle infeasiblity problem
    function [cin,ceq] = infeasConstr(input,varargin)
        x = input(1:end-1);
        cin = zeros(numNonlinIneqcstr,1);
        [cin(:),ceq] = feval(options.NonconFcn,x',options.NonconFcnArgs{:});
        cin = cin - input(end)*shift(1:numNonlinIneqcstr);
    end % End of infeasConstr
%-----------------------------------------------------------------
% Hybrid function
    function [xhybrid,fhybrid] = callHybridFunction
        xhybrid = X;
        fhybrid = FVAL;
        %Who is the hybrid function
        if isa(options.HybridFcn,'function_handle')
            hfunc = func2str(options.HybridFcn);
        else
            hfunc = options.HybridFcn;
        end
        %Inform about hybrid scheme
        if   verbosity > 1
            fprintf('%s%s%s\n','Switching to the hybrid optimization algorithm (',upper(hfunc),').');
        end
        % Create functions handle to be passed to hybrid function
        FitnessHybridFcn = @(x) FitnessFcn(x,options.FitnessFcnArgs{:});
        if ~isempty(constr)
            ConstrHybridFcn  = @(x) constr(x,conFcnArg{:});
        else
            ConstrHybridFcn = [];
        end
        %Determine which syntax to call
        switch hfunc
            case 'fmincon'
                [xx,ff,e,o] = feval(options.HybridFcn,FitnessHybridFcn,X,Aineq, ...
                    Bineq,Aeq,Beq,LB,UB,ConstrHybridFcn,options.HybridFcnArgs{:});
                OUTPUT.funccount = OUTPUT.funccount + o.funcCount;
                OUTPUT.message   = [OUTPUT.message sprintf('\nFMINCON: \n'), o.message];
            case 'patternsearch'
                [xx,ff,e,o] = feval(options.HybridFcn,FitnessHybridFcn,X,Aineq, ...
                    Bineq,Aeq,Beq,LB,UB,ConstrHybridFcn,options.HybridFcnArgs{:});
                OUTPUT.funccount = OUTPUT.funccount + o.funccount;
                OUTPUT.message   = [OUTPUT.message sprintf('\nPATTERNSEARCH: \n'), o.message];
            case {'fminsearch', 'fminunc'}
                msg = sprintf('%s is unconstrained optimization solver',upper(hfunc));
                msg = [msg, sprintf('\n%s',' using constrained solver FMINCON as hybrid function.')];
                warning('gads:GACON:unconstrainedHybridFcn',msg);
                % We need to store this warning state so that we can
                % display it in the GUI
                [lastmsg, lastid] = lastwarn; warning off;
                [xx,ff,e,o] = feval(@fmincon,FitnessHybridFcn,X,Aineq, ...
                    Bineq,Aeq,Beq,LB,UB,ConstrHybridFcn,options.HybridFcnArgs{:});
                warning on; lastwarn(lastmsg,lastid);
                OUTPUT.funccount = OUTPUT.funccount + o.funcCount;
                OUTPUT.message   = [OUTPUT.message sprintf('\nFMINCON: \n'), o.message];
           otherwise
                error('gads:GA:hybridFcnError','Hybrid function must be @FMINCON, or @PATTERNSEARCH.')
        end
        % Check for exitflag and fval
        if ff < fhybrid && e > 0
            fhybrid = ff;
            xhybrid = xx;
        end
        %Inform about hybrid scheme termination
        if  verbosity > 1
            fprintf('%s%s\n',upper(hfunc), ' terminated.');
        end
    end % End of callHybridFunction
end  % End of gaconstr.m

