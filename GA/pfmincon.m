function [X,FVAL,EXITFLAG,OUTPUT] = pfmincon(FUN,initialX,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options)
%PFMINCON Finds a constrained minimum of a function.
%   PFMINCON solves problems of the form:
%           min F(X)    subject to:    A*X <= b
%            X                         Aeq*X = beq
%                                      LB <= X <= UB
%                                      C(X) <= 0
%                                      Ceq(X) = 0
%

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/06/21 19:21:46 $
%    Rakesh Kumar

defaultopt = psoptimset;
% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(FUN,'defaults')
    X = defaultopt;
    return
end

% If FUN is a cell array with additional arguments, handle them
if iscell(FUN)
    objFcnArg = FUN(2:end);
    FUN = FUN{1};
else
    objFcnArg = {};
    FUN = FUN;
end

% If NONLCON is a cell array with additional arguments, handle them
if iscell(nonlcon)
    conFcnArg = nonlcon(2:end);
    nonlcon = nonlcon{1};
else
    conFcnArg = {};
end

% Only function_handle or inlines are allowed
if isempty(FUN) ||  ~(isa(FUN,'inline') || isa(FUN,'function_handle'))
    error('gads:PFMINCON:needHandleOrInline','Objective function must be a function handle.');
end

if isempty(nonlcon) ||  ~(isa(nonlcon,'inline') || isa(nonlcon,'function_handle'))
    error('gads:PFMINCON:needHandleOrInline','Constraint function must be a function handle.');
end

% Initialize output args
X = []; FVAL = []; EXITFLAG = [];
OUTPUT = struct('function',FUN,'iterations',0,'funccount',0,'message',' ');

if(~isempty(initialX))
    Iterate.x = initialX(:);
    X = initialX;
    numberOfVariables = length(Iterate.x);
    type ='nonlinearconstr';
    if ~isempty(Aineq) || ~isempty(Aeq)
        subtype = 'linearconstraints';
    elseif ~isempty(LB) || ~isempty(UB)
        subtype = 'boundconstraints';
    else
        subtype = 'unconstrained';
    end
else
    error('gads:PFMINCON:initialPoint','You must provide an initial point.');
end
% Retrieve all the options
[verbosity,MeshExpansion,MeshContraction,Completesearch,MeshAccelerator,minMesh,MaxMeshSize, ...
    maxIter,maxFun,TolCon,TolBind,TolFun,TolX,MeshSize,pollmethod, ...
    pollorder,Completepoll,outputTrue,OutputFcns,OutputFcnArgs,plotTrue,PlotFcns, ...
    PlotFcnArgs,PlotInterval,searchtype,searchFcnArg,Cache,Vectorized, ...
    NotVectorizedPoll,NotVectorizedSearch,cachetol,cachelimit,scaleMesh, ...
    RotatePattern,TimeLimit]  =  checkoptions(options,defaultopt,numberOfVariables);
% Bound correction
[LB,UB,msg,EXITFLAG] = checkbound(LB,UB,numberOfVariables,verbosity);
if EXITFLAG < 0
    FVAL     =    [];
    X(:)     =  Iterate.x;
    OUTPUT.message = msg;
    if verbosity > 0
        fprintf('%s\n',msg)
    end
    return
end
% Make A, L, U and check for constraints
[Iterate,A,L,U,nineqcstr,neqcstr,ncstr,IndIneqcstr,IndEqcstr,msg,EXITFLAG] = ...
    aluform(Iterate.x,Aineq,Bineq,Aeq,Beq,LB,UB,numberOfVariables,subtype, ...
    verbosity,nonlcon,X,conFcnArg);
if EXITFLAG < 0
    FVAL     =    [];
    X(:)     =  Iterate.x;
    OUTPUT.message = msg;
    if verbosity > 0
        fprintf('%s\n',msg)
    end
    return
end
% Get initial values of some parameters
[FUN,Iterate,Iter,FunEval,scale,Successdir,innerNextIterate,deltaF,deltaX, ...
    MeshCont,NewMeshSize,infMessage,how,stopOutput,stopPlot,run_inner,OUTPUT, ...
    EXITFLAG,X,FVAL,StartTime] = getinitial(FUN,X,Iterate,Vectorized,objFcnArg, ...
    type,neqcstr,MeshContraction,MeshSize,scaleMesh,numberOfVariables,LB,UB);
% Initialize augmented Lagrangian related parameters
[Iterate,run_outer,step,penalty,currentTolMesh,currentTolCon,currentOmega, ...
    lambda,shift,lambdabar,alphaL,penaltyFactor,betaconstr,betamesh, ...
    alphamesh,alphaconstr,startTolCon,startOmega,startTolMesh,numNonlinIneqcstr, ...
    numNonlinEqcstr,numNonlinCstr,OuterIter,reduceFunCount] = psAugInit(Iterate,options,verbosity);
% Objective funtion formulation as augmented Lagrangian
if NotVectorizedPoll || NotVectorizedSearch
    SubObjective = @augLagVecFun;
else
    SubObjective = @augLagFun;
end
% Set up output & plot function
if (outputTrue || plotTrue)
    callOutputPlotFunctions('init')
end
% Print some diagnostic information if verbosity > 2
if verbosity > 2
    psdiagnose(FUN,nonlcon,Iterate,initialX,type,nineqcstr,neqcstr,ncstr,options);
end
% Setup display header
if  verbosity>1
    fprintf('\n                                  max\n');
    fprintf('  Iter   f-count      f(x)      constraint   MeshSize      Method\n');
end
% Outer loop setup the augmeted Lagrangian formulation and inner loop
% solves the sub-problem
while run_outer
    % Check for augmented Lagrangian convergence
    [X,EXITFLAG,FVAL,maxConstr,msg,run_outer] = psAugConverged(stopOutput,stopPlot, ...
        verbosity,OuterIter,maxIter,FunEval,maxFun,currentTolMesh,minMesh, ...
        infMessage,Iterate,deltaX,deltaF,numNonlinCstr,numNonlinIneqcstr, ...
        lambdabar,TolCon,TolFun,TolX,X,EXITFLAG,FVAL,run_outer,step,how, ...
        StartTime,TimeLimit);
    if ~run_outer, break; end
    % Store X before inner iteration starts (used to calculate deltaX)
    outerX =  Iterate.x;
    % Create a new variable innerIterate for the sub-problem
    innerIterate.x = Iterate.x;
    innerIterate.f = SubObjective(reshapeinput(X,innerIterate.x));
    % Set current Fun and X tolerances same as currentTolMesh for
    % sub-problem
    currentTolX = currentTolMesh; currentTolFun = currentTolMesh;

    % Reset some parameters before solving sub-problem
    [Iter,innerMaxIter,Successdir,deltaF,deltaX,MeshSize,EXITFLAG, ...
        run_inner,step] = psAugReset(options,numberOfVariables,step);
    % Solve the sub-problem
    while run_inner
        FunEval = FunEval - reduceFunCount; reduceFunCount = 0;
        % Check for convergence of sub-problem (augmented Lagrangian
        % formulation)
        [X,EXITFLAG,FVAL,msg,run_inner] = isconverged(stopOutput,stopPlot,0,Iter, ...
            innerMaxIter,FunEval,maxFun,pollmethod,MeshSize,currentTolMesh, ...
            infMessage,innerNextIterate,how,deltaX,deltaF,TolCon,currentTolFun,currentTolX,X, ...
            EXITFLAG,FVAL,run_inner,StartTime,TimeLimit);
        if ~run_inner,  break; end
        % SEARCH.
        [successSearch,innerNextIterate,FunEval] = search(SubObjective,[],X, ...
            searchtype,Completesearch,innerIterate,Successdir,pollorder,MeshSize, ...
            scale,TolBind,A,L,U,IndIneqcstr,IndEqcstr,Iter,FunEval,maxFun, ...
            subtype,NotVectorizedSearch,Cache,cachetol,cachelimit, ...
            searchFcnArg,{},objFcnArg);
        % POLL
        if ~successSearch  % Unsuccessful search
            [successPoll,innerNextIterate,FunEval,Successdir] = poll(SubObjective, ...
                [],X,pollmethod,Completepoll,pollorder,innerIterate,Successdir, ...
                MeshSize,scale,TolBind,A,L,U,IndIneqcstr,IndEqcstr,Iter, ...
                FunEval,maxFun,subtype,NotVectorizedPoll,Cache,cachetol, ...
                cachelimit,{},objFcnArg);
        else
            successPoll = false;
        end
        %Scale the variables in every iterations
        if strcmpi(scaleMesh,'on') && ~neqcstr
            scale = logscale(LB,UB,mean(innerIterate.x,2));
        end
        % Update
        [NewMeshSize,MeshContraction,how,deltaX,deltaF,scale,innerIterate,X,Iter, ...
            infMessage] = updateparam(pollmethod,successPoll,successSearch, ...
            MeshAccelerator,RotatePattern,MaxMeshSize,minMesh,MeshExpansion, ...
            MeshCont,MeshContraction,MeshSize,scale,innerNextIterate,innerIterate,X, ...
            Iter,how,infMessage);
        % Update mesh size
        MeshSize = NewMeshSize;
    end % End of sub-problem 
    % Update lagrange parameters
    updateIntermediateSolution;
    % Call output/plot functions
    if (outputTrue || plotTrue)
        callOutputPlotFunctions('iter')
    end
end % End of outer while

% Call output/plot functions
if (outputTrue || plotTrue)
    callOutputPlotFunctions('done')
end

OUTPUT = struct('function',FUN,'problemtype',type, ...
    'pollmethod',pollmethod,'searchmethod',searchtype,'iterations',OuterIter, ...
    'funccount',FunEval,'meshsize',currentTolMesh,'maxconstraint',maxConstr,'message',msg);
%---------------------------------------------------------------
% Function to update Lagrange multipliers and penalty. Also create next
% sub-problem and make sure that intermediate solution is feasible for next
% sub-problem
    function updateIntermediateSolution
        % Update Iterate.x from inner iteration
        Iterate.x = innerIterate.x;
        % Compute the constraints at the intermediate point but no need to
        % calculate function value
        [Iterate.cineq(:),Iterate.ceq(:)] = feval(nonlcon,reshapeinput(X,Iterate.x),conFcnArg{:});
        % Update Lagrange multipliers and penalty parameters (based on
        % constraint values only)
        [lambda,lambdabar,penalty,currentOmega,currentTolMesh,currentTolCon,shift, ...
            OuterIter,how] = psAugUpdate(Iterate,lambda,penalty,startOmega,startTolMesh, ...
            startTolCon,currentOmega,currentTolMesh,currentTolCon, ...
            minMesh,TolCon,shift,penaltyFactor,betamesh,betaconstr,alphamesh,alphaconstr, ...
            alphaL,numNonlinIneqcstr,numNonlinCstr,OuterIter);
        if numNonlinIneqcstr
            shiftedConstr = -Iterate.cineq + shift(1:numNonlinIneqcstr);
            % Infeasible point for next sub-problem?
            if any(shiftedConstr <= 0) %&& strcmpi(how,'Increase penalty')
                e = -2; f= Inf; infMessage = [];
                % Find a feasible solution
                warnstate = warning;
                warning off;
                try
                    [x,f,e] = fmincon(@(x) x(end),[Iterate.x;penaltyFactor], ...
                        [Aineq zeros(nineqcstr,1)],Bineq,[Aeq zeros(neqcstr,1)],Beq,[LB;0],[UB;1], ...
                        @infeasConstr,optimset('LargeScale','off','FunValCheck','on','TolFun',eps,'Display','off'));
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
                if (e < 0 || abs(f) >= 1e-3)
                    step = 'Infeasible';
                    Iterate.f = feval(FUN,reshapeinput(X,Iterate.x),objFcnArg{:});
                    return;
                else % Update the current Iterate
                    Iterate.x = x(1:end-1);
                end
            end
        end
        % Compute the constraints at the intermediate point and function value
        [Iterate.cineq(:),Iterate.ceq(:)] = feval(nonlcon,reshapeinput(X,Iterate.x),conFcnArg{:});
        % This should not happen (sanity check)
        if numNonlinIneqcstr
            shiftedConstr = -Iterate.cineq + shift(1:numNonlinIneqcstr);
            if any (shiftedConstr <= 0)
                % Adjust shift
                shift(1:numNonlinIneqcstr) = max(0,Iterate.cineq) + 1e-4;
                how = 'Infeasible point';
            end
        end
        Iterate.f = feval(FUN,reshapeinput(X,Iterate.x),objFcnArg{:});
        % Compute deltaX 
        deltaX = norm(outerX - Iterate.x);
    end % End of updateIntermediateSolution
%---------------------------------------------------------------
% Augmented Lagrangian sub-problem objective function formulation
    function lagrangeFval = augLagFun(input,varargin)
        ceq = zeros(numNonlinEqcstr,1);
        cin = zeros(numNonlinIneqcstr,1);
        [cin(:),ceq(:)] = feval(nonlcon,input,conFcnArg{:});
        % Inequality constraint must satisfy this condition (see log term
        % in the augmented lagrangian formulation)
        if numNonlinIneqcstr
            shiftedConstr = -cin + shift(1:numNonlinIneqcstr); % must be > 0
            if any(shiftedConstr <= 0)
                lagrangeFval = NaN; % This is okay because PS and GA work on a population
                reduceFunCount = reduceFunCount + 1;
                return;
            end
            fval = feval(FUN,input,objFcnArg{:}); 
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
            fval = feval(FUN,input,objFcnArg{:});
            lagrangeFval = fval + sum(lambda.*ceq) ...
                + sum(ceq.^2 )*(penalty/2);
        end
    end % End of augLagFcn
%---------------------------------------------------------------
% Vectorized version of augmented lagrangian formulation
    function phi = augLagVecFun(input,varargin)
        % Determine the size of input; row or column major
        [input_row,input_col] = size(initialX);
        if input_row >= input_col
            % How many points are to be evaluated?
            [unused,nPoints] = size(input);
            phi = nan(nPoints,1);
            for i = 1:nPoints
                phi(i) =  augLagFun(input(:,i));
            end
        else
            % How many points are to be evaluated?
            [nPoints,unused] = size(input);
            phi = nan(nPoints,1);
            for i = 1:nPoints
                phi(i) =  augLagFun(input(i,:));
            end
        end
    end % End of augLagVecfun
%---------------------------------------------------------------
% Constraint formulation to handle infeasiblity problem
    function [cin,ceq] = infeasConstr(input,varargin)
        x = input(1:end-1);
        cin = zeros(numNonlinIneqcstr,1);
        [cin(:),ceq] = feval(nonlcon,reshapeinput(X,x),conFcnArg{:});
        cin = cin - input(end)*shift(1:numNonlinIneqcstr);
    end % End of infeasConstr
%---------------------------------------------------------------
% Nested function to call output/plot functions
    function callOutputPlotFunctions(state)
        switch state
            case 'init'
                optimvalues.x = X; optimvalues.iteration = OuterIter; optimvalues.fval = Iterate.f;
                optimvalues.meshsize = currentTolMesh;optimvalues.funccount = FunEval;
                optimvalues.method = how; optimvalues.TolFun = deltaF; optimvalues.TolX = deltaX;
                optimvalues.problemtype = type; optimvalues.nonlinineq = Iterate.cineq; optimvalues.nonlineq = Iterate.ceq;

                if (outputTrue)
                    [stopOutput,options,optchanged] = psoutput(OutputFcns,OutputFcnArgs, ...
                        optimvalues,options,state);
                end
                if(plotTrue)
                    stopPlot = psplot(PlotFcns,PlotFcnArgs,PlotInterval,optimvalues,state);
                end
            case 'iter'
                optimvalues.x = X; optimvalues.iteration = OuterIter; optimvalues.fval = Iterate.f;
                optimvalues.meshsize = currentTolMesh;optimvalues.funccount = FunEval;
                optimvalues.method = how; optimvalues.TolFun = deltaF; optimvalues.TolX = deltaX;
                optimvalues.problemtype = type; optimvalues.nonlinineq = Iterate.cineq; optimvalues.nonlineq = Iterate.ceq;
                if (outputTrue)
                    [stopOutput,options,optchanged] = psoutput(OutputFcns,OutputFcnArgs, ...
                        optimvalues,options,state);
                    if optchanged
                        [verbosity,MeshExpansion,MeshContraction,Completesearch,MeshAccelerator,minMesh,MaxMeshSize, ...
                            maxIter,maxFun,TolCon,TolBind,TolFun,TolX,MeshSize,pollmethod, ...
                            pollorder,Completepoll,outputTrue,OutputFcns,OutputFcnArgs,plotTrue,PlotFcns, ...
                            PlotFcnArgs,PlotInterval,searchtype,searchFcnArg,Cache,Vectorized, ...
                            NotVectorizedPoll,NotVectorizedSearch,cachetol,cachelimit,scaleMesh, ...
                            RotatePattern,TimeLimit] =  checkoptions(options,defaultopt,numberOfVariables);
                    end
                end
                if(plotTrue)
                    stopPlot = psplot(PlotFcns,PlotFcnArgs,PlotInterval,optimvalues,state);
                end
            case 'done'
                optimvalues.x = X; optimvalues.iteration = OuterIter; optimvalues.fval = Iterate.f;
                optimvalues.meshsize = currentTolMesh;optimvalues.funccount = FunEval;
                optimvalues.method = how; optimvalues.TolFun = deltaF; optimvalues.TolX = deltaX;
                optimvalues.problemtype = type; optimvalues.nonlinineq = Iterate.cineq; optimvalues.nonlineq = Iterate.ceq;
                if (outputTrue)
                    psoutput(OutputFcns,OutputFcnArgs,optimvalues,options,state);
                end
                if(plotTrue)
                    psplot(PlotFcns,PlotFcnArgs,PlotInterval,optimvalues,state);
                end
        end
    end % End of callOutputPlotFunctions
%------------------------------------------------------------------
end  % End of pfmincon

