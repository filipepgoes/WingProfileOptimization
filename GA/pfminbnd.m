function [X,FVAL,EXITFLAG,OUTPUT] = pfminbnd(FUN,initialX,LB,UB,options)
%PFMINBND Finds minimum of a function with bound constraints.
%   PFMINBND solves problems of the form:
%        min F(X)  subject to: LB <= X <= UB  (Bound constraints) 
%         X                  
%           
%   X = PFMINBND(FUN,X0,LB,UB) starts at X0 and finds a local minimum X to the
%   function FUN, subject to the bounds LB <= X <= UB. FUN accepts input X
%   and returns a scalar function value F evaluated at X. X0 may be a scalar,
%   or vector. If LB (or UB) is empty (or not provided) it is automatically 
%   expanded to -Inf (or Inf)
%
%   X = PFMINBND(FUN,X0,LB,UB,OPTIONS) minimizes with the default optimization 
%   parameters replaced by values in the structure OPTIONS. OPTIONS can be created 
%   with the PSOPTIMSET function.
%
%   [X,FVAL] = PFMINBND(FUN,X0,LB,UB,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = PFMINBND(FUN,X0,LB,UB,...) returns a string EXITFLAG that 
%   describes the exit condition of PFMINBND.  
%     If EXITFLAG is:
%        > 0 then PFMINBND converged to a solution X.
%        = 0 then the algorithm reached the maximum number of iterations or maximum number of function evaluations.
%        < 0 then PFMINBND did not converge to a solution.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = PFMINBND(FUN,X0,...) returns a structure
%   OUTPUT with the following information:
%          function: <Objective function>
%       problemtype: <Type of problem> (Unconstrained, Bound constrained or linear constrained)
%        pollmethod: <Polling technique>
%      searchmethod: <Search technique> used, if any
%        iterations: <Total iterations>
%         funcCount: <Total function evaluations>
%          meshsize: <Mesh size at X>



%   Copyright 2003-2005 The MathWorks, Inc. 
%   $Revision: 1.20.6.5 $  $Date: 2005/06/21 19:21:45 $


defaultopt = psoptimset;  

%If FUN is a cell array with additional arguments, handle them
if  iscell(FUN)
    objFcnArg = FUN(2:end);
    FUN = FUN{1};
else
    objFcnArg = {};
    FUN = FUN;
end

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(FUN,'defaults')
    X = defaultopt;
    return
end

%Only function_handle or inlines are allowed
if isempty(FUN) ||  ~(isa(FUN,'inline') || isa(FUN,'function_handle'))
    error('gads:PFMINBND:needHandleOrInline','Objective function must be a function handle.');
end

%Initialize output args
X = []; FVAL = []; EXITFLAG = [];
OUTPUT = struct('function',FUN,'iterations',0,'funccount',0,'message',' ');


if(~isempty(initialX))
    Iterate.x = initialX(:);
    X = initialX;
    numberOfVariables = length(Iterate.x);
    type = 'boundconstraints';
else
    error('gads:PFMINBND:initialPoint','You must provide an initial point.');
end
% Retrieve all the options
[verbosity,MeshExpansion,MeshContraction,Completesearch, MeshAccelerator,minMesh,MaxMeshSize, ...
 maxIter,maxFun,TolCon,TolBind,TolFun, TolX, MeshSize, pollmethod, ...
 pollorder,Completepoll,outputTrue,OutputFcns,OutputFcnArgs,plotTrue,PlotFcns,PlotFcnArgs, ...
 PlotInterval,searchtype ,searchFcnArg,Cache,Vectorized,NotVectorizedPoll, ...
 NotVectorizedSearch,cachetol,cachelimit,scaleMesh,RotatePattern,TimeLimit]  = ...
              checkoptions(options,defaultopt,numberOfVariables);

% Bound correction
[LB,UB,msg,EXITFLAG] = checkbound(LB,UB,numberOfVariables,verbosity);
if EXITFLAG < 0
    FVAL = [];
    X(:) = Iterate.x;
    OUTPUT.message = msg;
    if verbosity > 0
         fprintf('%s\n',msg)
    end
    return
end
% Make Iterate feasible and form A.
[Iterate,A] =aluform(Iterate.x,[],[],[],[],LB,UB,numberOfVariables,type,verbosity);
%Get some initial values
[FUN,Iterate,Iter,FunEval,scale,Successdir,nextIterate,deltaF,deltaX,MeshCont,NewMeshSize, ...
        infMessage,how,stopOutput,stopPlot,run,OUTPUT,EXITFLAG,X,FVAL,StartTime] = getinitial(FUN,X,Iterate,Vectorized, ...
    objFcnArg,type,0,MeshContraction,MeshSize,scaleMesh,numberOfVariables,LB,UB);

% Call output and plot functions
if outputTrue || plotTrue
    callOutputPlotFunctions('init')
end
%Print some more diagnostic information if verbosity > 2
if verbosity > 2
    psdiagnose(FUN,[],Iterate,initialX,type,[],[],[],options);
end
%Setup display header 
if  verbosity > 1
    fprintf('\n\nIter     f-count          f(x)      MeshSize     Method\n');
end

while run 
    
    %Check for convergence
    [X,EXITFLAG,FVAL,msg,run] = isconverged(stopOutput,stopPlot,verbosity,Iter,maxIter,FunEval,maxFun, ...
        pollmethod,MeshSize,minMesh,infMessage,nextIterate,how,deltaX,deltaF,[], ...
        TolFun,TolX,X,EXITFLAG,FVAL,run,StartTime,TimeLimit);
    if ~run
        continue;
    end
    %SEARCH.
    [successSearch,nextIterate,FunEval] = search(FUN,[],X,searchtype,Completesearch,Iterate, ...
        Successdir,pollorder,MeshSize,scale,TolBind,A,LB,UB,[],[],Iter,FunEval, ...
        maxFun,type,NotVectorizedSearch,Cache,cachetol,cachelimit,searchFcnArg,{},objFcnArg);
    %POLL
    if ~successSearch  %Unsuccessful search
        [successPoll,nextIterate,FunEval,Successdir] = poll(FUN,[],X,pollmethod,Completepoll,pollorder ...
            ,Iterate,Successdir,MeshSize,scale,TolBind,A,LB,UB,[],[],Iter,FunEval,maxFun, ...
            type,NotVectorizedPoll,Cache,cachetol,cachelimit,{},objFcnArg);
    else
        successPoll =0;
    end 
    
    %Scale the variables (if needed)
    if strcmpi(scaleMesh,'on') 
        meanX = mean([Iterate.x],2);
        scale = logscale(LB,UB,meanX);
    end
    
% Call output and plot functions
if outputTrue || plotTrue
    callOutputPlotFunctions('iter')
end
    %UPDATE
    [NewMeshSize,MeshContraction,how,deltaX,deltaF,scale,Iterate,X,Iter,infMessage] = ...
        updateparam(pollmethod,successPoll,successSearch,MeshAccelerator,RotatePattern, ...
        MaxMeshSize,minMesh,MeshExpansion,MeshCont,MeshContraction,MeshSize,scale, ...
        nextIterate,Iterate,X,Iter,how,infMessage);
   
    %Update mesh size 
    MeshSize = NewMeshSize;
end

% Call output and plot functions
if outputTrue || plotTrue
    callOutputPlotFunctions('done')
end
OUTPUT = struct('function',FUN,'problemtype',type,'pollmethod',pollmethod,'searchmethod',searchtype, ...
    'iterations',Iter,'funccount',FunEval,'meshsize',MeshSize,'maxconstraint',0.0,'message',msg);
%-----------------------------------------------------------------
% Nested function to call output/plot functions
    function callOutputPlotFunctions(state)
        switch state
            case 'init'
                optimvalues.x = X; optimvalues.iteration = Iter; optimvalues.fval = Iterate.f;
                optimvalues.problemtype = type; optimvalues.meshsize = MeshSize; optimvalues.funccount = FunEval;
                optimvalues.method = how; optimvalues.TolFun = deltaF; optimvalues.TolX = deltaX;
                if (outputTrue)
                    [stopOutput,options,optchanged] = psoutput(OutputFcns,OutputFcnArgs, ...
                        optimvalues,options,state);
                end
                if(plotTrue)
                    stopPlot = psplot(PlotFcns,PlotFcnArgs,PlotInterval,optimvalues,state);
                end
            case 'iter'
                optimvalues.x = X; optimvalues.iteration = Iter; optimvalues.fval = Iterate.f;
                optimvalues.problemtype = type; optimvalues.meshsize = MeshSize; optimvalues.funccount = FunEval;
                optimvalues.method = how; optimvalues.TolFun = deltaF; optimvalues.TolX = deltaX;
                if (outputTrue)
                    [stopOutput,options,optchanged] = psoutput(OutputFcns,OutputFcnArgs, ...
                        optimvalues,options,state);
                    if optchanged
                        [verbosity,MeshExpansion,MeshContraction,Completesearch, MeshAccelerator,minMesh, ...
                            MaxMeshSize,maxIter,maxFun,TolCon,TolBind,TolFun,TolX, ...
                            InitialMeshSize,pollmethod, pollorder,Completepoll,outputTrue,OutputFcns, ...
                            OutputFcnArgs,plotTrue,PlotFcns, PlotFcnArgs,PlotInterval,searchtype, ...
                            searchFcnArg,Cache,Vectorized,NotVectorizedPoll,NotVectorizedSearch,cachetol, ...
                            cachelimit,scaleMesh,RotatePattern,TimeLimit]  = ...
                            checkoptions(options,defaultopt,numberOfVariables);
                    end
                end
                if(plotTrue)
                    stopPlot = psplot(PlotFcns,PlotFcnArgs,PlotInterval,optimvalues,state);
                end
            case 'done'
                optimvalues.x = X; optimvalues.iteration = Iter; optimvalues.fval = Iterate.f;
                optimvalues.problemtype = type; optimvalues.meshsize = MeshSize; optimvalues.funccount = FunEval;
                optimvalues.method = how; optimvalues.TolFun = deltaF; optimvalues.TolX = deltaX;
                if (outputTrue)
                    psoutput(OutputFcns,OutputFcnArgs,optimvalues,options,state);
                end
                if(plotTrue)
                    psplot(PlotFcns,PlotFcnArgs,PlotInterval,optimvalues,state);
                end
        end
    end % End of callOutputPlotFunctions
%------------------------------------------------------------------
end  % End of pfminbnd
