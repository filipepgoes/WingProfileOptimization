function [X,FVAL,EXITFLAG,OUTPUT] = pfminunc(FUN,initialX,options)
%PFMINUNC Finds a box constrained minimum of a function.
%   PFMINUNC solves problems of the form:
%      min F(X)
%       X
%
%   X = PFMINUNC(FUN,X0) starts at X0 and finds a local minimum X to the
%   function FUN. FUN accepts input X and returns a scalar function value F evaluated at X.
%   X0 may be a scalar or vector.
%
%   X = PFMINUNC(FUN,X0,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument created
%   with the PSOPTIMSET function.
%
%   [X,FVAL] = PFMINUNC(FUN,X0,...) returns the value of the objective function
%   FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = PFMINUNC(FUN,X0,...) returns a string EXITFLAG that
%   describes the exit condition of PFMINUNC.
%   If EXITFLAG is:
%        > 0 then PFMINUNC converged to a solution X.
%        = 0 then the algorithm reached the maximum number of iterations or maximum number of function evaluations.
%        < 0 then PFMINUNC did not converge to a solution.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = PFMINUNC(FUN,X0,...) returns a structure
%   OUTPUT with the following information:
%          function: <Objective function>
%       problemtype: <Type of problem> (Unconstrained, Bound constrained or linear constrained)
%        pollmethod: <Polling technique>
%      searchmethod: <Search technique> used, if any
%        iterations: <Total iterations>
%         funcCount: <Total function evaluations>
%          meshsize: <Mesh size at X>
%


%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.19.6.4 $  $Date: 2005/05/31 16:30:06 $

%ALL LOWER CASE ONLY
defaultopt = psoptimset;

%If FUN is a cell array with additional arguments, handle them
if iscell(FUN)
    objFcnArg = FUN(2:end);
    FUN = FUN{1};
else
    objFcnArg = {};
end

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(FUN,'defaults')
    X = defaultopt;
    return
end

%Only function_handle or inlines are allowed
if isempty(FUN) ||  ~(isa(FUN,'inline') || isa(FUN,'function_handle'))
    error('gads:PFMINUNC:needHandleOrInline','Objective function must be a function handle.');
end

%Initialize output args
X = []; FVAL = []; EXITFLAG = [];
OUTPUT = struct('function',FUN,'iterations',0,'funccount',0,'message',' ');

if(~isempty(initialX))
    Iterate.x = initialX(:);
    X    = initialX;
    numberOfVariables = length(Iterate.x);
    type = 'unconstrained';
else
    error('gads:PFMINUNC:initialPoint','You must provide an initial point.');
end
%Retrieve all the options
[verbosity,MeshExpansion,MeshContraction,Completesearch, MeshAccelerator,minMesh,MaxMeshSize, ...
    maxIter,maxFun,TolCon,TolBind,TolFun,TolX,MeshSize,pollmethod, ...
    pollorder,Completepoll,outputTrue,OutputFcns,OutputFcnArgs,plotTrue,PlotFcns,PlotFcnArgs, ...
    PlotInterval,searchtype ,searchFcnArg,Cache,Vectorized,NotVectorizedPoll,NotVectorizedSearch, ...
    cachetol,cachelimit,scaleMesh,RotatePattern,TimeLimit]  = ...
    checkoptions(options,defaultopt,numberOfVariables);

%Get some initial values
[FUN,Iterate,Iter,FunEval,scale,Successdir,nextIterate,deltaF,deltaX,MeshCont,NewMeshSize, ...
    infMessage,how,stopOutput,stopPlot,run,OUTPUT,EXITFLAG,X,FVAL,StartTime] = ...
    getinitial(FUN,X,Iterate,Vectorized,objFcnArg,type,0,MeshContraction,MeshSize, ...
    scaleMesh,numberOfVariables,-Inf*ones(numberOfVariables,1),Inf*ones(numberOfVariables,1));
%Set up output function plot
if outputTrue || plotTrue
    callOutputPlotFunctions('init')
end
%Print some more diagnostic information if verbosity > 2
if verbosity > 2
    psdiagnose(FUN,[],Iterate,initialX,type,[],[],[],options);
end
%Setup display header
if  verbosity>1
    fprintf('\n\nIter     f-count          f(x)      MeshSize     Method\n');
end

while run
    %Check for convergence
    [X,EXITFLAG,FVAL,msg,run] = isconverged(stopOutput,stopPlot,verbosity,Iter, ...
        maxIter,FunEval,maxFun,pollmethod,MeshSize,minMesh,infMessage,nextIterate, ...
        how,deltaX,deltaF,[],TolFun,TolX,X,EXITFLAG,FVAL,run,StartTime,TimeLimit);
    if ~run, continue; end
    %SEARCH.
    [successSearch,nextIterate,FunEval] = search(FUN,[],X,searchtype,Completesearch,Iterate, ...
        Successdir,pollorder,MeshSize,scale,[],[],[],[],[],[],Iter,FunEval,maxFun, ...
        type,NotVectorizedSearch,Cache,cachetol,cachelimit,searchFcnArg,{},objFcnArg);
    %POLL.
    if ~successSearch  %Unsuccessful search
        [successPoll,nextIterate,FunEval,Successdir] = poll(FUN,[],X,pollmethod, ...
            Completepoll,pollorder,Iterate,Successdir,MeshSize,scale,[],[],[],[], ...
            [],[],Iter,FunEval,maxFun,type,NotVectorizedPoll,Cache,cachetol, ...
            cachelimit,{},objFcnArg);
    else
        successPoll =0;
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
end % End of while
% Call output and plot functions
if outputTrue || plotTrue
    callOutputPlotFunctions('done')
end

OUTPUT = struct('function',FUN,'problemtype',type,'pollmethod',pollmethod, ...
    'searchmethod',searchtype,'iterations',Iter,'funccount',FunEval, ...
    'meshsize',MeshSize,'message',msg);

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
end  % End of pfminunc


