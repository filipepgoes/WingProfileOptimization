function [verbosity,MeshExpansion,MeshContraction,Completesearch,MeshAccelerator,minMesh,MaxMeshSize, ...
        maxIter,maxFun,TolCon,TolBind,TolFun, TolX,MeshSize,pollmethod,pollorder, ...
        Completepoll,outputTrue,outputFcns,outputFcnsArg,plotTrue,plotFcns,plotFcnsArg,plotInterval, ...
        searchtype,searchFcnArg,Cache,Vectorized,NotVectorizedPoll,NotVectorizedSearch,cachetol, ...
        cachelimit,scaleMesh,RotatePattern,TimeLimit]  = ...
                        checkoptions(options,defaultopt,numberOfVariables)
%CHECKOPTIONS is private to pfminlcon, pfminbnd and pfminunc.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.10.4.2 $  $Date: 2005/05/31 16:29:35 $

% Define verbosity here (Later we can use options structure)
% Sanity check for the options structure
options = psoptimset(options);

switch  psoptimget(options,'Display',defaultopt,'fast')
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

% Retrieve options using PSOPTIMGET
MeshExpansion     = psoptimget(options,'MeshExpansion',defaultopt,'fast'); 
MeshContraction   = psoptimget(options,'MeshContraction',defaultopt,'fast'); 
Completesearch    = psoptimget(options,'CompleteSearch',defaultopt,'fast');
MeshAccelerator   = psoptimget(options,'MeshAccelerator',defaultopt,'fast');
minMesh           = psoptimget(options,'TolMesh',defaultopt,'fast');
TolCon            = psoptimget(options,'TolCon',defaultopt,'fast');
MaxMeshSize       = psoptimget(options,'MaxMeshSize',defaultopt,'fast');
maxIter           = psoptimget(options,'MaxIter',defaultopt,'fast');
maxFun            = psoptimget(options,'MaxFunEvals',defaultopt,'fast');
TimeLimit         = psoptimget(options,'TimeLimit',defaultopt,'fast');
TolBind           = psoptimget(options,'TolBind',defaultopt,'fast');
TolFun            = psoptimget(options,'TolFun',defaultopt,'fast');
TolX              = psoptimget(options,'TolX',defaultopt,'fast');
MeshSize          = psoptimget(options,'InitialMeshSize',defaultopt,'fast');
pollmethod        = psoptimget(options,'PollMethod',defaultopt,'fast');
pollorder         = psoptimget(options,'PollingOrder',defaultopt,'fast');
Completepoll      = psoptimget(options,'CompletePoll',defaultopt,'fast');
outputFcns        = psoptimget(options,'OutputFcns',defaultopt,'fast');
plotFcns          = psoptimget(options,'PlotFcns',defaultopt,'fast');
plotInterval      = psoptimget(options,'PlotInterval',defaultopt,'fast');
searchtype        = psoptimget(options,'SearchMethod',defaultopt,'fast');
Vectorized        = psoptimget(options,'Vectorized',defaultopt,'fast');
Cache             = psoptimget(options,'Cache',defaultopt,'fast');
cachetol          = psoptimget(options,'CacheTol',defaultopt,'fast');
cachelimit        = psoptimget(options,'CacheSize',defaultopt,'fast');
scaleMesh         = psoptimget(options,'ScaleMesh',defaultopt,'fast');
RotatePattern     = psoptimget(options,'MeshRotate',defaultopt,'fast');

% Modify some fields if they are not yet assigned
if ischar(maxFun)
    if isequal(lower(maxFun),'2000*numberofvariables')
        maxFun = 2000*numberOfVariables;
    else
        error('gads:CHECKOPTIONS:maxFunEvals','Option ''MaxFunEvals'' must be a positive numeric if not the default.')
    end
end
if ischar(maxIter)
    if isequal(lower(maxIter),'100*numberofvariables')
        maxIter = 100*numberOfVariables;
    else
        error('gads:CHECKOPTIONS:maxIter','Option ''MaxIter'' must be a positive numeric if not the default.')
    end
end

maxFun  = floor(maxFun);
maxIter = floor(maxIter);

% If searchtype is a cell array with additional arguments, handle them
if iscell(searchtype)
    searchFcnArg = searchtype(2:end);
    searchtype = searchtype{1};
else
    searchFcnArg = {};
    searchtype = searchtype;
end
% Search technique could be any poll methods, customSearch or None
if isempty(searchtype)
    searchtype = [];
elseif isa(searchtype,'function_handle')
    [searchtype,msg] = fcnchk(searchtype);
    if ~isempty(msg)
        error('gads:CHECKOPTIONS:searchFcnArgCheck',msg);
    end
elseif ischar(searchtype) && ~any(strcmpi(searchtype,{'positivebasisnp1', 'positivebasis2n', ...
        'gpspositivebasisnp1', 'gpspositivebasis2n','madspositivebasisnp1', 'madspositivebasis2n'}))  
    [searchtype,msg] = fcnchk(searchtype);
    if ~isempty(msg)
        error('gads:CHECKOPTIONS:searchFcnArgCheck',msg);
    end
elseif ischar(searchtype) && any(strcmpi(searchtype,{'positivebasisnp1', 'positivebasis2n', ...
        'gpspositivebasisnp1', 'gpspositivebasis2n','madspositivebasisnp1', 'madspositivebasis2n'}))  
    [searchtype] = fcnchk(searchtype);
else
    error('gads:CHECKOPTIONS:invalidSearchType','Invalid choice of Search technique: See psoptimset for SearchMethod.\n');
end

% Make sure that search method is different from poll method
if ~isempty(searchtype) && strcmpi(func2str(searchtype),pollmethod)
   searchtype = [];
end
% If MaxMeshSize is less than Meshsize (This should not happen)
if MaxMeshSize <MeshSize
    if verbosity > 0
        warning('gads:CHECKOPTIONS:MaxMeshSize','MaxMeshSize should be greater than InitialMeshSize.\n');
    end
    MeshSize = MaxMeshSize;
end

% It is NOT vectorized in these conditions
NotVectorizedPoll   = (strcmpi(Vectorized,'off') || (strcmpi(Vectorized, 'on') && strcmpi(Completepoll,'off')));
NotVectorizedSearch = (strcmpi(Vectorized,'off') || (strcmpi(Vectorized, 'on') && strcmpi(Completesearch,'off')));

% If using 2N basis or MADS RotatePattern has no effect.
if any(strcmpi(pollmethod,{'positivebasis2n','gpspositivebasis2n','madspositivebasisnp1','madspositivebasis2n'}))
    RotatePattern = 'off';
end

% If outputFcns is a cell array with additional arguments, handle them
[outputFcns,outputFcnsArg] = functionHandleOrCellArray(outputFcns,'OutputFcns');

 if isempty(outputFcns)
     outputTrue = false;
 else
     outputTrue = true;
 end
 
 %If plotFcns is a cell array with additional arguments, handle them
[plotFcns,plotFcnsArg] = functionHandleOrCellArray(plotFcns,'PlotFcns');

 if isempty(plotFcns)
     plotTrue = false;
 else
     plotTrue = true;
 end
 
%-----------------------------------------------------------------------

function [FUN,ARGS] = functionHandleOrCellArray(value,property)

FUN  = [];
ARGS = [];

%if a scalar  ~cell  is passed convert to cell (for clarity, not speed)
if ~iscell(value) && length(value) == 1
     value = {value};
end

% If value is an array of functions, it must be a cell array
for i = 1:length(value)
    candidate = value(i);
    %If any element is also a cell array
    if iscell(candidate)
        if isempty(candidate{1})
            continue;
        end
        %Sometimes the variable 'candidate' might have nested cell array 
        %e.g. {{@outputfcn, p1,p2}} instead of just
        %{@outputfcn,p1,p2}. The following code gets rid of extra braces,
        %which are typically introduced by GUI import/export options.
        temp = candidate{1};
        while iscell(temp) && isscalar(temp)
            candidate = temp(1);
            temp = candidate{1};
        end
        [handle,args] = isFcn(candidate{:});
    else
        [handle,args] = isFcn(candidate);
    end
    if(~isempty(handle)) && (isa(handle,'inline') || isa(handle,'function_handle'))
        FUN{i} = handle;
        ARGS{i} = args;
    else
        msg = sprintf('The field ''%s'' must contain a function handle.',property);
        error('gads:CHECKOPTIONS:FUNCTIONHANDLEORCELLARRAY:needHandleOrInline',msg);
    end
end

%-------------------------------------------------------------------------
% if it's a scalar fcn handle or a cellarray starting with a fcn handle and
% followed by something other than a fcn handle, return parts, else empty
function [handle,args] =  isFcn(x)
  handle = [];
  args = {};
  %If x is a cell array with additional arguments, handle them
  if iscell(x) 
      if ~isempty(x)
          args = x(2:end);
          handle = x{1};
      else  %Cell could be empty too
          args = {};
      end
  else % Not a cell
      args = {};
      handle = x;
  end
  
  if ~isempty(handle)
      [handle,msg] = fcnchk(handle);
      if ~isempty(msg)
          handle =[];
      end
  end

  
