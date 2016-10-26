function  [msg,nextIterate,FunEval,Successdir] = poll(ObjFunc,nonlcon,Xin,pollmethod,Completepoll,pollorder ...
    ,Iterate,Successdir,MeshSize,scale,ToL,A,LB,UB,IndIneqcstr,IndEqcstr,Iter,FunEval,maxFun,type, ...
    NotVectorized,Cache,cachetol,cachelimit,conFcnArg,objFcnArg)
%POLL: Performs the poll step in GPS. This step gaurantees the local convergence.
%
%  OBJFUNC: The objective function on which POLL step is implemented.
%
%  POLLMETHOD: Search directions in POLL step is obtained according to
%  different pollmethod.
%
%  COMPLETEPOLL: If 'off' indicates that POLL can be called off as soon
%  as a better point is found i.e. no sufficient decrease condion is imposed;
%  If 'on' then ALL the points are evaluated and point with least function value
%  is returned. Default is 'off'. If function is expensive, make this 'off'.
%
%  POLLORDER: Ordering of poll directions.
%
%  ITERATE: Incumbent point around which polling will be done. Iterate Stores
%  the current point 'x' and function value 'f' at this point.
%
%  SUCCESSDIR: Last successful POLL/SEARCH direction. This information can be used
%  by the POLL step in ordering the search direction (last successful
%  direction is polled first).
%
%  MESHSIZE: Current mesh size used in POLL step.
%
%  SCALE: Scale factor used to scale the design points.
%
%  TOL: Tolerance used for determining whether constraints are active or not.
%
%  A,LB,UB: Defines the feasible region in case of linear/bound constraints
%  as LB<=A*X<=UB.
%
%  IndIneqcstr: Logical indices of inequality constraints. A(IndIneqcstr), LB(IndIneqcstr)
%  UB(IndIneqcstr) represents inequality constraints.
%
%  IndEqcstr: Logical indices of equality constraints. A(IndEqcstr), LB(IndEqcstr)
%  UB(IndEqcstr) represents equality constraints.
%
%  FUNEVAL: Counter for number of function evaluations. FunEval is always less than
%  'MAXFUN', which is maximum number of function evaluation.
%
%  MAXFUN: Limit on number of function evaluation.
%
%  FLAG: This flag is passed to the SEARCH routines, indicating that the
%  problem is 'unconstrained', 'boundconstraints', 'linearconstraints'.
%
%  NotVectorized: A flag indicating ObjFunc is not evaluated as vectorized
%
%  CACHE: A flag for using CACHE. If 'off', no cache is used.
%
%  CACHETOL: Tolerance used in cache in order to determine whether two points
%  are same or not.
%
%  CACHELIMIT: Limit the cache size to 'cachelimit'.
%
%  NEXTITERATE: Successful iterate after polling is done. If POLL is NOT
%  successful, NEXTITERATE is same as ITERATE.
%
%  MSG: A binary flag indicating, POLL is successful or not.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.15.4.3 $  $Date: 2005/06/21 19:21:48 $
%   Rakesh Kumar

nextIterate = Iterate;
AdaptiveMesh = any(strcmpi(pollmethod,{'madspositivebasisnp1','madspositivebasis2n'}));
%Get the directons which forms the positive basis(minimal or maximal basis)
if  strcmpi(type, 'unconstrained')
    Dirvector = uncondirections(pollmethod,AdaptiveMesh,MeshSize,Iterate.x);
elseif strcmpi(type,'boundconstraints')    % Only Box constraints
    Dirvector = boxdirections(pollmethod,AdaptiveMesh,MeshSize,Iterate.x,A,LB,UB,ToL);
    %If the point is on the constraint boundary, we use scale = 1
    if any(activecstr(Iterate.x,A,LB,UB,ToL))
        scale = 1;
    end
elseif any(strcmpi(type, {'nonlinearconstr','linearconstraints'}))  % Constraints
    Dirvector = lcondirections(pollmethod,AdaptiveMesh,MeshSize,Iterate.x, ...
        A,LB,UB,ToL,IndEqcstr,IndIneqcstr);
    %If the point is on the constraint boundary, we use scale = 1
    if any(activecstr(Iterate.x,A,LB,UB,ToL))
        scale = 1;
    end
end

%make sure that diretcions does not have any empty entries (safegaurd)
Dirvector(:,(~any(Dirvector))) = [];
span = size(Dirvector,2);
%Get order of direction vector based on polling order choice. (Mark Abramson)
% MADS algorithms are not affected by poll order
if AdaptiveMesh
    OrderVec = 1:span;
elseif strcmpi(pollorder,'consecutive')
    OrderVec = 1:span;
elseif strcmpi(pollorder, 'success')
    Successdir = min(Successdir,span);
    OrderVec = 1:span;
    OrderVec = [Successdir OrderVec(1:Successdir-1) OrderVec(Successdir+1:end)];
elseif strcmpi(pollorder, 'random')
    OrderVec = randperm(span);
else
    msg = 'Invalid choice of polling order\nSee PSOPTIMSET for valid choice choosing consecutive order';
    warning('gads:POLL:invalidPollOrder',msg);
    OrderVec = 1:span;
end

%Get the trial points along the direction vectors;
sites = struct('x',cell(span,1,1),'f',cell(span,1,1));
for k = 1:span
    sites(k).x = Iterate.x + MeshSize*scale.*Dirvector(:,OrderVec(k));
end
%Find an iterate with lower objective
[msg,nextIterate,direction,FunEval] = nextfeasible(ObjFunc,nonlcon,Xin,sites,Iterate,A,LB,UB,FunEval, ...
    maxFun,Completepoll,ToL,IndEqcstr,NotVectorized,Cache,cachetol,cachelimit,conFcnArg,objFcnArg{:});
if (msg)
    Successdir = direction;
    if  AdaptiveMesh
        Iterate = nextIterate;
        % Search again along the successful direction
        site.x = Iterate.x + min(1,4*MeshSize)*scale.*Dirvector(:,direction);
        % Find an iterate with lower objective
        [msg1,nextIterate,direction1,FunEval] = nextfeasible(ObjFunc,nonlcon,Xin,site,Iterate,A,LB,UB,FunEval, ...
            maxFun,Completepoll,ToL,IndEqcstr,NotVectorized,Cache,cachetol,cachelimit,conFcnArg,objFcnArg{:});
    end
else
    Successdir = 1;
end
