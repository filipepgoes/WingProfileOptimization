function directions = boxdirections(pollmethod,AdaptiveMesh,MeshSize,x,A,LB,UB,tol)
%BOXDIRECTIONS finds search vectors when bound constraints are present.
%   POLLMETHOD: Poll method used to get search vectors.
% 	
% 	X: Point at which polling is done (usually the best point found so
% 	far).
% 	
% 	A,LB,UB: Defines the feasible region in case of linear constraints.
% 	L<=A*X<=U.
% 	
% 	TOL: Tolerance used for determining whether constraints are active or not.
% 	
% 	DIRECTIONS:  Returns direction vectors that positively span the tangent
% 	cone at the current iterate, with respect to bound and linear constraints.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.8.4.3 $  $Date: 2005/05/31 16:29:34 $
%   Rakesh Kumar

vars = length(x);
I = eye(vars);
%Check which constraints are active for LB <= AX <= UB at 'x'
[lowerbounds,upperbounds] = checkconstraints(x,A,LB,UB,tol);
active  = lowerbounds | upperbounds; 
%Include all directions parallel to active constraints
TangentCone = I(:,active);

% N linearly independent vectors 
if ~AdaptiveMesh
   Basis = I(:,~active);
else
    pollParam = 1/sqrt(MeshSize);
    lowerT = tril((round((pollParam+1)*rand(vars)-0.5)),-1);
    diagtemp = pollParam*sign(rand(vars,1)-0.5);
    diagtemp(~diagtemp) = pollParam*sign(0.5-rand);
    diagT  = diag(diagtemp);
    Basis = lowerT + diagT;
    order = randperm(vars);
    Basis = Basis(order,order);
end

% Form directions that positively span the tangent cone at x
switch lower(pollmethod)
    case {'positivebasisnp1','gpspositivebasisnp1','madspositivebasisnp1'}
        directions = [-sum(Basis,2) Basis  TangentCone -TangentCone];
    case {'positivebasis2n','gpspositivebasis2n','madspositivebasis2n'}
        directions = [Basis -Basis TangentCone -TangentCone];
    otherwise
        error('gads:BOXDIRECTIONS:pollmethod','Invalid choice of Poll method.');
end

