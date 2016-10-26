function directions = lcondirections(pollmethod,AdaptiveMesh,MeshSize,x,A,LB,UB,tol,IndEqcstr,IndIneqcstr)
%LCONDIRECTIONS finds search vectors when linear constraints and bounds are present.
% 	POLLMETHOD: Poll method used to get search vectors.
%
% 	X: Point at which polling is done (usually the best point found so far)
%
% 	A,LB,UB: Defines the feasible region in case of linear/bound constraints as L<=A*X<=U.
%
% 	TOL: Tolerance used for determining whether constraints are active or not.
%
% 	IndIneqcstr: Logical indices of inequality constraints. A(IndIneqcstr), LB(IndIneqcstr)
% 	UB(IndIneqcstr) represents inequality constraints.
%
% 	IndEqcstr: Logical indices of equality constraints. A(IndEqcstr), LB(IndEqcstr)
% 	UB(IndEqcstr) represents equality constraints.
%
% 	DIRECTIONS:  Returns direction vectors that positively span the tangent
% 	cone at the current iterate, with respect to bound and linear constraints.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.11.4.4 $  $Date: 2005/06/21 19:21:44 $
%   Rakesh Kumar


%Initialization
LB(IndEqcstr) = -inf;
TangentCone = [];
vars = length(x);
% N linearly independent vectors
if ~AdaptiveMesh
    Basis = eye(vars);
else
    % Create random generating directions
    pollParam = 1/sqrt(MeshSize);
    lowerT = tril((round((pollParam+1)*rand(vars)-0.5)),-1);
    diagtemp = pollParam*sign(rand(vars,1)-0.5);
    diagtemp(~diagtemp) = pollParam*sign(0.5-rand);
    diagT  = diag(diagtemp);
    Basis = lowerT + diagT;
    order = randperm(vars);
    Basis = Basis(order,order);
end
% Separate constraints and bounds 
Acstr = [A(IndEqcstr,:); A(IndIneqcstr,:)];
Lcstr = [LB(IndEqcstr); LB(IndIneqcstr)];
Ucstr = [UB(IndEqcstr); UB(IndIneqcstr)];
Abnd  = A((end-vars)+1:end,:);
Lbnd  = LB((end-vars)+1:end);
Ubnd  = UB((end-vars)+1:end);

Normals = zeros(vars,1);
tolDep = 100*vars*eps;
Tol = tol;
% The cone generators for minimumm epsilon is in active set
% (Lewis & Torczon section 8.2)
if ~isempty(Acstr)
    while rank(Normals) ~= min(size(Normals))
        if tol < tolDep
            error('gads:LCONDIRECTIONS:degenconstr','Constraints are dependent at current iterate\nTry increasing OPTIONS.TolBind (<eps).');
        end
        [lowerbounds,upperbounds] = checkconstraints(x,Acstr,Lcstr,Ucstr,tol);
        Normals = [Acstr(upperbounds,:); -Acstr(lowerbounds,:)]';
        tol = tol/2;
    end

    % Lewis & Torczon section 8.2. T = V*inv(V'V), which is computed using QR
    % decomposition
    if (~isempty(Normals))
        [Q,R] = qr(Normals,0);
        TangentCone = Q/R';
        Basis = Basis - TangentCone*Normals';
    end
end

% Add active bounds in search directions
I = eye(vars);
% Check which bounds are active
[lowerbounds,upperbounds] = checkconstraints(x,Abnd,Lbnd,Ubnd,Tol);
active  = lowerbounds | upperbounds; 
% Include all directions parallel to active bounds
TangentCone = [TangentCone, I(:,active)];

% Form directions that positively span the tangent cone at x
switch lower(pollmethod)
    case {'positivebasisnp1','gpspositivebasisnp1','madspositivebasisnp1'}
        directions = [-sum(Basis,2) Basis  TangentCone -TangentCone];
    case {'positivebasis2n','gpspositivebasis2n','madspositivebasis2n'}
        directions = [Basis -Basis TangentCone -TangentCone];
    otherwise
        error('gads:LCONDIRECTIONS:pollmethod','Invalid choice of Poll method.');
end
