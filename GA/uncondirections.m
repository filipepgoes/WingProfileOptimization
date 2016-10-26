function directions = uncondirections(pollmethod,AdaptiveMesh,MeshSize,x)
%UNCONDIRECTIONS: finds search vectors when no constraints are present.
%   POLLMETHOD: Poll method used to get search vectors.
%
% 	X: Point at which polling is done (usually the best point found so
% 	far).
% 	
% 	DIRECTIONS:  Returns direction vectors that positively span the tangent
% 	cone at the current iterate, with respect to bound and linear constraints.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.8.4.3 $  $Date: 2005/05/31 16:30:25 $
%   Rakesh Kumar
vars = length(x);
% N linearly independent vectors
if ~AdaptiveMesh
    Basis  = eye(vars);
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

% Form directions that forms the positive basis
switch lower(pollmethod)
    case {'positivebasisnp1','gpspositivebasisnp1','madspositivebasisnp1'} % Minimal positive basis (n+1 vectors)
        directions = [-1*ones(vars,1) Basis];
case {'positivebasis2n','gpspositivebasis2n','madspositivebasis2n'}   % Maximal positive basis (2n vectors)
        directions = [Basis -Basis];
    otherwise
        error('gads:UNCONDIRECTIONS:pollmethod','Invalid choice of Poll method.');
end
