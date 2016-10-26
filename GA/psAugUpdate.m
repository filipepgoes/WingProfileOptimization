function  [lambda,lambdabar,penalty,currentOmega,currentTolMesh,currentTolCon,shift,OuterIter,how] = ...
    psAugUpdate(Iterate,lambda,penalty,startOmega,startTolMesh,startTolCon,currentOmega,currentTolMesh, ...
    currentTolCon,minMesh,TolCon,shift,penaltyFactor,betamesh,betaconstr,alphamesh,alphaconstr,alphaL, ...
    numNonlinIneqcstr,numNonlinCstr,OuterIter)
%PSAUGUPDATE Updates values of parameters for augmented Lagragian PS 
%   Private to PATTERNSEARCH

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2.2.1 $  $Date: 2005/07/24 20:56:57 $
%    Rakesh Kumar

how = ' ';
OuterIter = OuterIter + 1;
stallTol = min(minMesh,eps);
% Initialize new values of lambda
lambdabar = zeros(numNonlinCstr,1);
C1 = 0;
C2 = 0;
% Calculate switching condition which will decide whether to
% update multipliers or increase penalty. We also update lambda.

if numNonlinIneqcstr % Nonlinear inequality constraints
    shiftedConstr = -Iterate.cineq + shift(1:numNonlinIneqcstr);
    lambdabar(1:numNonlinIneqcstr) = (lambda(1:numNonlinIneqcstr).*shift(1:numNonlinIneqcstr))./shiftedConstr;
    if numNonlinCstr > numNonlinIneqcstr % Eualities too
        lambdabar(numNonlinIneqcstr+1:numNonlinCstr) = ...
            lambda(numNonlinIneqcstr+1:numNonlinCstr) + Iterate.ceq.*penalty;
    end
    for i = 1:numNonlinIneqcstr
        if lambda(i) > eps
          C1 =  C1 + (Iterate.cineq.*lambdabar(i))./(lambda(i).^alphaL);
        end
    end
    C1 = norm(C1);
    C2 = norm(Iterate.ceq);
else % Only nonlinear equality constraints
    lambdabar(numNonlinIneqcstr+1:numNonlinCstr) = ...
        lambda(numNonlinIneqcstr+1:numNonlinCstr) + Iterate.ceq.*penalty;
    C2 = norm(Iterate.ceq);
end
% Update multipliers?
updateLang = max(C1 , C2) <= currentTolCon;
if updateLang % Update Lagrange multipliers estimate
    how = 'Update multipliers';
    invPenalty = max(stallTol,1/penalty);
    % Update these three quantities.
    lambda = lambdabar;
    currentOmega = min(1e-1,currentOmega*(invPenalty^betamesh));
    currentTolMesh = currentOmega;
    currentTolCon  = min(1e-1,currentTolCon*(invPenalty^betaconstr));
else % Increase penalty
    how = 'Increase penalty';
    penalty = penaltyFactor*penalty;
    invPenalty = max(stallTol,1/penalty);
    currentOmega = min(1e-1,startOmega*(invPenalty^alphamesh));
    currentTolMesh = currentOmega; % Lewis & Torczon
    currentTolCon = min(1e-1,startTolCon*(invPenalty^alphaconstr));
end
% If currentTolMesh, currentTolCon is too low then use penalty values to
% limit accuracy of inner sub-problem's solution and constraint tolerance
if currentTolMesh <= minMesh
    currentTolMesh = max(currentOmega,currentTolMesh);
end
if currentTolCon <= TolCon
    currentTolCon = max(TolCon,currentTolCon);
end

% Compute shift
shift = invPenalty*(lambda.^alphaL);
