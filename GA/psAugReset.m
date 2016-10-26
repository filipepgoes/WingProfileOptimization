function [Iter,innerMaxIter,Successdir,deltaF,deltaX,MeshSize, ...
    EXITFLAG,run,step] = psAugReset(options,numberOfVariables,step)
%PSAUGRESET Reset some variables before solving a new sub-problem
%   Private to PATTERNSEARCH.

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/05/31 16:30:11 $

% Reset values before inner iteration
how = '';
Iter = 0;
Successdir = 1;
deltaF = NaN;
deltaX = NaN;
MeshSize = psoptimget(options,'InitialMeshSize',psoptimset,'fast');
EXITFLAG = -1;
if isempty(step)
    run = true;
else
    run = false;
    step = '';
end
% Setting maxiter to Inf is okay because maxfuneval is still finite
innerMaxIter = Inf;%200*numberOfVariables;
% Reset the cache by calling funevaluate
funevaluate('reset');