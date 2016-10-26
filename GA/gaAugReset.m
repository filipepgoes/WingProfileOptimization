function  [exitFlag,innerState,options,step] = gaAugReset(Iterate, ...
    state,options,step,currentTolFun)
%GAAUGRESET Reset some variables before solving a new sub-problem
%   Private to GA.

% Reset values before inner iteration

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/05/31 16:29:42 $

% Reset values before inner iteration
if isempty(step)
    exitFlag = '';
else
    exitFlag = 'Infeasible problem';
    step = '';
end

% Population is created for sub-problem
innerState.Population = repmat(Iterate.x',options.PopulationSize,1);
innerState.Score = repmat(Iterate.f,options.PopulationSize,1);
innerState.FunEval = 0;
% a variety of data used in various places
innerState.Generation = 0;		% current generation counter
innerState.StartTime  = state.StartTime;
innerState.StopFlag   = state.StopFlag;
innerState.LastImprovement = 1;	% generation stall counter
innerState.LastImprovementTime = cputime;	% time stall counter
innerState.Selection = [];       % selection indices
innerState.Expectation = [];     % expection of individuals
innerState.Best = [];            % best score in every generation
options.TolFun = currentTolFun;  % Use the current value of tolerance

