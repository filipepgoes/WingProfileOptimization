function state = makeState(GenomeLength,FitnessFcn,options)
%MAKESTATE Create an initial population and fitness scores
%   state = makeStates(GenomeLength,FitnessFcn,options) Creates an initial 
%   state structure using the information in the options structure.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.6.4.2 $  $Date: 2005/05/31 16:30:00 $

if isempty(options.InitialPopulation) % we have to make the initial pop.
    state.Population = feval(options.CreationFcn,GenomeLength,FitnessFcn,options,options.CreationFcnArgs{:});
else % the initial pop was passed in!
    state.Population = options.InitialPopulation;
end

if(isempty(options.InitialScores))
    % score each member of the population
    if strcmpi(options.Vectorized, 'off')
        try
            state.Score = feval(@fcnvectorizer,state.Population,FitnessFcn,options.FitnessFcnArgs{:});
        catch
            error('gads:MAKESTATE:fitnessCheck', ...
                'GA cannot continue because user supplied fitness function failed with the following error:\n%s', lasterr)
        end
        if numel(state.Score) ~=size(state.Population,1)
            msg = sprintf('%s\n', ...
                'Your fitness function must return a scalar value.');
            error('gads:MAKESTATE:fitnessCheck',msg);
        end
    else
        try
            state.Score = feval(FitnessFcn,state.Population,options.FitnessFcnArgs{:});
        catch
            error('gads:MAKESTATE:fitnessCheck', ...
                'GA cannot continue because user supplied fitness function failed with the following error:\n%s', lasterr)
        end
        if numel(state.Score) ~=size(state.Population,1)
            msg = sprintf('%s\n', ...
                ['When ''Vectorized'' is ''on'', your fitness function must ' ...
                'return a vector of length equal to the size of the population.']);
            error('gads:MAKESTATE:fitnessCheck',msg);
        end
    end
    state.FunEval = length(state.Score);          % number of function evaluations
else
    state.Score = options.InitialScores;
    if numel(state.Score) ~=size(state.Population,1)
            msg = sprintf('%s\n', ...
                'Length of InitialScores must be equal to the size of the population.');
            error('gads:MAKESTATE:scoreCheck',msg);
    end
    state.FunEval = 0;          % number of function evaluations
end

% Partial population is allowed for 'doubleVector' population type so make
% population of appropriate size
if strcmpi(options.PopulationType,'doubleVector')
    lens = size(state.Population,1);
    npop = sum(options.PopulationSize);
    if npop > lens
        population = zeros(npop,GenomeLength);
        population(1:lens,:) = state.Population;
        population(lens+1:end,:) = repmat(state.Population(end,:),(npop-lens),1);
        scores = zeros(npop,1);
        scores(1:lens) = state.Score;
        scores(lens+1:end) = repmat(state.Score(end),(npop-lens),1);
        state.Population = population;
        state.Score = scores;
    else
        state.Population(npop+1:end,:) = [];
        state.Score(npop+1:end) = [];
    end
end

% a variety of data used in various places
state.Generation = 0;		% current generation counter
state.StartTime = cputime;	% start time
state.StopFlag = []; 		% reason for termination
state.LastImprovement = 1;	% generation stall counter
state.LastImprovementTime = state.StartTime;	% time stall counter
state.Selection = [];       % selection indices
state.Expectation = [];     % expection of individuals
state.Best = [];            % best score in every generation
