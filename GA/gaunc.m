function [x,fval,exitFlag,output,population,scores] = gaunc(FitnessFcn,GenomeLength,options)
%GAUNC Genetic algorithm unconstrained solver.
%   X = GAUNC(FITNESSFCN,NVARS) finds the minimum of FITNESSFCN using
%   GAUNC. NVARS is the dimension (number of design variables) of the
%   FITNESSFCN. FITNESSFCN accepts a vector X of size 1-by-NAVRS,
%   and returns a scalar evaluated at X. 
% 		
%   X = GAUNC(FITNESSFCN,NAVRS,OPTIONS) finds the minimum for
%   FITNESSFCN with the default optimization parameters replaced by values
%   in the structure OPTIONS. OPTIONS can be created with the GAOPTIMSET
%   function.
% 		
%   [X, FVAL] = GAUNC(FITNESSFCN, ...) returns FVAL, the value of the fitness
%   function FITNESSFCN at the solution X.
% 		
%   [X,FVAL,REASON] = GAUNC(FITNESSFCN, ...) returns the REASON for stopping.
%
%   [X,FVAL,REASON,OUTPUT] = GAUNC(FITNESSFCN, ...) returns a
%   structure OUTPUT with the following information: 
%            randstate: <State of the function RAND used before GAUNC started>
%           randnstate: <State of the function RANDN used before GAUNC started>
%          generations: <Total generations, excluding HybridFcn iterations>
%            funccount: <Total function evaluations>
%              message: <GAUNC termination message>
%
%   [X,FVAL,REASON,OUTPUT,POPULATION] = GAUNC(FITNESSFCN, ...) returns the final
%   POPULATION at termination.
% 		
%   [X,FVAL,REASON,OUTPUT,POPULATION,SCORES] = GAUNC(FITNESSFCN, ...) returns the
%   SCORES of the final POPULATION.
% 		
%   See also GAOPTIMSET, FITNESSFUNCTION, PATTERNSEARCH, @.

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2.2.1 $  $Date: 2005/07/17 06:06:35 $

% If the first arg is not a gaoptimset, then it's a fitness function followed by a genome
% length. Here we make a gaoptimset from the args.
defaultopt = gaoptimset;   

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(FitnessFcn,'defaults') 
    x = defaultopt;
    return
end
% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
    error('gads:GAUNC:tenthInputNotStruct','Invalid input argument to GAUNC.');
elseif isempty(options)
    options = gaoptimset;
end

x =[];fval =[];exitFlag='';population=[];scores=[];user_options = options;
type = 'unconstrained';
%Remember the random number states used
output.randstate  = rand('state');
output.randnstate = randn('state');
output.generations = 0;
output.funccount   = 0;
output.message   = '';   
        
%Validate arguments
[GenomeLength,FitnessFcn,options] = validate(GenomeLength,FitnessFcn,options);
%Create initial state: population, scores, status data
state = makeState(GenomeLength,FitnessFcn,options);
%Give the plot/output Fcns a chance to do any initialization they need.
state = gaplot(FitnessFcn,options,state,'init');
[state,options] = gaoutput(FitnessFcn,options,state,'init');

%Print some diagnostic information if asked for
if strcmpi(options.Display,'diagnose')
   gadiagnose(FitnessFcn,[],GenomeLength,type,[],[],[],user_options); 
end
%Setup display header 
if  any(strcmpi(options.Display, {'iter','diagnose'}))
    fprintf('\n                               Best           Mean      Stall\n');
    fprintf('Generation      f-count        f(x)           f(x)    Generations\n');
end

exitFlag = '';
% run the main loop until some termination condition becomes true
while isempty(exitFlag)
        state.Generation = state.Generation + 1;
        %Repeat for each subpopulation (element of the populationSize vector)
        offset = 0;
        totalPop = options.PopulationSize;
        % each sub-population loop
        for pop = 1:length(totalPop)
            populationSize =  totalPop(pop);
            thisPopulation = 1 + (offset:(offset + populationSize - 1));
            population = state.Population(thisPopulation,:);
            score = state.Score( thisPopulation );
            %Empty population is also possible
            if isempty(thisPopulation)
                continue;
            end
            [score,population,state] = stepGA(score,population,options,state,GenomeLength,FitnessFcn);
            
            % store the results for this sub-population
            state.Population(thisPopulation,:) = population;
            state.Score(thisPopulation) = score;
            offset = offset + populationSize;
        end 
        
        % remember the best score
        scores = state.Score;
        best = min(state.Score);
        generation = state.Generation;
        state.Best(generation) = best;
        
        % keep track of improvement in the best
        if((generation > 1) && finite(best))
            if(state.Best(generation-1) > best)
                state.LastImprovement = generation;
                state.LastImprovementTime = cputime;
            end
        end
        
        % do any migration
        state = migrate(FitnessFcn,GenomeLength,options,state);
        % update the Output
        state = gaplot(FitnessFcn,options,state,'iter');
        [state,options] = gaoutput(FitnessFcn,options,state,'iter');

        % check to see if any stopping criteria have been met
        exitFlag = isItTimeToStop(options,state);
   end %End while loop

% find and return the best solution
[fval,best] = min(state.Score);
x = state.Population(best,:);

%Update output structure
output.generations = state.Generation;
output.message     = exitFlag;
output.funccount   = state.Generation*length(state.Score);

% load up outputs
if(nargout > 4)
    population = state.Population;
    if(nargout > 5)
        scores = state.Score;
    end
end

%A hybrid scheme. Try another minimization method if there is one.
if(strcmpi(options.PopulationType,'doubleVector') && ~isempty(options.HybridFcn))
    %Who is the hybrid function
    if isa(options.HybridFcn,'function_handle')
        hfunc = func2str(options.HybridFcn);
    else 
        hfunc = options.HybridFcn;
    end
    % Create functions handle to be passed to hybrid function
    FitnessHybridFcn = @(x) FitnessFcn(x,options.FitnessFcnArgs{:});
    %Inform about hybrid scheme
    if  any(strcmpi(options.Display, {'iter','diagnose'}))
        fprintf('%s%s%s\n','Switching to the hybrid optimization algorithm (',upper(hfunc),').');
    end

    %Determine which syntax to call
    switch hfunc
        case 'fminsearch'
            [xx,ff,e,o] = feval(options.HybridFcn,FitnessHybridFcn,x,options.HybridFcnArgs{:});
            output.funccount = output.funccount + o.funcCount;
            output.message   = [output.message sprintf('\nFMINSEARCH:\n'), o.message];
        case 'patternsearch'
            [xx,ff,e,o] = feval(options.HybridFcn,FitnessHybridFcn,x,[],[],[],[],[],[],[],options.HybridFcnArgs{:});
            output.funccount = output.funccount + o.funccount;
            output.message   = [output.message sprintf('\nPATTERNSEARCH: \n'), o.message];
        case 'fminunc'
            [xx,ff,e,o] = feval(options.HybridFcn,FitnessHybridFcn,x,options.HybridFcnArgs{:});
            output.funccount = output.funccount + o.funcCount;
            output.message   = [output.message sprintf('\nFMINUNC: \n'), o.message];
        case {'fmincon'}
             msg = sprintf('%s is a constrained optimization solver',upper(hfunc));
             msg = [msg, sprintf('\n%s',' using unconstrained solver FMINUNC as hybrid function.')];
             warning('gads:GAUNC:constrainedHybridFcn',msg);
             % We need to store this warning state so that we can
             % display it in the GUI
             [lastmsg, lastid] = lastwarn; warning off;
            [xx,ff,e,o] = feval(@fminunc,FitnessHybridFcn,x,options.HybridFcnArgs{:});
            warning on; lastwarn(lastmsg,lastid);
            output.funccount = output.funccount + o.funcCount;
            output.message   = [output.message sprintf('\nFMINUNC: \n'), o.message];
        otherwise
            error('gads:GAUNC:hybridFcnError','Hybrid function must be one of the following:\n@FMINSEARCH, @FMINUNC, @PATTERNSEARCH.')
    end
    if e > 0 && ff < fval
        fval = ff;
        x = xx;
    end
    %Inform about hybrid scheme termination
    if  any(strcmpi(options.Display, {'iter','diagnose'}))
        fprintf('%s%s\n',upper(hfunc), ' terminated.');
    end
end

% give the Output functions a chance to finish up
gaplot(FitnessFcn,options,state,'done');
gaoutput(FitnessFcn,options,state,'done');
