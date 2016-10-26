function [x,fval,exitFlag,output,population,scores] = ga(FUN,GenomeLength,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options)
%GA    finds a constrained minimum of a function of several variables. 
%   GA attempts to solve problems of the form: 
%       min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%        X                     C(X) <= 0, Ceq(X) = 0 (nonlinear constraints)
%                              LB <= X <= UB
%           
%   X = GA(FITNESSFCN,NVARS) finds a local unconstrained minimum X to the
%   FITNESSFCN using GA. NVARS is the dimension (number of design
%   variables) of the FITNESSFCN. FITNESSFCN accepts a vector X of size
%   1-by-NVARS, and returns a scalar evaluated at X.
%
%   X = GA(FITNESSFCN,NVARS,A,b) finds a local minimum X to the function
%   FITNESSFCN, subject to the linear inequalities A*X <= B. 
%
%   X = GA(FITNESSFCN,NVARS,A,b,Aeq,beq) finds a local minimum X to the
%   function FITNESSFCN, subject to the linear equalities Aeq*X = Beq as
%   well as A*X <= B. (Set A=[] and B=[] if no inequalities exist.)
%
%   X = GA(FITNESSFCN,NVARS,A,b,Aeq,beq,LB,UB) defines a set of lower and
%   upper bounds on the design variables, X, so that a solution is found in
%   the range LB <= X <= UB. Use empty matrices for LB and UB if no bounds
%   exist. Set LB(i) = -Inf if X(i) is unbounded below;  set UB(i) = Inf if
%   X(i) is unbounded above.
%
%   X = GA(FITNESSFCN,NVARS,A,b,Aeq,beq,LB,UB,NONLCON) subjects the
%   minimization to the constraints defined in NONLCON. The function
%   NONLCON accepts X and returns the vectors C and Ceq, representing the
%   nonlinear inequalities and equalities respectively. GA minimizes
%   FUN such that C(X)<=0 and Ceq(X)=0. (Set LB=[] and/or UB=[] if no
%   bounds exist.) 
%
%   X = GA(FITNESSFCN,NVARS,A,b,Aeq,beq,LB,UB,NONLCON,options) minimizes
%   with the default optimization parameters replaced by values in the
%   structure OPTIONS. OPTIONS can be created with the GAOPTIMSET function.
%   See GAOPTIMSET for details.
%
%   X = GA(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a structure
%   that has the following fields:
%       fitnessfcn: <Fitness Function>
%            nvars: <Number of design variables>
%          options: <Options structure created with GAOPTIMSET>
%            Aineq: <A matrix for inequality constraints>
%            Bineq: <B vector for inequality constraints>
%              Aeq: <A matrix for equality constraints>
%              Beq: <B vector for equality constraints>
%               LB: <Lower bound on X>
%               UB: <Upper bound on X>
%          nonlcon: <nonlinear constraint function>
%        randstate: <Optional field to reset rand state>
%       randnstate: <Optional field to reset randn state>
%
%   [X,FVAL] = GA(FITNESSFCN, ...) returns FVAL, the value of the fitness
%   function FITNESSFCN at the solution X.
%
%   [X,FVAL,REASON] = GA(FITNESSFCN, ...) returns the REASON for stopping.
%
%   [X,FVAL,REASON,OUTPUT] = GA(FITNESSFCN, ...) returns a
%   structure OUTPUT with the following information:
%            randstate: <State of the function RAND used before GA started>
%           randnstate: <State of the function RANDN used before GA started>
%          generations: <Total generations, excluding HybridFcn iterations>
%            funccount: <Total function evaluations>
%        maxconstraint: <Maximum constraint violation>, if any
%              message: <GA termination message>
%
%   [X,FVAL,REASON,OUTPUT,POPULATION] = GA(FITNESSFCN, ...) returns the
%   final POPULATION at termination.
%
%   [X,FVAL,REASON,OUTPUT,POPULATION,SCORES] = GA(FITNESSFCN, ...) returns
%   the SCORES of the final POPULATION.
%
%
%   Example:
%     Unconstrained minimization of 'rastriginsfcn' fitness function of
%     numberOfVariables = 2 
%     x = ga(@rastriginsfcn,2)
%
%     Display plotting functions while GA minimizes
% 	      options = gaoptimset('PlotFcns',...
% 	      {@gaplotbestf,@gaplotbestindiv,@gaplotexpectation,@gaplotstopping});
% 	      [x,fval,reason,output] = ga(@rastriginsfcn,2,[],[],[],[],[],[],[],options)
%
%   An example with inequality constraints and lower bounds
%    A = [1 1; -1 2; 2 1];  b = [2; 2; 3];  lb = zeros(2,1);
%    % Use mutation function which can handle constraints
%    options = gaoptimset('MutationFcn',@mutationadaptfeasible);
%    [X,FVAL,EXITFLAG] = ga(@lincontest6,2,A,b,[],[],lb,[],[],options);
%
%     FUN can also be an anonymous function:
%        X = ga(@(x) 3*sin(x(1))+exp(x(2)),2)
%
%   If FUN or NONLCON are parameterized, you can use anonymous functions to
%   capture the problem-dependent parameters. Suppose you want to minimize
%   the fitness given in the function myfit, subject to the nonlinear
%   constraint myconstr, where these two functions are parameterized by their
%   second argument a1 and a2, respectively. Here myfit and myconstr are
%   M-file functions such as 
%
%        function f = myfit(x,a1)
%        f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + a1);
%
%   and
%
%        function [c,ceq] = myconstr(x,a2)
%        c = [1.5 + x(1)*x(2) - x(1) - x(2); 
%              -x(1)*x(2) - a2];
%        % No nonlinear equality constraints:
%         ceq = [];
%
%   To optimize for specific values of a1 and a2, first assign the values
%   to these two parameters. Then create two one-argument anonymous
%   functions that capture the values of a1 and a2, and call myfit and
%   myconstr with two arguments. Finally, pass these anonymous functions to
%   GA: 
%
%     a1 = 1; a2 = 10; % define parameters first
%     % Mutation function for constrained minimization
%     options = gaoptimset('MutationFcn',@mutationadaptfeasible);
%     x = ga(@(x)myfit(x,a1),2,[],[],[],[],[],[],@(x)myconstr(x,a2),options)
%
%   See also GAOPTIMSET, FITNESSFUNCTION, GAOUTPUTFCNTEMPLATE,
%   PATTERNSEARCH, @.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.28.4.4.2.1 $  $Date: 2005/07/17 06:06:29 $

% If the first arg is not a gaoptimset, then it's a fitness function followed by a genome
% length. Here we make a gaoptimset from the args.
defaultopt = struct('PopulationType', 'doubleVector', ...
    'PopInitRange', [0;1], ...
    'PopulationSize', 20, ...
    'EliteCount', 2, ...
    'CrossoverFraction', 0.8, ...
    'MigrationDirection','forward', ...
    'MigrationInterval',20, ...
    'MigrationFraction',0.2, ...
    'Generations', 100, ...
    'TimeLimit', inf, ...
    'FitnessLimit', -inf, ...
    'StallGenLimit', 50, ...
    'StallTimeLimit', 20, ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-6, ...
    'InitialPopulation',[], ...
    'InitialScores', [], ...
    'InitialPenalty', 10, ...
    'PenaltyFactor', 100, ...
    'PlotInterval',1, ...
    'CreationFcn',@gacreationuniform, ...
    'FitnessScalingFcn', @fitscalingrank, ...
    'SelectionFcn', @selectionstochunif, ...
    'CrossoverFcn',@crossoverscattered, ...
    'MutationFcn',@mutationgaussian, ...
    'HybridFcn',[], ...
    'Display', 'final', ...
    'PlotFcns', [], ...
    'OutputFcns', [], ...
    'Vectorized','off');

% Check number of input arguments
errmsg = nargchk(1,10,nargin);
if nargin <1
    error('gads:GA:numberOfInputs',[errmsg,' GA requires at least 1 input argument.']);
elseif ~isempty(errmsg)
    error('gads:GA:numberOfInputs',[errmsg,' GA takes at most 10 input arguments.']);
end
% Check number of output arguments
errmsg = nargoutchk(0,6,nargin);
if nargout > 6 && ~isempty(errmsg)
    error('gads:GA:outputArg',[errmsg,' GA returns at most 6 output arguments.']);
end

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(FUN,'defaults')
    x = defaultopt;
    return
end
if nargin < 10,  options = [];
    if nargin < 9,  nonlcon = [];
        if nargin < 8, UB = [];
            if nargin < 7, LB = [];
                if nargin <6, Beq = [];
                    if nargin <5, Aeq = [];
                        if nargin < 4, Bineq = [];
                            if nargin < 3, Aineq = [];
                            end, end, end, end, end, end, end, end

% Is third argument a structure
if nargin == 3 && isstruct(Aineq) % Old syntax
    options = Aineq; Aineq = [];
    if ~isempty(options) && ~isstruct(options)
        error('gads:GA:thirdInputNotStruct','Invalid Input for GA.');
    end
end
% Input can be a problem structure
if nargin == 1
    try
        options = FUN.options;
        GenomeLength = FUN.nvars;
        % If using new syntax then must have all the fields; check one
        % field
        if isfield(FUN,'Aineq')
            Aineq    = FUN.Aineq;
            Bineq    = FUN.Bineq;
            Aeq      = FUN.Aeq;
            Beq      = FUN.Beq;
            LB       = FUN.LB;
            UB       = FUN.UB;
            nonlcon  = FUN.nonlcon;
        else
            Aineq = []; Bineq = [];
            Aeq = []; Beq = [];
            LB = []; UB = [];
            nonlcon = [];
        end
        % optional fields
        if isfield(FUN, 'randstate') && isfield(FUN, 'randnstate') && ...
                isa(FUN.randstate, 'double') && isequal(size(FUN.randstate),[35, 1]) && ...
                isa(FUN.randnstate, 'double') && isequal(size(FUN.randnstate),[2, 1])
            rand('state',FUN.randstate);
            randn('state',FUN.randnstate);
        end
         FUN = FUN.fitnessfcn;
    catch
        error('gads:GA:invalidStructInput','The input should be a structure with valid fields or provide two arguments to GA.' );
    end
end
% We need to check the GenomeLength here before we call any solver
valid =  isnumeric(GenomeLength) && isscalar(GenomeLength)&& (GenomeLength > 0) ...
         && (GenomeLength == floor(GenomeLength));
if(~valid)
   msg = sprintf('Number of variables (NVARS) must be a positive integer.');
   error('gads:GA:validNumberofVariables:notValidNvars',msg);
end
% All inputs should be double
try
    dataType = superiorfloat(GenomeLength,Aineq,Bineq,Aeq,Beq,LB,UB);
    if ~isequal('double', dataType)
        error('gads:GA:dataType', ...
            'GA only accepts inputs of data type double.')
    end
catch
    error('gads:GA:dataType', ...
        'GA only accepts inputs of data type double.')
end

% Call appropriate solver
% If constraints/bounds then call constraint solver
if ~isempty(nonlcon)
    [x,fval,exitFlag,output,population,scores] = gacon(FUN,GenomeLength,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options);
% If Aeq or Aineq is not empty, then problem has linear constraints. 
elseif ~isempty(Aeq) || ~isempty(Aineq) 
    [x,fval,exitFlag,output,population,scores] = galinconf(FUN,GenomeLength,Aineq,Bineq,Aeq,Beq,LB,UB,options);    
    % This condition satisfies bound constraints 
elseif (isempty(Aeq) && isempty(Aineq) && isempty(Bineq) && isempty(Beq)) && ( ~isempty(LB) || ~isempty(UB) ) 
    [x,fval,exitFlag,output,population,scores] = galinconf(FUN,GenomeLength,Aineq,Bineq,Aeq,Beq,LB,UB,options);
% Unconstrained problem
elseif (isempty(Aeq) && isempty(Aineq) && isempty(Bineq) && isempty(Beq) && isempty(LB) && isempty(UB))
    [x,fval,exitFlag,output,population,scores] = gaunc(FUN,GenomeLength,options);
   % Try nonlinear constrained solver
else  
    try  
        [x,fval,exitFlag,output,population,scores] = gacon(FUN,GenomeLength,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options);  
    catch
        error('gads:GA:entryPoint','%s','Invalid input to GA'); 
    end
end
