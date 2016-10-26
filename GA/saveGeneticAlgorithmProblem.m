function saveGeneticAlgorithmProblem(problem)
%private to gatool, gaguiimportproblem

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.2.4.3 $  $Date: 2005/06/21 19:21:54 $

%Create a temporary structure to save in appdata
if validrandstates(problem)
    tempstruct.randstate = problem.randstate;
    tempstruct.randnstate = problem.randnstate;
else
    tempstruct.randstate = [];
    tempstruct.randnstate = [];
end

if isfield(problem,'fitnessfcn') 
    tempstruct.fitnessfcn = problem.fitnessfcn;
else
    tempstruct.fitnessfcn = [];
end
if isfield(problem,'nvars')
    tempstruct.nvars = problem.nvars;
else
    tempstruct.nvars = [];
end
if isfield(problem,'Aineq')
    tempstruct.Aineq = problem.Aineq;    
else
    tempstruct.Aineq = [];    
end
if isfield(problem,'Bineq')
    tempstruct.Bineq = problem.Bineq;        
else
    tempstruct.Bineq = [];        
end
if isfield(problem,'Aeq')
    tempstruct.Aeq = problem.Aeq;   
else
    tempstruct.Aeq = [];    
end
if isfield(problem,'Beq')
    tempstruct.Beq = problem.Beq;    
else
    tempstruct.Beq = [];
end
if isfield(problem,'LB')
    tempstruct.LB = problem.LB;    
else
    tempstruct.LB = [];    
end
if isfield(problem,'UB')
    tempstruct.UB = problem.UB;    
else
    tempstruct.UB = [];    
end
if isfield(problem,'nonlcon')
    tempstruct.nonlcon = problem.nonlcon;    
else
    tempstruct.nonlcon = [];    
end
setappdata(0,'gads_gatool_problem_data',tempstruct);
%Save options;
if isfield(problem, 'options')
    setappdata(0,'gads_gatool_options_data',problem.options);
end
%------------------------------------------------------------------------
function valid = validrandstates(problem)
    valid = false;
    if isfield(problem, 'randstate') && isfield(problem, 'randnstate') && ...
       isa(problem.randstate, 'double') && isequal(size(problem.randstate),[35, 1]) && ...
       isa(problem.randnstate, 'double') && isequal(size(problem.randnstate),[2, 1])
        valid = true;
    end
