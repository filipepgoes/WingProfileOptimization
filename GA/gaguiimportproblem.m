function [selection, problemModel, optionModel] = gaguiimportproblem()
%GAGUIIMPORT GUI helper  

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.14.4.3 $  $Date: 2005/06/21 19:21:38 $


selection = '';
optionModel = '';
problemModel = '';
names = {};
whoslist =evalin('base','whos');

for i = 1:length(whoslist)
    if strcmp(whoslist(i).class, 'struct') && strcmp(num2str(whoslist(i).size), '1  1')
        s = evalin('base', whoslist(i).name);
        if isfield(s, 'fitnessfcn') && isfield(s, 'nvars') ...
            && isfield(s, 'options') && validOptions(s.options)
            names{end + 1 } = whoslist(i).name;
        end
    end
end
   
 
if isempty(names) 
    msgbox('There are no problem structures in the workspace.', 'Genetic Algorithm Tool');
else
    [value, OK] = listdlg('ListString', names, 'SelectionMode', 'Single', ...
            'ListSize', [250 200], 'Name', 'Import GA Problem', ...
            'PromptString', 'Select a problem structure to import:', ...
            'OKString', 'Import');
    if OK == 1
        selection = names{value};
        %stuff all the fields into the hashtable.
        problemModel = java.util.Hashtable;
        problem = evalin('base', selection);
        problemModel.put('fitnessfcn', value2RHS(problem.fitnessfcn));
        problemModel.put('nvars', value2RHS(problem.nvars));
        if validrandstates(problem)
            problemModel.put('randchoice', true);
        else
            problemModel.put('randchoice', false);
        end
        % Add constraint related fields if present
        if isfield(problem, 'Aineq')
            problemModel.put('Aineq', value2RHS(problem.Aineq));
        else
            aineq = '';
        end
        
        if isfield(problem, 'Bineq')
            problemModel.put('Bineq', value2RHS(problem.Bineq));
        else
            bineq = '';
        end
        if isfield(problem, 'Aeq')
            problemModel.put('Aeq', value2RHS(problem.Aeq));
        else
            aeq = '';
        end
        if isfield(problem, 'Beq')
            problemModel.put('Beq', value2RHS(problem.Beq));
        else
            beq = '';
        end
        if isfield(problem, 'LB')
            problemModel.put('LB', value2RHS(problem.LB));
        else
            lb = '';
        end
        if isfield(problem, 'UB')
            problemModel.put('UB', value2RHS(problem.UB));
        else
            ub = '';
        end
        if isfield(problem, 'nonlcon')
            problemModel.put('nonlcon', value2RHS(problem.nonlcon));
        else
            nonlcon = '';
        end
        options = problem.options;
        s = struct(options);
        f = fieldnames(s);
        optionModel = java.util.Hashtable;
        gafieldnames = fieldnames(gaoptimset);
        for i = 1:length(f);
            n = f{i};
            if ismember(n, gafieldnames)
                rhs = value2RHS(s.(n));
                % remove string quotes
                q = find(rhs == '''');
                rhs(q) = [];
                optionModel.put(n,rhs);
            end
        end
        saveGeneticAlgorithmProblem(problem);
     end
end    

 %--------------------------------------------------------------------------
 function valid = validOptions(options)
    valid = false;
    gafieldnames = fieldnames(gaoptimset);
    ofieldnames = fieldnames(options);
    if nnz(ismember(gafieldnames, ofieldnames)) >= 25
        valid = true;
        return;
    end
 
 %--------------------------------------------------------------------------
 function valid = validrandstates(problem)
    valid = false;
    if isfield(problem, 'randstate') && isfield(problem, 'randnstate') && ...
       isa(problem.randstate, 'double') && isequal(size(problem.randstate),[35, 1]) && ...
       isa(problem.randnstate, 'double') && isequal(size(problem.randnstate),[2, 1])
        valid = true;
    end
     
 
    
    
