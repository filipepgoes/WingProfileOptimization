function gadiagnose(FUN,nonlcon,GenomeLength,type,nineqcstr,neqcstr,ncstr,options)
%GADIAGNOSE prints some diagnostic information about the problem
%   private to GA

%   Copyright 2004-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/05/31 16:29:45 $

properties =  optionsList('ga');
defaultOpt = gaoptimset;
Output_String = sprintf('\nDiagnostic information.');

Output_String = [Output_String sprintf('\n\tFitness function = %s',value2RHS(FUN))];
if ~isempty(GenomeLength)
    Output_String = [Output_String sprintf('\n\tNumber of variables = %d',GenomeLength)];
end

%print some information about constraints
if ~isempty(nonlcon)
    Output_String = [Output_String sprintf('\n\tnonlinear constraint function = %s',value2RHS(nonlcon))];
end
if ~isempty(nineqcstr)
    Output_String = [Output_String sprintf('\n\t%d Inequality constraints',nineqcstr)];
end
if ~isempty(neqcstr)
    Output_String = [Output_String sprintf('\n\t%d Equality constraints',neqcstr)];
end
if ~isempty(ncstr)
    Output_String = [Output_String sprintf('\n\t%d Total number of linear constraints\n',ncstr)];
end

Output_String = [Output_String sprintf('\n%s','Modified options:')];
for i = 1:length(properties)
    prop = properties{i};
    if(~isempty(prop)) % the property list has blank lines, ignore them
        value = options.(prop);
        if ~(isequal(value,defaultOpt.(prop)) || isempty(value)) 
            Output_String = [Output_String sprintf('\n\toptions.%s = %s',prop,value2RHS(value))];
        end
    end
end
Output_String = [Output_String sprintf('\nEnd of diagnostic information.')];
fprintf('%s',Output_String)
