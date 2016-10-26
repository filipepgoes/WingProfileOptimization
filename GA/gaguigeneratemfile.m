function gaguigeneratemfile(probmodel, optmodel, randchoice)
%GAGUIGENERATEMFILE generates an mfile from gatool. 

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.6.4.2.2.1 $ $Date: 2005/07/17 06:06:33 $

[fitnessFcn,nvars,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,randstate,randnstate] =  gaguiReadProblem(probmodel);

options = gaguiReadHashTable(optmodel);

%remove special gui outputfcn which is the first in the list
if ~isempty(options.OutputFcns) 
    temp = options.OutputFcns{1};
    temp = temp{1};
    if strcmp(func2str(temp), 'gatooloutput')
        options.OutputFcns(1) = [];
    end
end
%Create a struct for generateMfile
tempstruct = struct;
tempstruct.fitnessfcn = fitnessFcn;
tempstruct.nvars = nvars;
tempstruct.LB = LB;
tempstruct.UB = UB;
tempstruct.Aineq = Aineq;
tempstruct.Bineq = Bineq;
tempstruct.Aeq = Aeq;
tempstruct.Beq = Beq;
tempstruct.nonlcon = nonlcon;
if randchoice
    tempstruct.randstate = randstate;
    tempstruct.randnstate = randnstate;
end
tempstruct.options=options;
%Call generate Mfile code
generateMfile(tempstruct, 'ga');
