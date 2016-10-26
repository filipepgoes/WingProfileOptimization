function p = propertyList
%These are the properties that display, and generate operate on.
%   p = propertyList; returns a cell array of strings that list all of the
%   user visibnle properties in the gaoptimset struct. There are empty
%   strings in the list to seperate functional areas from one another.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.7.4.3 $  $Date: 2005/05/31 16:30:08 $

p = {
    'CreationFcn'
    'PopulationType'
    'PopInitRange'
    'InitialPenalty'
    'PenaltyFactor'
    ''
    'OutputFcns'
    'PlotFcns'
    'Display'
    ''
    'FitnessScalingFcn'
    'SelectionFcn'
    'CrossoverFcn'
    'PopulationSize'
    'EliteCount'
    'CrossoverFraction'
    ''
    'MutationFcn'
    ''
    'MigrationInterval'
    'MigrationFraction'
    'MigrationDirection'
    ''
    'TolFun'
    'TolCon'
    'Generations'
    'TimeLimit'
    'FitnessLimit'
    'StallGenLimit'
    'StallTimeLimit'
    ''
    'HybridFcn'
};
