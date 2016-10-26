function [lowerbounds,upperbounds] = checkconstraints(x,A,LB,UB,ToL)
%CHECKCONSTRINTS determines the active lower and upper consraints with
% 	respect to A, LB and UB with a a specified tolerance 'tol'
% 	
% 	LOWERBOUNDS, UPPERBOUNDS are indices of active constraints w.r.t. lower
% 	and upper bounds (LB and UB)

%   Copyright 2003-2004 The MathWorks, Inc.
%   $Revision: 1.4.4.1 $  $Date: 2004/08/20 19:48:53 $
%   Rakesh Kumar

%setup the costraint status A*x; we already have LB and UB.
Ax = A*x;
%Check the tolerance with respect to each constraints;
lowerbounds = (abs(LB-Ax)<=ToL);
upperbounds = (abs(Ax-UB)<=ToL);

