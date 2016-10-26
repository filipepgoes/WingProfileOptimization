function flag = isequalpoint(X,CachePoint,cachetol)
%ISEQUALPOINT: Checks if two points are same within a tolerance
% 	Element wise comparison. Using AbsTol for speed only. For RelTol, choose a
% 	basis; %basis = min(abs(X(:)),abs(Xc(:))); %basis(abs(basis)<eps) = cachetol.

%   Copyright 2003-2004 The MathWorks, Inc.
%   $Revision: 1.6.4.1 $  $Date: 2004/08/20 19:49:16 $

flag = true;
for i = 1:length(X)
    if (abs( X(i) - CachePoint(i) ))  > cachetol
        flag = false;
        return;
    end
end

