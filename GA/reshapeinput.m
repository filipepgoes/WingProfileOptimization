function [X] = reshapeinput(Xin, X)
% RESHAPEINPUT reshape X to match the shape of Xin 

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.2 $  $Date: 2005/05/31 16:30:20 $


[m,n] = size(Xin);
% scalar input
if m == n && m == 1 
    return; % Retain shape
end

% Single point evaluation
if isvector(X) % X is a vector so use shape information
    Xin(:) = X;
    X = Xin;
    return;
elseif isvector(Xin)
    if  m > n && n == 1  % col major input
        return;
    elseif n > m && m == 1 % row major input
        X = X';
    end
else % X is a matrix; not a documented feature
    p = size(X,2);
    if p > 1          % Matrix Xin with vectorized 'on' 
        X = reshape(X,m,n,p);
    else
        X = reshape(X,m,n);
    end
end
