function [X,EXITFLAG,FVAL,msg,run] = isconverged(stopOutput,stopPlot,verbosity,Iter,maxIter,FunEval,maxFun, ...
    pollmethod,MeshSize,minMesh,infMessage,nextIterate,how,deltaX,deltaF,TolCon,TolFun, ...
    TolX,X,EXITFLAG,FVAL,run,StartTime,TimeLimit)
%ISCONVERGED Checks several conditions of convergence.
%
% 	STOP: A flag passed by user to stop the iteration (Used from OutPutFcn)
%
% 	VERBOSITY: Level of display
%
% 	ITER, MAXITER: Current Iteration and maximum iteraion allowed respectivley
%
% 	FUNEVAL,MAXFUN: Number of function evaluation and maximum iteraion
% 	allowed respectively
%
% 	MESHSIZE,MINMESH; Current mesh size used and minimum mesh size
% 	allowed respectively
%
% 	INFMESSAGE: Has function reached -inf value (if yes then stop)
%
% 	NEXTITERATE: Next iterate is stored in this structure nextIterate.x
%   and nextIterate.f

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.15.6.5 $  $Date: 2005/06/21 19:21:43 $
%   Rakesh Kumar

X(:) = nextIterate.x;
FVAL = nextIterate.f;
msg = '';
if verbosity > 1
    fprintf('%5.0f    %5.0f   %12.6g  %12.4g     %s\n',Iter, FunEval, nextIterate.f, MeshSize, how);
end

% User interruption
if stopOutput || stopPlot
    msg = sprintf('%s','Stop requested.');
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    EXITFLAG = -1;
    run = false;
    return;
end
if ~isempty(infMessage)
    msg = sprintf('%s','Optimization terminated: ');
    msg = [msg, sprintf('%s','objective function has reached -Inf value (objective function is unbounded below).')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    EXITFLAG = 1;
    run = false;
    return;
end

% Convergence check is different for adaptive mesh and fixed mesh
% algorithms
AdaptiveMesh = any(strcmpi(pollmethod,{'madspositivebasisnp1','madspositivebasis2n'}));

% Check mesh size parameter for fixed mesh direct search
if MeshSize < minMesh && (deltaF < TolFun || deltaX < TolX) && ...
        ~AdaptiveMesh
    EXITFLAG = 1;
    run  = false;
    msg = sprintf('%s','Optimization terminated: ');
    msg = [msg,sprintf('%s', 'mesh size less than options.TolMesh.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% X and Fun tolerances will be used only when iteration is successful and
% Meshsize is of the order of TolX/TolFun.
if ~strcmpi(how,'Mesh refined') && ~AdaptiveMesh && ...
        ((MeshSize < TolX || MeshSize < TolFun)&& (deltaF < TolFun || deltaX < TolX))
    EXITFLAG = 1;
    run  = false;
    msg = sprintf('%s','Optimization terminated: ');
    if deltaX < TolX && MeshSize < TolX
        msg = [msg, sprintf('%s', 'change in X less than options.TolX.')];
    else
        msg = [msg, sprintf('%s', 'change in the function value less than options.TolFun.')];
    end
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end
% Check poll size parameter for mesh adaptive direct search
if AdaptiveMesh
    if strcmpi(pollmethod,'madspositivebasisnp1')
        framesize = length(X)*sqrt(MeshSize);
    else
        framesize = sqrt(MeshSize);
    end
    if framesize < minMesh && (deltaF < TolFun || deltaX < TolX)
        EXITFLAG = 1;
        run  = false;
        msg = sprintf('%s','Optimization terminated: ');
        if strcmpi(pollmethod,'madspositivebasisnp1')
            msg = [msg,sprintf('%s', 'mesh size less than ''numberOfVariables*sqrt(TolMesh)''.')];
        else
            msg = [msg,sprintf('%s', 'mesh size less than ''sqrt(TolMesh)''.')];
        end
        if verbosity > 0
            fprintf('%s\n',msg);
        end
        return;
    end

    % X and Fun tolerances will be used only when iteration is successful and
    % Meshsize is of the order of TolX/TolFun.
    if ~strcmpi(how,'Mesh refined') && ((framesize < TolX || framesize < TolFun) ...
            && (deltaF < TolFun || deltaX < TolX))
        EXITFLAG = 1;
        run  = false;
        msg = sprintf('%s','Optimization terminated: ');
        if deltaX < TolX && framesize < TolX
            msg = [msg, sprintf('%s', 'change in X less than options.TolX.')];
        else
            msg = [msg, sprintf('%s', 'change in the function value less than options.TolFun.')];
        end
        if verbosity > 0
            fprintf('%s\n',msg);
        end
        return;
    end

end
% Maximum iteration limit
if Iter > maxIter
    EXITFLAG = 0;
    run  = false;
    msg = sprintf('%s', 'Maximum number of iterations exceeded: ');
    msg = [msg,sprintf('%s', 'increase options.MaxIter.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end
% Maximum function evaluation limit
if FunEval >= maxFun
    EXITFLAG = 0;
    run  = false;
    msg = sprintf('%s', 'Maximum number of function evaluations exceeded: ');
    msg = [msg, sprintf('%s', 'increase options.MaxFunEvals.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end
% Maximum time limit
if(cputime-StartTime) > TimeLimit
    EXITFLAG = 0;
    run  = false;
    msg = sprintf('%s', 'Time limit reached: ');
    msg = [msg, sprintf('%s', 'increase options.TimeLimit.')];
    if verbosity > 0
        fprintf('%s\n',msg);
    end
    return;
end

% Setup display header every thirty iterations
if verbosity > 1 && rem(Iter,30)== 0 && Iter >0 && Iter < maxIter
    fprintf('\nIter     f-count        f(x)       MeshSize      Method\n');
end
