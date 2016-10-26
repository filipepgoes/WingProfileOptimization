% pso_filipe.m
% a generic particle swarm optimizer
% to find the minimum or maximum of any 
% MISO matlab function
%
% Implements Common, Trelea type 1 and 2, and Clerc's class 1". It will
% also automatically try to track to a changing environment (with varied
% success - BKB 3/18/05)
%
% This vectorized version removes the for loop associated with particle
% number. It also *requires* that the cost function have a single input
% that represents all dimensions of search (i.e., for a function that has 2
% inputs then make a wrapper that passes a matrix of ps x 2 as a single
% variable)
%
% Usage:
%  [optOUT]=PSO(functname,D)
% or:
%  [optOUT,tr,te]=...
%        PSO(functname,D,mv,VarRange,minmax,PSOparams,plotfcn,PSOseedValue)
%
% Inputs:
%    functname - string of matlab function to optimize
%    D - # of inputs to the function (dimension of problem)
%    
% Optional Inputs:
%    mv - max particle velocity, either a scalar or a vector of length D
%           (this allows each component to have it's own max velocity), 
%           default = 4, set if not input or input as NaN
%
%    VarRange - matrix of ranges for each input variable, 
%      default -100 to 100, of form:
%       [ min1 max1 
%         min2 max2
%            ...
%         minD maxD ]
%
%    minmax = 0, funct minimized (default)
%           = 1, funct maximized
%           = 2, funct is targeted to P(12) (minimizes distance to errgoal)
%
%    PSOparams - PSO parameters
%      P(1) - Maximum number of iterations (epochs) to train, default = 2000.
%      P(2) - population size, default = 24
%
%      P(3) - acceleration const 1 (local best influence), default = 2
%      P(4) - acceleration const 2 (global best influence), default = 2
%      P(5) - Initial inertia weight, default = 0.9
%      P(6) - Final inertia weight, default = 0.4
%      P(7) - Epoch when inertial weight at final value, default = 1500
%      P(8)- minimum global error gradient, 
%                 if abs(Gbest(i+1)-Gbest(i)) < gradient over 
%                 certain length of epochs, terminate run, default = 1e-25
%      P(9)- epochs before error gradient criterion terminates run, 
%                 default = 150, if the SSE does not change over 250 epochs
%                               then exit
%      P(10)- error goal, if NaN then unconstrained min or max, default=NaN
%      P(11)- PSOseed, default=0
%               = 0 for initial positions all random
%               = 1 for initial particles as user input
%
%
%    PSOseedValue - initial particle position, depends on P(13), must be
%                   set if P(13) is 1 or 2, not used for P(13)=0, needs to
%                   be nXm where n<=ps, and m<=D
%                   If n<ps and/or m<D then remaining values are set random
%                   on Varrange
% Outputs:
%    optOUT - optimal inputs and associated min/max output of function, of form:
%        [ bestin1
%          bestin2
%            ...
%          bestinD
%          bestOUT ]
%
% Optional Outputs:
%    tr    - Gbest at every iteration, traces flight of swarm
%    te    - epochs to train
%
% Example:  out=pso_Trelea_vectorized('f6',2)

% Brian Birge
% Rev 3.3
% 2/18/06
function [OUT,varargout]=pso_filipe(functname,D,varargin)

rand('state',sum(100*clock));
if nargin < 2
   error('Not enough arguments.');
end

% PSO PARAMETERS
if nargin == 2      % only specified functname and D
   VRmin=ones(D,1)*-100;
   VRmax=ones(D,1)*100 ;   
   VR=[VRmin,VRmax];
   minmax = 0;
   P = [];
   mv = 4;
elseif nargin == 3  % specified functname, D, and mv
   VRmin=ones(D,1)*-100; 
   VRmax=ones(D,1)*100;    
   VR=[VRmin,VRmax];
   minmax = 0;
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end
   P = [];
elseif nargin == 4  % specified functname, D, mv, Varrange
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end
   VR=varargin{2}; 
   minmax = 0;
   P = [];
elseif nargin == 5  % Functname, D, mv, Varrange, and minmax
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = [];
elseif nargin == 6  % Functname, D, mv, Varrange, minmax, and psoparams
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = varargin{4}; % psoparams  
elseif nargin == 7  % Functname, D, mv, Varrange, minmax, and psoparams, PSOseedValue
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = varargin{4}; % psoparams
   PSOseedValue = varargin{5};
else    
   error('Wrong # of input arguments.');
end

% sets up default pso params
Pdef = [2000 24 2 2 0.9 0.4 1500 1e-25 250 NaN 0];
Plen = length(P);
P    = [P,Pdef(Plen+1:end)];

me      = P(1);
ps      = P(2);
ac1     = P(3);
ac2     = P(4);
iw1     = P(5);
iw2     = P(6);
iwe     = P(7);
ergrd   = P(8);
ergrdep = P(9);
errgoal = P(10);
PSOseed = P(11);

% error checking
 if ((minmax==2) & isnan(errgoal))
     error('minmax= 2, errgoal= NaN: choose an error goal or set minmax to 0 or 1');
 end

 if ( (PSOseed==1) & ~exist('PSOseedValue') )
     error('PSOseed flag set but no PSOseedValue was input');
 end

 if exist('PSOseedValue')
     tmpsz=size(PSOseedValue);
     if D < tmpsz(2)
         error('PSOseedValue column size must be D or less');
     end
     if ps < tmpsz(1)
         error('PSOseedValue row length must be # of particles or less');
     end
 end
 
% take care of setting max velocity and position params here
if length(mv)==1
 velmaskmin = -mv*ones(ps,D);     % min vel, psXD matrix
 velmaskmax = mv*ones(ps,D);      % max vel
elseif length(mv)==D     
 velmaskmin = repmat(forcerow(-mv),ps,1); % min vel
 velmaskmax = repmat(forcerow( mv),ps,1); % max vel
else
 error('Max vel must be either a scalar or same length as prob dimension D');
end
posmaskmin  = repmat(VR(1:D,1)',ps,1);  % min pos, psXD matrix
posmaskmax  = repmat(VR(1:D,2)',ps,1);  % max pos
posmaskmeth = 3; % 3=bounce method (see comments below inside epoch loop)

% PLOTTING
message = sprintf('PSO: %%g/%g iterations, GBest = %%20.20g.\n',me);
 
% INITIALIZE INITIALIZE INITIALIZE INITIALIZE INITIALIZE INITIALIZE
 
% initialize population of particles and their velocities at time zero,
% format of pos= (particle#, dimension)
 % construct random population positions bounded by VR
pos(1:ps,1:D) = normmat(rand([ps,D]),VR',1);
  
if PSOseed == 1         % initial positions user input, see comments above
    tmpsz                      = size(PSOseedValue);
    pos(1:tmpsz(1),1:tmpsz(2)) = PSOseedValue;  
end

 % construct initial random velocities between -mv,mv
vel(1:ps,1:D) = normmat(rand([ps,D]),[forcecol(-mv),forcecol(mv)]',1);

% initial pbest positions vals
pbest = pos;

% VECTORIZE THIS, or at least vectorize cost funct call 
for k=1:ps
    out(k,1) = feval(functname,pos(k,:));  % returns column of cost values (1 for each particle)
end
%---------------------------
 
pbestval=out;   % initially, pbest is same as pos

% assign initial gbest here also (gbest and gbestval)
if minmax==1
   % this picks gbestval when we want to maximize the function
    [gbestval,idx1] = max(pbestval);
elseif minmax==0
   % this works for straight minimization
    [gbestval,idx1] = min(pbestval);
elseif minmax==2
   % this works when you know target but not direction you need to go
   % good for a cost function that returns distance to target that can be either
   % negative or positive (direction info)
    [temp,idx1] = min((pbestval-ones(size(pbestval))*errgoal).^2);
    gbestval    = pbestval(idx1);
end

 % preallocate a variable to keep track of gbest for all iters
bestpos        = zeros(me,D+1)*NaN;
gbest          = pbest(idx1,:);  % this is gbest position
   % used with trainpso, for neural net training
   % assign gbest to net at each iteration, these interim assignments
   % are for plotting mostly
bestpos(1,1:D) = gbest;
 
% INITIALIZE END INITIALIZE END INITIALIZE END INITIALIZE END
% start PSO iterative procedures
cnt2   = 0; % counter used for the stopping subroutine based on error convergence
iwt(1) = iw1;
fid=fopen('tese.res','w');
fprintf(fid,'AQUI COMEÇA O PROGRAMA\n\n');
fprintf(fid,'MAXIMUM VELOCITY:                ');
fprintf(fid,'%4.2f  ', mv);
fprintf(fid,'\nMAX EPOCHS:                      %5.2f\n', P(1));
fprintf(fid,'POPULATION SIZE:                 %5.2f\n', P(2));
fprintf(fid,'ACC. CONSTANT 1 (LOCAL):         %5.2f\n', P(3));
fprintf(fid,'ACC. CONSTANT 2 (LOCAL):         %5.2f\n', P(4));
fprintf(fid,'INITIAL INERTIA CONSTANT:        %5.2f\n', P(5));
fprintf(fid,'FINAL INERTIA CONSTANT:          %5.2f\n', P(6));
fprintf(fid,'ITERATION FINAL IN. CONST.:      %5.2f\n', P(7));
fprintf(fid,'MINIMAL CONVERGENCE SPEED:       %5.2f\n', P(8));
fprintf(fid,'MINIMAL EPOCHS WITHOUT CHANGE:   %5.2f\n', P(9));
fprintf(fid,'ERROR GOAL:                      %5.2f\n', P(10));
fprintf(fid,'MANUAL INITIAL:                  %5.2f\n\n\n', P(11));
fclose(fid);
for i=1:me  % start epoch loop (iterations)
    fid=fopen('tese.res','a');
    fprintf(fid, 'SCORES\n');
    fprintf(fid, '%10.5f  \n', out);
    fprintf(fid, '\nGBESTVAL\n%10.6f\n\nGBEST\n', gbestval);
    fprintf(fid, '%10.5f\n',gbest);
    fclose(fid);
    disp('---------------ITEROU!!!---------------');
    disp('---------------ITEROU!!!---------------');
    disp('---------------ITEROU!!!---------------');
    for k=1:ps
        out(k,1) = feval(functname,pos(k,:));
    end
    outbestval = feval(functname,gbest);

    tr(i+1)          = gbestval; % keep track of global best val
    te               = i; % returns epoch number to calling program when done
    bestpos(i,1:D+1) = [gbest,gbestval];
     
     %assignin('base','bestpos',bestpos(i,1:D+1));

   
    % find particles where we have new pbest, depending on minmax choice 
    % then find gbest and gbestval
     %[size(out),size(pbestval)]

     if minmax == 0
        [tempi]            = find(pbestval>=out); % new min pbestvals
        pbestval(tempi,1)  = out(tempi);   % update pbestvals
        pbest(tempi,:)     = pos(tempi,:); % update pbest positions
       
        [iterbestval,idx1] = min(pbestval);
        
        if gbestval >= iterbestval
            gbestval = iterbestval;
            gbest    = pbest(idx1,:);
        end
     elseif minmax == 1
        [tempi,dum]        = find(pbestval<=out); % new max pbestvals
        pbestval(tempi,1)  = out(tempi,1); % update pbestvals
        pbest(tempi,:)     = pos(tempi,:); % update pbest positions
 
        [iterbestval,idx1] = max(pbestval);
        if gbestval <= iterbestval
            gbestval = iterbestval;
            gbest    = pbest(idx1,:);
        end
     elseif minmax == 2  % this won't work as it is, fix it later
        egones            = errgoal*ones(ps,1); % vector of errgoals
        sqrerr2           = ((pbestval-egones).^2);
        sqrerr1           = ((out-egones).^2);
        [tempi,dum]       = find(sqerr1 <= sqrerr2); % find particles closest to targ
        pbestval(tempi,1) = out(tempi,1); % update pbestvals
        pbest(tempi,:)    = pos(tempi,:); % update pbest positions

        sqrerr            = ((pbestval-egones).^2); % need to do this to reflect new pbests
        [temp,idx1]       = min(sqrerr);
        iterbestval       = pbestval(idx1);
        
        if (iterbestval-errgoal)^2 <= (gbestval-errgoal)^2
           gbestval = iterbestval;
           gbest    = pbest(idx1,:);
        end
     end
    
     %PSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSO

      % get new velocities, positions (this is the heart of the PSO algorithm)     
      % each epoch get new set of random numbers
       rannum1 = rand([ps,D]); % for Trelea and Clerc types
       rannum2 = rand([ps,D]);       

        % common PSO algo with inertia wt 
        % get inertia weight, just a linear funct w.r.t. epoch parameter iwe
         if i<=iwe
            iwt(i) = ((iw2-iw1)/(iwe-1))*(i-1)+iw1;
         else
            iwt(i) = iw2;
         end
        % random number including acceleration constants
         ac11 = rannum1.*ac1;    % for common PSO w/inertia
         ac22 = rannum2.*ac2;
         
         vel = iwt(i).*vel...                             % prev vel
               +ac11.*(pbest-pos)...                      % independent
               +ac22.*(repmat(gbest,ps,1)-pos);           % social                  
       
       % limit velocities here using masking
        vel = ( (vel <= velmaskmin).*velmaskmin ) + ( (vel > velmaskmin).*vel );
        vel = ( (vel >= velmaskmax).*velmaskmax ) + ( (vel < velmaskmax).*vel );     
        
       % update new position (PSO algo)    
        pos = pos + vel;
    
       % position masking, limits positions to desired search space
       % method: 0) no position limiting, 1) saturation at limit,
       %         2) wraparound at limit , 3) bounce off limit
        minposmask_throwaway = pos <= posmaskmin;  % these are psXD matrices
        minposmask_keep      = pos >  posmaskmin;     
        maxposmask_throwaway = pos >= posmaskmax;
        maxposmask_keep      = pos <  posmaskmax;
     
        if     posmaskmeth == 1
         % this is the saturation method
          pos = ( minposmask_throwaway.*posmaskmin ) + ( minposmask_keep.*pos );
          pos = ( maxposmask_throwaway.*posmaskmax ) + ( maxposmask_keep.*pos );      
        elseif posmaskmeth == 2
         % this is the wraparound method
          pos = ( minposmask_throwaway.*posmaskmax ) + ( minposmask_keep.*pos );
          pos = ( maxposmask_throwaway.*posmaskmin ) + ( maxposmask_keep.*pos );                
        elseif posmaskmeth == 3
         % this is the bounce method, particles bounce off the boundaries with -vel      
          pos = ( minposmask_throwaway.*posmaskmin ) + ( minposmask_keep.*pos );
          pos = ( maxposmask_throwaway.*posmaskmax ) + ( maxposmask_keep.*pos );

          vel = (vel.*minposmask_keep) + (-vel.*minposmask_throwaway);
          vel = (vel.*maxposmask_keep) + (-vel.*maxposmask_throwaway);
        else
         % no change, this is the original Eberhart, Kennedy method, 
         % it lets the particles grow beyond bounds if psoparams (P)
         % especially Vmax, aren't set correctly, see the literature
        end

     %PSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSO
% check for stopping criterion based on speed of convergence to desired 
   % error   
    tmp1 = abs(tr(i) - gbestval);
    if tmp1 > ergrd
       cnt2 = 0;
    elseif tmp1 <= ergrd
       cnt2 = cnt2+1;
       if cnt2 >= ergrdep
         
          fprintf(message,i,gbestval);           
          disp(' ');
          disp(['--> Solution likely, GBest hasn''t changed by at least ',...
              num2str(ergrd),' for ',...
                  num2str(cnt2),' epochs.']);  
                
         break
       end
    end
    
   % this stops if using constrained optimization and goal is reached
    if ~isnan(errgoal)
     if ((gbestval<=errgoal) & (minmax==0)) | ((gbestval>=errgoal) & (minmax==1))  

         
             fprintf(message,i,gbestval);
             disp(' ');            
             disp(['--> Error Goal reached, successful termination!']);
             
         
         break
     end
     
    % this is stopping criterion for constrained from both sides    
     if minmax == 2
       if ((tr(i)<errgoal) & (gbestval>=errgoal)) | ((tr(i)>errgoal) ...
               & (gbestval <= errgoal))        
         
             fprintf(message,i,gbestval);
             disp(' ');            
             disp(['--> Error Goal reached, successful termination!']);            
             
         
         break              
       end
     end % end if minmax==2
    end  % end ~isnan if

 %    % convert back to inertial frame
 %     pos = pos - repmat(gbestoffset,ps,1);
 %     pbest = pbest - repmat(gbestoffset,ps,1);
 %     gbest = gbest + gbestoffset;
  

end  % end epoch loop

%% clear temp outputs
% evalin('base','clear temp_pso_out temp_te temp_tr;');

% output & return
 OUT=[gbest';gbestval];
 varargout{1}=te;
 varargout{2}=[tr(find(~isnan(tr)))];
 
 return