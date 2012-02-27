function temperature_Out = tempCalcEquillibrium(tissue,bloodT,airT,nt,tmax,pastCalc,printprogress)
% tempCalcEquillibrium  Find the equillibrium values
%   tissue: holds all of the strucual information
%   bloodT: Temperature of the blood
%   airT:   Temperature of the surrounding ait
%   nt:     Max number of time steps
%   tmax:   Total amount of time the simulation should run over
%
%   This is based off of tempCalc() but loops until the rate of 
% change of a each voxel is sufficiently small then outputs 
% what's calculated.  If if takes too long to do all at once, 
% split it up into smaller time chunks and use the last step 
% from the previous dataset as pastCalc in order to resume.
%
%   Note: This does not save the time corse because it can take 
% a lot of step to find the equillibrium.  It outputs the last 
% time step.
%
%   Writen by Greggory Rothmeier (greggroth@gmail.com)
%   Georgia State University Dept. Physics and Astronomy
%   May, 2011
tic
%%   Default Values
if nargin<2, bloodT = 37;       end
if nargin<3, airT = 24;         end
if nargin<4, nt = 100;          end
if nargin<5, tmax = 50;         end   
if nargin<6, pastCalc = 0;      end
if nargin<7, printprogress = 1; end

% These rescue the data if the calculation is interrupted.
global temperature
global dirty

c = onCleanup(@InterCatch);
dirty = 1;

dx = 2*10^-3;        %  Voxel size (m)

if nt<(2*tmax),
   warning('Time step size is not large enough.  Results will be unreliable.  Consider increasing the number of steps or reducing tmax.')
end


% Constants used that aren't already stored in tissue
[xmax ymax zmax t] = size(tissue); 
clear t;
dt = tmax/(nt-1);
% rhoBlood = 1057;
% wBlood = 1000;
% cBlood = 3600;

% =========
% = Setup =
% =========
%   Starts all tissue voxels at bloodT (default 37) and maintains 
% air at airT (default 24)
%   The condition squeeze(tissue(:,:,:,)~=airIndex picks out the 
% elements that are tissue

temperature = ones(3,xmax,ymax,zmax,'single')*airT;
if pastCalc == 0
    temperature(1,squeeze(tissue(:,:,:,1))~=1) = bloodT;
else
    temperature(1,:,:,:) = pastCalc;
end
numElements = numel(temperature(1,:,:,:));

% ===========
% = Do Work =
% ===========
%   This is a vectorized version of the next section.  For the love 
% of god don't make any changes to this without first looking below
% to make sure you know what you're changing.  This is [nearly] 
% impossible to understand, so take your time and don't break it.
%   data is stored in 'tissue' as such :
%   [tissuetype 0 Qm c rho k w]; <-- second element is blank for all.
%   [    1      2  3 4  5  6 7

averagedk = (circshift(tissue(:,:,:,6),[1 0 0])+circshift(tissue(:,:,:,6),[-1 0 0])+circshift(tissue(:,:,:,6),[0 1 0])+circshift(tissue(:,:,:,6),[0 -1 0])+circshift(tissue(:,:,:,6),[0 0 1])+circshift(tissue(:,:,:,6),[0 0 -1])+tissue(:,:,:,6))/7;
rhoblood = 1057;
cblood = 3600;

%%  Specify Percision Goal
tolerence = 1;     % fraction of voxels have a slope less than 'zeropoint'
zeropoint = 2.5e-7;  % point at which the slope between two *steps* is considered essentially zero


goal = numElements - tolerence*numElements;
goon = numElements; %  Forces the while loop to run the first time
format shortG;
% temperature(1,:,:,:) = Current Temperature
% temperature(2,:,:,:) = Next Temperature
% Resets after each update
if printprogress
    disp(['Goal:  ', num2str(goal),' remaining voxels'])
end
t2 = 1;
while goon(1)>goal && t2<=nt  % runs until either 'goal' elements have a slope greater than 'zeropoint' or it exceeds nt
   if printprogress
    disp([t2 goon(1) ((numElements-goon(1))/numElements)*100]) % progress
   end
   temperature(2,:,:,:) = squeeze(temperature(1,:,:,:)) + ...
        dt/(tissue(:,:,:,5).*tissue(:,:,:,4)).* ...
        ((averagedk/dx^2).*...
        (circshift(squeeze(temperature(1,:,:,:)),[1 0 0])-2*squeeze(temperature(1,:,:,:))+circshift(squeeze(temperature(1,:,:,:)),[-1 0 0])+...  % shift along x
         circshift(squeeze(temperature(1,:,:,:)),[0 1 0])-2*squeeze(temperature(1,:,:,:))+circshift(squeeze(temperature(1,:,:,:)),[0 -1 0])+...  % shift along y
         circshift(squeeze(temperature(1,:,:,:)),[0 0 1])-2*squeeze(temperature(1,:,:,:))+circshift(squeeze(temperature(1,:,:,:)),[0 0 -1]))...  % shift along z
            -(1/6000)*rhoblood*tissue(:,:,:,7)*cblood.*(squeeze(temperature(1,:,:,:))-bloodT)+tissue(:,:,:,3));
    %   resets the air temperature back since it's also modified 
    % above, but it needs to be kept constant throughout the 
    % calculations
    temperature(2,squeeze(tissue(:,:,:,1))==1) = airT;
    %   checks how quickly the temperature is changing and if it is 
    % close enough to zero to be considered stopped ('zeropoint')
    goon = size(temperature(abs(squeeze(temperature(2,:,:,:)-temperature(1,:,:,:)))>zeropoint));
    temperature(1,:,:,:) = temperature(2,:,:,:);    
    t2 = t2 + 1;
end

temperature_Out = temperature(2,:,:,:);
dirty = 0;

% equilTemperature = temperature_Out;
% save('equil.mat','equilTemperature');

%% To Combine Datasets
%  use this technique if there are seperate datasets that need 
% combining
%    vertcat(squeeze(res1(:,:,:,:)),squeeze(res2(2:end,:,:,:)))
%  Where for all by the first dataset, you need to do the time from 
% 2:end so that there are no repeats (remember that the last 
% timestep from the previous dataset serves as the first for the 
% new one)


time = toc;
end

% Recovers the data if calculation was interrupted
function InterCatch
global dirty
if dirty
    disp('Interupt Intercepted.  Inprepretating Interworkspace Data Interfaces.')
    global temperature
    equillibriumT = temperature;
    save('equiltempAbortDump.mat','equillibriumT');
    % setappdata(0,'InterpOut',temperature);
end
end

