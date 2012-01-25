function temperatureOut = tempCalcDynMF(tissue,bloodT,airT,nt,tmax,pastCalc,metab,flow,savesteps,region)
% tempCalcChaning Metabolism  How does changin metabolism 
% affect things?
% 
%   tissue: holds all of the strucual information
%   bloodT: Temperature of the blood
%   airT:   Temperature of the surrounding ait
%   nt:     Number of time steps
%   tmax:   Total amount of time the simulation should run over
%
%   region: logical matrix same size as head
%
%   Writen by Greggory Rothmeier (greggroth@gmail.com)
%   Georgia State University Dept. Physics and Astronomy
%   May, 2011

statusbar = waitbar(0,'Initializing');

%%   Default Values
if nargin<2,  bloodT = 37;          end
if nargin<3,  airT = 24;            end
if nargin<4,  nt = 3;               end
if nargin<5,  tmax = 1;             end   
if nargin<6,  pastCalc = 0;         end 
  

% Length of one side of a voxel (m)
dx = 2*10^-3;

if nt<(2*tmax),
   warning('Time step size is not large enough.  Results will be unreliable.  Consider increasing the number of steps or reducing tmax.')
end


% Constants used that aren't already stored in tissue
[xmax ymax zmax t] = size(tissue); 
clear t;
dt = ones([xmax ymax zmax])*(tmax/(nt-1));
% rhoBlood = 1057;
% wBlood = 1000;
% cBlood = 3600;

%%  Determine Metab/Flow Data Storage System
if ischar(metab)&&ischar(flow)  
  % if file locations are given rather than data
    option = 1;
else
  % Preallocate matrices for holding metabolism and blood flow data
    metabMulti = ones([xmax ymax zmax],'single');
    flowMulti = ones([xmax ymax zmax],'single');
    option = 0;
end

%%  Maps
% Creates a map that identifies where there is tissue
% the condition squeeze(tissue(:,:,:,)~=airIndex picks out the 
% elements that are tissue

tmax = ceil((nt-1)/savesteps);
temperatureOut = ones(tmax,xmax,ymax,zmax,'single');
temperature = ones(2,xmax,ymax,zmax,'single')*airT;
if pastCalc == 0
    temperature(1,squeeze(tissue(:,:,:,1))~=1) = bloodT;
else
  % Starts everything off at the pre-determined equilibium temperatures
    temperature(1,:,:,:) = pastCalc(end,:,:,:);
end
temperatureOut(1,:,:,:) = temperature(1,:,:,:);


% ===========
% = Do Work =
% ===========
%   This is a vectorized version of the next section.  For the love of
%   god don't make any changes to this without first looking below to
%   make sure you know what you're changing.  This is [nearly]
%   impossible to understand, so take your time and don't break it.
%   data is stored in 'tissue' as such : 
%  [tissuetype 0 Qm c rho k w] <--  second element is blank for all. 
%  [    1      2  3 4  5  6 7]

% This makes an array that has smoothed out variations in k by 
% averaging all of the k's around each voxel (including itself). This 
% is a hap-hazard solution to the problem that if you only take the 
% value of k for the voxel without considering what surrounds it, it 
% doesn't matter whether the head is surrounded by air, water or 
% anything else. Since water is a better thermal conductor than air, we 
% need a way of accounting for this. This is one way:

averagedk = (circshift(tissue(:,:,:,6),[1 0 0])+circshift(tissue(:,:,:,6),[-1 0 0])+circshift(tissue(:,:,:,6),[0 1 0])+circshift(tissue(:,:,:,6),[0 -1 0])+circshift(tissue(:,:,:,6),[0 0 1])+circshift(tissue(:,:,:,6),[0 0 -1])+tissue(:,:,:,6))/7;
rhoblood = 1057;
cblood = 3600;

%%  Only saves every 4 steps to reduce the final matrix size
for t2 = 1:nt-1
   waitbar(t2/(nt-1),statusbar,sprintf('%d%%',round(t2/(nt-1)*100)));
   
% if a variable needs to be used multiple times for the same time step.
   t3 = floor((t2-1)/4)+1;  % 1 1 1 1 2 2 2 2 3 3 . . .
   
   % if a file is specified, pulls the data from the file for each step
   if option  
       eval(strcat('load(fullfile(metab),''-mat'',''m',sprintf('%04d',t3),''');'));
       eval(strcat('load(fullfile(flow),''-mat'',''f',sprintf('%04d',t3),''');'));
       eval(strcat('metabMulti = m',sprintf('%04d',t3),';'));
       eval(strcat('flowMulti = f',sprintf('%04d',t3),';'));
       eval(strcat('clear f', sprintf('%04d',t3),' m',sprintf('%04d',t3)))
   else
       metabMulti(region) = metab(t2);   %  region is hardcoded here
       flowMulti(region) = flow(t2); 
   end

   temperature(2,:,:,:) = squeeze(temperature(1,:,:,:)) + ...
        dt./(tissue(:,:,:,5).*tissue(:,:,:,4)).* ...
        ((averagedk/dx^2).*...
        (circshift(squeeze(temperature(1,:,:,:)),[1 0 0])-2*squeeze(temperature(1,:,:,:))+circshift(squeeze(temperature(1,:,:,:)),[-1 0 0])+...  % shift along x
         circshift(squeeze(temperature(1,:,:,:)),[0 1 0])-2*squeeze(temperature(1,:,:,:))+circshift(squeeze(temperature(1,:,:,:)),[0 -1 0])+...  % shift along y
         circshift(squeeze(temperature(1,:,:,:)),[0 0 1])-2*squeeze(temperature(1,:,:,:))+circshift(squeeze(temperature(1,:,:,:)),[0 0 -1]))...  % shift along z
            -(1/6000)*rhoblood*flowMulti.*tissue(:,:,:,7)*cblood.*(squeeze(temperature(1,:,:,:))-bloodT)+metabMulti.*tissue(:,:,:,3));
    % resets the air temperature back since it's also modified above, 
    % but it needs to be kept constant throughout the calculations
    temperature(2,squeeze(tissue(:,:,:,1))==1) = airT; 
    temperatureOut(ceil(t2/savesteps),:,:,:) = temperature(2,:,:,:);
    temperature(1,:,:,:) = temperature(2,:,:,:);  % moves 2 back to 1
    clear metabMulti flowMulti
end
close(statusbar);

% ============
% = Old Code =
% ============
% This is what used to be used. It's much slower (~60 times slower), 
% but it's much easier to understand compared to the above code. If any
% changes need to be made above, first look through this code to ensure
% you understand what's happening before making changes. It's really 
% easy to mess up the code above and nearly impossible to figure out 
% where.
% 
%  good luck.

% for t2 = 1:nt-1
%     for x2 = 2:xmax-1
%         for y2 = 2:ymax-1
%             for z2 = 2:zmax-1
%                 if tissue(x2,y2,z2,1) ~= 1,
%                     temperature(t2+1,x2,y2,z2) = temperature(t2,x2,y2,z2) + (dt/(tissue(x2,y2,z2,5)*tissue(x2,y2,z2,4)))*((tissue(x2,y2,z2,6)/dx^2)*...
%                       (temperature(t2,x2+1,y2,z2)-2*temperature(t2,x2,y2,z2)+temperature(t2,x2-1,y2,z2)+...
%                       temperature(t2,x2,y2+1,z2)-2*temperature(t2,x2,y2,z2)+temperature(t2,x2,y2-1,z2)+...
%                       temperature(t2,x2,y2,z2+1)-2*temperature(t2,x2,y2,z2)+temperature(t2,x2,y2,z2-1))...
%                       -(1/6000)*rhoBlood*wBlood*cBlood*(temperature(t2,x2,y2,z2)-bloodT)+tissue(x2,y2,z2,3));
%                 end
%             end
%         end
%     end
% end

end