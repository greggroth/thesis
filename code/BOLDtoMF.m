function [ ] = BOLDtoMF( varargin)
%BOLDtoMF Calculate metabolism and blood from from BOLD reponse
%
%   Input: Directory containing a series of *.nii files of the BOLD
%   response.  
%
%   Output: Two new files will be created in a new subdirectory with a
%   variable for each time step.  
%
%   Usage:                              
%       BOLDtoMF 
%       BOLDtoMF(directory)
%
%   If a directory is not provided, one will be requested.  
%
%   From Sotero, et. al. 2010

%%  Setup

%  Check input
switch length(varargin)   % if a folder isn't an argument, it'll prompt for one
    case 0
        fold_name = uigetdir;
        if ~fold_name  % Cancel Button
            return
        end
    case 1
        fold_name = varargin{1};
    otherwise
        error('Input is not understood')
end

%  Go to the folder containing the files
oldfold = cd(fold_name);
file_list = dir('*.nii');
file_count = length(file_list);

%  Set up a directory for the outputs
newFolder = ['Output_',datestr(clock,1)];
mkdir(newFolder)

%  Make *.mat files to append the data to
m0001 = 0; f0001 = 0;
save(['./' newFolder '/m.mat'],'m0001');
save(['./' newFolder '/f.mat'],'f0001');



%%  Norm
s = loadNII(file_list(1).name);
norm = ones(size(s));

%%  Calculate
%
% This will calculate the metabolism and blood flow.  The output is
% appended to 'm.mat' and 'f.mat' within a new folder created within the
% directory containing the data.

statusbar = waitbar(0,'Initializing');

maxBOLD = 0.22;
%{
%%  Find the max BOLD response
for j=1:file_count
    try 
        waitbar(j/file_count, statusbar, sprintf('Finding max change in BOLD: %d%%',round((j/file_count)*100)));
    catch err
        return
    end
    s = loadNII(file_list(j).name);  %  Load up the file
    if max(s(:)) > maxBOLD           %  if the max value beats the current max, take it
        maxBOLD = max(s(:));
        disp([j maxBOLD])
    end
end
%}
%  Required Parameters
p = [0.4 1.5 0.1870 0.1572 -0.6041 maxBOLD]; % [alpha beta a b c A]

%%  Calc flow and metabolism (when BOLD = 1)
%   I thought that the equations should work out so that an input of s = 1
%   returns f and m = 1, but until I sort that out here is a cheating work
%   around.  Make sure this is valid before publishing.

s = 0;
y = -((p(4)*p(2))/(p(1)+p(2)*p(5)))*((p(6)-s)/(p(6)*p(3)^p(2)))^(1/(p(1)+p(2)*p(5)));
fNOACT = -((p(1)+p(2)*p(5))/(p(4)*p(2)))*lambertw(y);
mNOACT = p(3)*fNOACT^(p(5)+1)*exp(-p(4)*fNOACT);


%%  Calc flow and metabolism
disp(fold_name)
for j=1:file_count
    tic
    %avgtime = mean(timelist);
    %disp(avgtime)
    %timeremaining = (file_count-j)*avgtime;
    try
        waitbar(j/file_count,statusbar,sprintf('%d%%',round((j/file_count)*100)));
    catch err
        return
    end
    s = loadNII(file_list(j).name);  %  Load up the file
    s(isnan(s)) = 1; %what to do with NaNs and INFS?  Not sure. maybe set to zero for now.
    s(isinf(s)) = 1;
    y = -((p(4)*p(2))/(p(1)+p(2)*p(5))).*((p(6)-s)./(p(6)*p(3)^p(2))).^(1/(p(1)+p(2)*p(5)));
    if (size(y,1)==91)&&(size(y,2)==109)&&(size(y,3)==91)
        f = -((p(1)+p(2)*p(5))/(p(4)*p(2))).*lambw_mex(real(y));  % <-- compiled version: runs faster
    else
        f = -((p(1)+p(2)*p(5))/(p(4)*p(2))).*lambw(y);  % <-- not compiled, but still pretty fast
    end
    m = p(3)*f.^(p(5)+1).*exp(-p(4)*f);
    %  Clean up NaNs
    m(isnan(m))=1;
    f(isnan(f))=1;
    %  make sure that if the BOLD was 1 then the meabolism/flow is 1
    %  DOUBLE CHECK THAT THIS IS OK!!!!!!!
    m = m./mNOACT;
    f = f./fNOACT;
    
    eval(['m' sprintf('%04d',j) ' = m;']);
    eval(['f' sprintf('%04d',j) ' = f;']);
    eval(['save(''./' newFolder '/m.mat'', ''m' sprintf('%04d',j) ''',''-append'');']);
    eval(['save(''./' newFolder '/f.mat'', ''f' sprintf('%04d',j) ''',''-append'');']);
    clear m0* f0*   % prevent holding onto variables after they're done being used.
    
    t = toc;
    rem = ((file_count-j)*t)/60;
    disp([file_list(j).name '    ' num2str(rem,4) ' minutes remaining'])
end

close(statusbar)
cd(oldfold)


%%  OLD METHOD

%{
tic
y = -((p(4)*p(2))/(p(1)+p(2)*p(5))).*((p(6)-s)./(p(6)*p(3)^p(2))).^(1/(p(1)+p(2)*p(5)));
f = -((p(1)+p(2)*p(5))/(p(4)*p(2))).*lambertw(y);
m = p(3)*f.^(p(5)+1).*exp(-p(4)*f);
toc
%}

%f_out = f;
%m_out = m;
% m = s;
% f = s;
%%  Output
% In order to make it easier when calculating the change in temperature,
% this function will create one *.mat file with a seperate variable for
% each time step.  It's a little anoying but since it's such a large file
% when combined, it's the only way to do it.  

%{
newFolder = ['Output ',datestr(clock)];
mkdir(newFolder)
oldFolder = cd(newFolder);
varsM = cell(size(s,1), 1);
varsF = cell(size(s,1), 1);
for k = 1:size(s,1);
    eval(strcat('m',num2str(k),' = squeeze(m(k,:,:,:));'));
    eval(strcat('f',num2str(k),' = squeeze(f(k,:,:,:));'));
    varsM{k} = strcat(',''m',num2str(k),''' ');
    varsF{k} = strcat(',''f',num2str(k),''' ');
end
mfin = strcat(cell2mat(varsM'));
ffin = strcat(cell2mat(varsF'));
eval(strcat('save(''m_BOLD.mat''',mfin,');'));
eval(strcat('save(''f_BOLD.mat''',ffin,');'));
cd(oldFolder);
%}

end

