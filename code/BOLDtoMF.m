function [ ] = BOLDtoMF( varargin)
%BOLDtoMF Calculate metabolism and blood from from BOLD reponse
%
%   Input: Directory containing a series of *.nii files of the BOLD
%   response.  
%
%   Output: Two new files will be created in a new subdirectory 
%   with a variable for each time step.  
%
%   Usage:                              
%       BOLDtoMF 
%       BOLDtoMF(directory)
%
%   If a directory is not provided, one will be requested.  
%
%   Method from Sotero, et. al. 2011

% =========
% = Setup =
% =========
% if a folder isn't an argument, it'll prompt for one
switch length(varargin)   
    case 0
        fold_name = uigetdir;
        if ~fold_name  % Cancel Button pressed
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

s = loadNII(file_list(1).name);
norm = ones(size(s));

% ===========
% = Do Work =
% ===========
% This will calculate the metabolism and blood flow.  The output is
% appended to 'm.mat' and 'f.mat' within a new folder created 
% within the directory containing the data.

statusbar = waitbar(0,'Initializing');

maxBOLD = 0.22;

% Required Parameters
%   [alpha beta a      b      c      A    ]
p = [0.4 1.5 0.1870 0.1572 -0.6041 maxBOLD]; 

% Calc flow and metabolism for when BOLD = 1
s = 0;
y = -((p(4)*p(2))/(p(1)+p(2)*p(5)))*((p(6)-s)/(p(6)*p(3)^p(2)))^(1/(p(1)+p(2)*p(5)));
fNOACT = -((p(1)+p(2)*p(5))/(p(4)*p(2)))*lambertw(y);
mNOACT = p(3)*fNOACT^(p(5)+1)*exp(-p(4)*fNOACT);


%%  Calc flow and metabolism
disp(fold_name)
for j=1:file_count
    try 
      waitbar(j/file_count, statusbar, sprintf('%d%%', round((j/file_count)*100)));
    catch err
        return
    end
    s = loadNII(file_list(j).name);  %  Load up the file
    s(isnan(s)) = 1;
    s(isinf(s)) = 1;
    y = -((p(4)*p(2))/(p(1)+p(2)*p(5))).*((p(6)-s)./(p(6)*p(3)^p(2))).^(1/(p(1)+p(2)*p(5)));
    if (size(y,1)==91)&&(size(y,2)==109)&&(size(y,3)==91)
        f = -((p(1)+p(2)*p(5))/(p(4)*p(2))).*lambw_mex(real(y));
    else
        f = -((p(1)+p(2)*p(5))/(p(4)*p(2))).*lambw(y);
    end
    m = p(3)*f.^(p(5)+1).*exp(-p(4)*f);
    % Clean up NaNs that may have popped up
    m(isnan(m))=1;
    f(isnan(f))=1;
    % Normalize to resting m and f
    m = m./mNOACT;
    f = f./fNOACT;
    
    % Rename and save the data
    eval(['m' sprintf('%04d',j) ' = m;']);
    eval(['f' sprintf('%04d',j) ' = f;']);
    eval(['save(''./' newFolder '/m.mat'', ''m' sprintf('%04d',j) ''',''-append'');']);
    eval(['save(''./' newFolder '/f.mat'', ''f' sprintf('%04d',j) ''',''-append'');']);
    clear m0* f0*
end

close(statusbar)
cd(oldfold)
end