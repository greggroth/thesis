function [ total ] = ImportSegmentedT1(varargin)
%  ImportSegmentedT1 Import NII files from a directory
%   Must be run within the directory containing the files
%
%   Output: head data as single with variables stored in the 4th dimension.
%
%   Author:  Greggory Rothmeier (greggroth@gmail.com)
%   Georgia State University
%   Created:  5/31/11

statusbar = waitbar(0,'Initializing');

if size(varargin) == 1
    oldFolder = cd(varargin{1});
end


% =====================
% = Tissue Parameters =
% =====================
% Each tissue type is assigned an integer index (i.e. gray matter -> 11) such that
% tissue-specific parameters can be found by looking at that element within the 
% corresponding storage matrix (i.e. QmSTORE(11) -> gray matter Qm)

% Parameters taken from Colins, 2004

tisorder = [11 15 5 13 3];  %  Using:  [GM WM CSF Muscle Bone]

QmSTORE = [0 0 26.1 11600 0 26.1 697 0 0 302 15575 0 697 1100 5192];
cSTORE = [1006 4600 2110 3640 3800 1300 3720 3000 4200 2300 3680 3500 3720 3150 3600];
rhoSTORE = [1.3 1057 1080 1035.5 1007 1850 1126 1076 1009 916 1035.5 1151 1041 1100 1027.4];
kSTORE = [0.026 0.51 0.65 0.534 0.5 0.65 0.527 0.4 0.594 0.25 0.565 0.4975 0.4975 .342 .503];
wSTORE = [0 1000 3 45.2 0 1.35 40 0 0 2.8 67.1 3.8 3.8 12 23.7];

% =====================================
% = Import the pre-segmented T1 files =
% =====================================
% The T1 contrast image sould be segmented using SPM8.
%   This loop needs to complete before the next one can begin
%  Import all of the datat and store as 'cdat1','cdat2', etc.
for i = 1:5
    eval(strcat('dat',num2str(i),' = loadNII(''rc', num2str(i), 'single_subj_T1.nii'');'))  
    %  Preallocate 
    eval(strcat('out', num2str(i),' = zeros(cat(2,size(dat', num2str(i),'),7));'))
end

% ============================
% = Populate the head matrix =
% ============================
%   For each data file, it fills in the data from the data storage arrays
%   for that particular type of tissue.  It picks which ever tissue is the 
%   most likely candidate for that voxel based on the segmented data

%   PROBLEM:  It returns 0 (later filled with air) if there is equal
%   probability of a voxel being two or more different types of tissue.  
%   SOLVED BY fillHoles()


for i = 1:5
    %  Preallocate
    holder = zeros(cat(2,size(dat1),7),'single');         
    mask = zeros(size(dat1));
    final = zeros(size(holder),'single');
    
    %  Create a mask that indicates whether it is the mostly likely tissue type
    guide = [1 2 3 4 5 1 2 3 4 5];  %  This guides it through the data correctly
    eval(strcat('mask = (dat',num2str(i),'>dat',num2str(guide(i+1)),') & (dat',num2str(i),'>dat',num2str(guide(i+2)),') & (dat',num2str(i),'>dat',num2str(guide(i+3)),') & (dat',num2str(i),'>dat',num2str(guide(i+4)),') & (dat',num2str(i),'~=0);'))
    holder(:,:,:,1) = mask;                      %  move structure data to new matrix
    a = find(holder(:,:,:,1) == 1);              %  get indicies of tissues
    [x y z t] = ind2sub(size(holder),a);         %  gets coordinates from index
   
    %  go to each tissue point and store the info
    for j = 1:length(a)                        
        final(x(j),y(j),z(j),:) = [tisorder(i) 0 QmSTORE(tisorder(i)) cSTORE(tisorder(i)) rhoSTORE(tisorder(i)) kSTORE(tisorder(i)) wSTORE(tisorder(i))];
    end
    
    %  Saves the result to a unique output variable (out1, out2, etc)
    eval(strcat('out',num2str(i),'= final;'))  
    
    clearvars a x y z t holder final;
    waitbar(i/6,statusbar,sprintf(['File ',num2str(i),' Import Compete']));
end    

% The filleAir() function checks for any voxels which were not assigned a 
% tissue type and fills them in with air
almostthere = fillAir(out1+out2+out3+out4+out5);
% The fillHoles() function corrects for a voxel having two equally-probable 
% tissue types
total = single(buildskin(fillHoles(dat1,dat2,dat3,dat4,dat5,almostthere)));
waitbar(1,statusbar,'Saving Data')

cd(oldFolder);
close(statusbar);

end

