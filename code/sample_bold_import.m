%%================================================================
%% How to process preprocessed BOLD data to calculate temperature
%%================================================================

% This Matlab script was used to automate the the process of using 
% BOLD data stored in NIFTI (*.nii) format to calculate temperature 
% changes.  The particulars of the code may be specific to this 
% case, but the procedure should be the same when doing these 
% calculations on other datasets.  All required functions are 
% included as an attachment to my thesis and are available on my 
% github (https://github.com/greggroth/tempcalc)

cd('/Users/Greggory/Documents/Data/fmri_rhythmic_tapping01/NIFTI')

directories = dir('*01');

%%  Move coregistered files to new Directory
for i = 1:length(directories)
    dir_name = directories(i).name;
    main_path = cd( [dir_name filesep dir_name '_NIFTI_1'] );
    mkdir 'Coregistered'
    movefile('r*.nii','Coregistered')
    main_path = cd( [dir_name filesep dir_name '_NIFTI_2'] );
    mkdir 'Coregistered'
    movefile('r*.nii','Coregistered')
    cd(main_path)
end

%%  Calculate Rest State
disp('Calculating Rest State')
for i = 1:length(directories)
    dir_name = directories(i).name;
    avg_NII_rest([dir_name filesep dir_name '_NIFTI_1' filesep 'Coregistered']);
    avg_NII_rest([dir_name filesep dir_name '_NIFTI_2' filesep 'Coregistered']);
end


%%  Normalize to Rest and Mask
disp('Normalize to Rest and Mask')
for i = 1:length(directories)
    dir_name = directories(i).name;
    avg_NII_normalize([dir_name filesep dir_name '_NIFTI_1' filesep 'Coregistered'], fullfile(dir_name, [dir_name '_NIFTI_1'], 'Coregistered', 'RestState', 'RestStateAvg.nii'), 'fullBrainMask.nii');
    avg_NII_normalize([dir_name filesep dir_name '_NIFTI_2' filesep 'Coregistered'], fullfile(dir_name, [dir_name '_NIFTI_2'], 'Coregistered', 'RestState', 'RestStateAvg.nii'), 'fullBrainMask.nii');
end


%%  Calculate metabolism and blood flow change
disp('Calculate metabolism and blood flow change')
for i = 1:length(directories)
    dir_1 = [ directories(i).name filesep  directories(i).name '_NIFTI_1' filesep 'Coregistered' filesep 'Normalized_to_rest'];
    dir_2 = [ directories(i).name filesep  directories(i).name '_NIFTI_2' filesep 'Coregistered' filesep 'Normalized_to_rest'];
    BOLDtoMF(dir_1);
    BOLDtoMF(dir_2);
end


%%  Calculate the change in temperature based on metabolism and 
%   blood flow

% load('equil.mat');  % equillibriumT
% load('tt_headdata.mat');  % headdata
mask = loadNII('fullBrainMask.nii');

for i = 1:length(directories)
    disp([int2str(i) '-1 started'])
    tic
    % Part I
    actResult.dat = tempCalcDynMF(headdata, 37, 24, 720, 360, equillibriumT, ...
        fullfile(directories(i).name,[directories(i).name '_NIFTI_1'],'Coregistered', 'Normalized_to_rest', 'Output_18-Sep-2011', 'm.mat'), ...
        fullfile(directories(i).name,[directories(i).name '_NIFTI_1'],'Coregistered', 'Normalized_to_rest', 'Output_18-Sep-2011', 'f.mat'), ...
        4, mask);
    % Store the parameters used for the calculations for reference in the future
    [c lmax] = max(actResult.dat(:));
    [likelymax x y z] = ind2sub(size(actResult.dat),lmax);
    actResult.likelymaxslice = round(likelymax/2);
    actResult.bloodT = 37;
    actResult.airT = 24;
    actResult.tmax = 360;
    actResult.stepf = 2;
    actResult.savestepf = 4;
    actResult.metabandflowdata = 'From Dataset';
    save(fullfile(directories(i).name,[directories(i).name '_NIFTI_1'],'Coregistered', 'Normalized_to_rest', 'Output_18-Sep-2011','tt_act_res.mat'), 'actResult');
    old = cd([directories(i).name,filesep,[directories(i).name '_NIFTI_1'],filesep,'Coregistered', filesep,'Normalized_to_rest', filesep,'Output_18-Sep-2011']);
    writeT_to_nii(actResult, equillibriumT, exp_nii);
    cd(old)
    clear actResult
    % Part II
    disp([int2str(i) '-2 started'])
    actResult.dat = tempCalcDynMF(headdata, 37, 24, 720, 360, equillibriumT, ...
        fullfile(directories(i).name,[directories(i).name '_NIFTI_2'],'Coregistered', 'Normalized_to_rest', 'Output_18-Sep-2011', 'm.mat'), ...
        fullfile(directories(i).name,[directories(i).name '_NIFTI_2'],'Coregistered', 'Normalized_to_rest', 'Output_18-Sep-2011', 'f.mat'), ...
        4, mask);
    [c lmax] = max(actResult.dat(:));
    [likelymax x y z] = ind2sub(size(actResult.dat),lmax);
    actResult.likelymaxslice = round(likelymax/2);
    actResult.bloodT = 37;
    actResult.airT = 24;
    actResult.tmax = 360;
    actResult.stepf = 2;
    actResult.savestepf = 4;
    actResult.metabandflowdata = 'From Dataset';
    save(fullfile(directories(i).name,[directories(i).name '_NIFTI_2'],'Coregistered', 'Normalized_to_rest', 'Output_18-Sep-2011','tt_act_res.mat'), 'actResult');
    old = cd([directories(i).name,filesep,[directories(i).name '_NIFTI_2'],filesep,'Coregistered', filesep,'Normalized_to_rest', filesep,'Output_18-Sep-2011']);
    writeT_to_nii(actResult, equillibriumT, exp_nii);
    cd(old)
    clear actResult
    disp([int2str(i) ' finished in ' num2str(toc)])
end