function [  ] = avg_NII_normalize( varargin )
%  Uses the resting-state image calculated using
%  avg_NII_rest() to normalize the rest of the data

% If no inputs are given, the "open file..." UI will
% prompt for the required information.

%% Setup
switch length(varargin)
    case 0
        fold_name = uigetdir('Directory Containing Data');
        if ~fold_name  % Cancel Button
            return
        end
        
        [rest_file rest_path rest_index]= uigetfile('*.nii','Resting State NIFTI File');
        switch rest_index
            case 0
                return
            case 1
                rest_dat = load_nii(fullfile(rest_path,rest_file));
                rest_dat = double(rest_dat.img);
            otherwise
                error('An error has occured loading the resting state data')
        end
        
        [mask_file mask_path mask_index] = uigetfile('*.nii','Mask');
        switch mask_index
            case 0
                return
            case 1
                mask_dat = load_nii(fullfile(mask_path, mask_file));
                mask_dat = logical(mask_dat.img);
                if max(size(mask_dat) ~= size(rest_dat))
                    error('The Mask and Resting State files must have the same size')
                end
            otherwise
                error('An error has occured loading the resting state data')
        end
    case 1
        fold_name = varargin{1};
        [rest_file rest_path rest_index]= uigetfile('*.nii','Resting State NIFTI File');
        switch rest_index
            case 0
                return
            case 1
                rest_dat = load_nii(fullfile(rest_path,rest_file));
                rest_dat = double(rest_dat.img);
            otherwise
                error('An error has occured loading the resting state data')
        end
    case 2
        fold_name = varargin{1};
        rest_dat = loadNII(varargin{2});
        [mask_file mask_path mask_index] = uigetfile('*.nii','Mask');
        switch mask_index
            case 0
                return
            case 1
                mask_dat = load_nii(fullfile(mask_path, mask_file));
                mask_dat = logical(mask_dat.img);
                if max(size(mask_dat) ~= size(rest_dat))
                    error('The Mask and Resting State files must have the same size')
                end
            otherwise
                error('An error has occured loading the resting state data')
        end
    case 3
        fold_name = varargin{1};
        rest_dat = loadNII(varargin{2});
        mask_dat = loadNII(varargin{3});
    otherwise
        return
end

%  Go to the folder containing the files
oldfold = cd(fold_name);
file_list = dir('*.nii');
file_count = length(file_list);

%  Make a directoy to save the normalized data to
saveDir = 'Normalized_to_rest';
if ~isdir(saveDir)
    mkdir(saveDir);
end

statusbar = waitbar(0,'Initializing');

% for each file: load it, devide by the rest state and save it
for i=1:file_count
    try
        waitbar(i/file_count,statusbar,[fold_name sprintf('%d%%',round((i/file_count)*100))]);
    catch err
        return
    end
    [file_path file_name file_ext] = fileparts(file_list(i).name);
    file_hold = load_nii(file_list(i).name);
    file_hold.img = double(file_hold.img)./rest_dat - 1;
    file_hold.img(~mask_dat) = 0;             % set everything outside the mask to 0
    file_hold.img(isnan(file_hold.img)) = 0;  % set all NaN's to 0
    file_hold.img(isinf(file_hold.img)) = 0;  % set all inf's to 0
    file_hold.img(file_hold.img == -1) = 0;   % correct these for voxels that are giving me problems
    file_hold.hdr.dime.datatype = 16;  % set the datatype to single
    file_hold.hdr.dime.bitpix = 16;
    save_nii(file_hold,fullfile(saveDir,[file_name '_rn' file_ext]))
end

close(statusbar)
cd(oldfold)

end

