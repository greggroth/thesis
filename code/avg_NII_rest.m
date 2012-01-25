function [ ] = avg_NII_rest( varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Setup
switch length(varargin)
    case 0
        fold_name = uigetdir;
        if ~fold_name  % Cancel Button
            return
        end
    case 1
        fold_name = varargin{1};
    otherwise
end

%  Go to the folder containing the files
oldfold = cd(fold_name);
file_list = dir('*.nii');

%  We're only interested in the rest period 
%  (first and last 10 steps in this case)
file_list = file_list([1:10 170:180]);
file_count = length(file_list);

%  Cell array to store all of the datasets in.
datHolder = cell(file_count,1);

statusbar = waitbar(0,'Initializing');

for j=1:file_count
    try
        waitbar(j/file_count,statusbar,sprintf('%d%%',round((j/file_count)*100)));
    catch err
        return
    end
    fi = load_nii(file_list(j).name);
    datHolder{j} = fi.img;
end

%%  Calculate the mean
ymax = size(datHolder{1},2);
zmax = size(datHolder{1},3);
output = zeros(size(datHolder{1}));

for i=1:ymax
    try
        waitbar(i/ymax,statusbar,sprintf('%d%%',round((i/ymax)*100)));
    catch err
        return
    end
    for k=1:zmax
        excStr = cell(length(datHolder),1);
        for l=1:length(datHolder)
            excStr{l} = [',datHolder{' int2str(l) '}(:,' int2str(i) ',' int2str(k) ')'''];
        end
        comb = eval(['cat(1' cell2mat(excStr') ')']);
        output(:,i,k) = mean(comb);
    end
end

close(statusbar)

fi.img = output;
mkdir('RestState')
save_nii(fi,fullfile('RestState','RestStateAvg.nii'));

cd(oldfold)

end

