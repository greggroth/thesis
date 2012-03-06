function [ ] = TempPlot( head, tempdata, highlightRegion, slice, equil,threshold,point)
%TempPlot Plot data from tempCalc() or BulkImportNII()
%   INPUT TempPlot(structuredata)
%         TempPlot(structuredata,temperaturedata)
%         TempPlot(structuredata,temperaturedata,highlightRegion)
%         TempPlot(structuredata,temperaturedata,highlightRegion,slice)
%         TempPlot(structuredata,temperaturedata,highlightRegion,slice,EquillibriumData)
%         
%   This function with determine which type of data it is and then
% plot it appropiately.  
%
%   equil - Equillibrium state data
%   threshold - threshold value for being displayed as an overlay
%   REQUIRES:  SliceBrowser (http://www.mathworks.com/matlabcentral/fileexchange/20604)
%%  Error checking and data restructuring where necessary
if ndims(head) == 4
    head = head(:,:,:,1);
elseif ndims(head) ~= 3
    error('Input ''head'' must have either 3 or 4 dimensions');
end

if nargin > 1
    if ndims(tempdata) == 3  % should only happen when comparing two equilibrium datasets
      temp = tempdata;
      tempdata = zeros([1 size(temp)]);
      tempdata(1,:,:,:) = temp;
    elseif ndims(tempdata) ~= 4
      error('Input ''tempdata'' must have either 3 or 4 dimensions');
    end
    tempdataShort = squeeze(tempdata(end,:,:,:));
end

if nargin > 2
    if ndims(highlightRegion) ~= 3
      error('Input ''highlightRegion'' must have 3 dimensions');
    end
    if size(highlightRegion) ~= size(head)
      error('Input ''highlightRegion'' must be of the same size as ''head''');
    end
    tempdataShort = squeeze(tempdata(end,:,:,:));
end

if nargin > 3
    if slice > size(tempdata,1)
      error('Input ''slice'' must be less or equal to the length of the first dimension of ''tempdata''');
    end
    tempdataShort = squeeze(tempdata(slice,:,:,:));
end

if nargin > 4
    if ndims(equil) == 3
        eq = equil;
    elseif ndims(equil) == 4
        eq = squeeze(equil(1,:,:,:));
    else
        error('Input ''equil'' must have either 3 or 4 dimensions');
    end
    clear 'equil';
end

%%  Pick how to format the call of SliceBrowser()
switch nargin
    case 1
    SliceBrowser(head,1,head);
    colormap(gray);
    case 2
    SliceBrowser(tempdataShort,tempdata,head);
    case 3
    SliceBrowser(tempdataShort,tempdata,head,highlightRegion);
    case 4
    SliceBrowser(tempdataShort,tempdata,head,highlightRegion);
    case 5
    SliceBrowser(tempdataShort-eq,tempdata,head,highlightRegion);
    case 6
    SliceBrowserOverlay(tempdataShort-eq,tempdata,head,highlightRegion,threshold);
    case 7
    imgoverlay(head,tempdataShort-eq,point,threshold)
end  

end

