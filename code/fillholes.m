function [ out_head ] = fillHoles( in1,in2,in3,in4,in5,headin)
% fillHoles() checks for misassigned voxels
% 
% Solves an issue where a voxel with two equally probable tissue
% types resulted in being assigned as air.  This checks for air 
% voxels that are surrounded by tissue and decides a tissue it
% it would be best suited as

head = squeeze(headin(:,:,:,1));  %  I only need the tissue indices so this makes things easier down the line

%%  Data Storage
QmSTORE = [0 0 26.1 11600 0 26.1 697 0 0 302 15575 0 697 1100 5192];
cSTORE = [1006 4600 2110 3640 3800 1300 3720 3000 4200 2300 3680 3500 3720 3150 3600];
rhoSTORE = [1.3 1057 1080 1035.5 1007 1850 1126 1076 1009 916 1035.5 1151 1041 1100 1027.4];
kSTORE = [0.026 0.51 0.65 0.534 0.5 0.65 0.527 0.4 0.594 0.25 0.565 0.4975 0.4975 .342 .503];
wSTORE = [0 1000 3 45.2 0 1.35 40 0 0 2.8 67.1 3.8 3.8 12 23.7];

%%  Get locations of holes
%   Where two tissue types have the same probability

idx1 = (in1==in2 | in1 == in3 | in1==in4 | in1==in5) & logical(in1);
idx2 = (in1==in2 | in2 == in3 | in2==in4 | in2==in5) & logical(in2);
idx3 = (in1==in3 | in2 == in3 | in3==in4 | in3==in5) & logical(in3);
idx4 = (in1==in4 | in2 == in4 | in3==in4 | in4==in5) & logical(in4);
idx5 = (in1==in5 | in2 == in5 | in3==in5 | in4==in5) & logical(in5);
%  This array will have a zero anywhere there were two or more common
%  elements between any of the five arrays.  
idx = idx1|idx2|idx3|idx4|idx5;   

[xmax ymax zmax] = size(in1)
[x y z] = ind2sub(size(in1),find(idx));  %  get x, y and z coordinates of the holes

for i = 1:length(x)  %  go to each hole and do work
    if (x(i)~=1)&&(y(i)~=1)&&(z(i)~=1)&&(x(i)~=xmax)&&(y(i)~=ymax)&&(z(i)~=zmax)&&(headin(x(i),y(i),z(i),1)==1)  %  keeps away from the edge and only looks at voxels that were assigned air
        [commonesttissue nouse secondbest] = mode([head(x(i)+1,y(i),z(i)) head(x(i)-1,y(i),z(i)) head(x(i),y(i)+1,z(i)) head(x(i),y(i)-1,z(i)) head(x(i),y(i),z(i)+1) head(x(i),y(i),z(i)-1)]);
        if commonesttissue == 1 && length(secondbest{1})>=2  % if air and something else are equally common, it'll choose air.  This forces it to pick the tissue if possible.
            commonesttissue = secondbest{1}(2);
        end
        headin(x(i),y(i),z(i),:) = [commonesttissue 0 QmSTORE(commonesttissue) cSTORE(commonesttissue) rhoSTORE(commonesttissue) kSTORE(commonesttissue) wSTORE(commonesttissue)];
    end
end

out_head = headin;

end

