function [ head_out ] = build_skin( head_in )
% build_skin() Creates a layer of skin around the head
% 
% This will check all voxels that were previously labeled
% as soft tissue and checks if it has a neighbor which is air.
% If so, then it is reassigned as skin.

if ndims(head_in)==4
    head_in = head_in(:,:,:,1);
end

% Git a list of all voxels labeled as muscle
muscle_voxels = find(head_in==13);

% Go through each of them and check for neighboring air voxels
for i=1:length(muscle_voxels)
   [x y z] = ind2sub(size(head_in), muscle_voxels(i));
   % makes sure we're not at a voxel at the boundry of the dataset
   if (x~=1) && (x~=size(head_in,1)) && (y~=1) && (y~=size(head_in,2)) && (z~=1) && (z~=size(head_in,3))
     % Looks for neighboring voxels that are air
     if ((head_in(x+1,y,z)==1) || (head_in(x-1,y,z)==1) || (head_in(x,y+1,z)==1) || (head_in(x,y-1,z)==1) || (head_in(x,y,z+1)==1) || (head_in(x,y,z-1)==1))
         head_in(x,y,z) = 14;
     end
   end
end

head_out = repair_headdata(head_in);

end