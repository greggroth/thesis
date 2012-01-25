function [ head_out ] = build_skin( head_in )
%build_skin Summary of this function goes here
%   Detailed explanation goes here

if ndims(head_in)==4
    head_in = head_in(:,:,:,1);
end

muscle_voxels = find(head_in==13);

for i=1:length(muscle_voxels)
   [x y z] = ind2sub(size(head_in), muscle_voxels(i));
   % if a muscle voxel borders any air voxels, it's set to skin
   if (x~=1)&&(x~=size(head_in,1))&&(y~=1)&&(y~=size(head_in,2))&&(z~=1)&&(z~=size(head_in,3))
       if ((head_in(x+1,y,z)==1)||(head_in(x-1,y,z)==1)||(head_in(x,y+1,z)==1)||(head_in(x,y-1,z)==1)||(head_in(x,y,z+1)==1)||(head_in(x,y,z-1)==1))
           head_in(x,y,z) = 14;
       end
   end
end

head_out = repair_headdata(head_in);


end

