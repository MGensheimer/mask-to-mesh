% Convert binary mask in DICOM format to smoothed STL mesh
% Created by Michael Gensheimer, Dec. 2016. Contact: michael.gensheimer@gmail.com
% Requires the "Smooth Triangulated Mesh" package from Dirk-Jan Kroon:
% https://www.mathworks.com/matlabcentral/fileexchange/26710-smooth-triangulated-mesh
%
% Input: directory with Dicom images that create a binary mask. Image value
% should be 1 inside the mask, 0 outside the mask.
% Output: STL polygon file

clear all;
resampVoxelSize=2; % resampled isotropic resolution in mm. 1 for 3D printed dental stent, 2-3 for foam cradle
smoothAmt=0.5; % can be customized
dicomDir=''; % replace with directory with Dicom files
outFile='mesh.stl';

padding=20;
cd(dicomDir);
dicomlist=dir(fullfile(pwd,'*.dcm'));
info=dicominfo(fullfile(pwd,dicomlist(1).name));
voxelSize=zeros(1,3);
voxelSize(1:2)=info.PixelSpacing;
outPlaneVec=cross(info.ImageOrientationPatient(1:3),info.ImageOrientationPatient(4:6));
clear instances; clear sliceSpace;
pos=zeros(numel(dicomlist),3);
for i=1:numel(dicomlist)
    info=dicominfo(fullfile(pwd,dicomlist(i).name));
    instances(i)=info.InstanceNumber;
    pos(instances(i),:)=info.ImagePositionPatient;
end
maskUncropped=zeros(info.Width,info.Height,numel(dicomlist));
for i=1:numel(dicomlist)
    maskUncropped(:,:,instances(i))=dicomread(fullfile(pwd,dicomlist(i).name));
    if i>1
        sliceDiff=pos(i,:)-pos(i-1,:);
        sliceSpace(i-1)=abs(dot(sliceDiff,outPlaneVec./norm(outPlaneVec)));
    end
end
if (max(sliceSpace)-min(sliceSpace))>0.1
    disp 'Slice spacing varies. Aborting.'
    keyboard;
end
voxelSize(3)=min(sliceSpace);
[dimx, dimy, dimz]=size(maskUncropped);
if max(maskUncropped(:))==0 && min(maskUncropped(:))==-1 %MIM sometimes saves masks with negative values
   maskUncropped=maskUncropped+1; 
end
if max(maskUncropped(:))>1 %change to 1's and 0's
   maskUncropped=maskUncropped+1-max(maskUncropped(:)); 
end
whereTumor=find(maskUncropped);
[whereTumorX,whereTumorY,whereTumorZ]=ind2sub([dimx dimy dimz],whereTumor);
minx=min(whereTumorX);
maxx=max(whereTumorX);
miny=min(whereTumorY);
maxy=max(whereTumorY);
minz=min(whereTumorZ);
maxz=max(whereTumorZ);

maskCropped=zeros(dimx+padding*2,dimy+padding*2,dimz+padding*2);
maskCropped(padding+1:end-padding,padding+1:end-padding,padding+1:end-padding) = maskUncropped(:,:,:);

[dimx, dimy, dimz]=size(maskCropped);

[x y z]=ndgrid((1:dimx)*voxelSize(1),(1:dimy)*voxelSize(2),(1:dimz)*voxelSize(3));
[xq yq zq]=ndgrid(1:resampVoxelSize:dimx*voxelSize(1),1:resampVoxelSize:dimy*voxelSize(2),1:resampVoxelSize:dimz*voxelSize(3));
resampMask = interpn(x,y,z,maskCropped,xq,yq,zq);
fv=isosurface(xq,yq,zq,resampMask,0.5);
fv2=smoothpatch(fv,1,5,smoothAmt);
stlwrite(outFile,fv2);