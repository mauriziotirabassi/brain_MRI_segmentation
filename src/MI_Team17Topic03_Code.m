%% DATA EXTRACTION
clear variables; close all; clc

% Setting the path to the resources relatively to the running OS
if isunix
    path = '/Users/mattiapezzano/Documents/GitHub/proj-mi-2023/src';
else
    path = '';
end

load('MRIdata.mat');
imshow(imread('memref.jpg'))

%% DATA PREPROCESSING
% Removal of singleton dimensions
brian = squeeze(vol);

%% VOLUME VISUALIZATION
% https://www.mathworks.com/help/images/exploring-slices-from-a-3-dimensional-mri-data-set.html

% % Displaying the MRI volume across the three planes
% figure
% orthosliceViewer(brian)
% 
% %% 
% % Displaying single planes
% figure, sliceViewer(brian, SliceDirection="X"), title("Coronal Plane")
% figure, sliceViewer(brian, SliceDirection="Y"), title("Sagittal Plane")
% figure, sliceViewer(brian, SliceDirection="Z"), title("Axial Plane")
% 
% %% 
% % Slice planes
% sx = 1; sy = 1; sz = 2.5;
% A = [sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1];
% tform = affinetform3d(A);
% brian2 = volshow(brian,renderingStyle="SlicePlanes",Transformation=tform);
% 
% %%
% % Volume rendering
% volshow(brian, RenderingStyle="VolumeRendering")

% Volume viewer
volumeViewer(brian)

%% SEGMENTATION OF SLICE 135
% Extracting the 135th sagittal slice
sliceNum = 135;
im = imrotate(squeeze(brian(sliceNum, :, :)), 90);

% Displaying the 135th sagittal slice with its histogram
figure
subplot(2, 1, 1), imshow(im, [], 'InitialMagnification', 'fit')
colorbar, drawnow
subplot(2, 1, 2), imhist(im); ylim([0 1000]), drawnow;

%% FILTERING
% TODO: Smooth out the image such that the tumor appears as a homogeneous mass

%% THRESHOLDING

thr = mean(im); % TODO: Verify optimal threshold
im_thr = imresize(im, [size(im, 1) size(im, 2)]); % TODO: Find a better way to copy images
for i = 1:1:size(im, 1)
    for j = 1:1:size(im, 2)
        if im(i, j) > thr
            im_thr(i, j) = 1;
        else
            im_thr(i, j) = 0;
        end
    end
end
figure, imshow(im_thr, [], 'InitialMagnification', 'fit')

%% SEGMENTATION

% Label connected components in the image
label = bwlabel(im_thr);

% Defining density and area of regions of the image
stats = regionprops(logical(im_thr), 'Solidity', 'Area');
density = [stats.Solidity]; area = [stats.Area];

% Isolating relatively dense areas
denseArea = density > 0.5; % TODO: Verify optimal density for tumor recognition

% Isolating the tumor as the largest area relatively to density threshold
maxArea = max(area(denseArea));

% Getting the label corresponding to the tumor area
tumorLabel = find(area == maxArea);

% Extracting the tumor out of the labeled image
tumor = ismember(label, tumorLabel);

% Displaying the segmentation
figure, imshow(tumor, [], 'InitialMagnification', 'fit')
title('Tumor Segmentation')

%% CROSS-SECTIONAL AREA
disp(['Cross-sectional area: ' int2str(maxArea)])

