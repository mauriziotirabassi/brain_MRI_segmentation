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

% Displaying the MRI volume across the three planes
figure
orthosliceViewer(brian)

%% 
% Displaying single planes
figure, sliceViewer(brian, SliceDirection="X"), title("Coronal Plane")
figure, sliceViewer(brian, SliceDirection="Y"), title("Sagittal Plane")
figure, sliceViewer(brian, SliceDirection="Z"), title("Axial Plane")

%% 
% Slice planes
sx = 1; sy = 1; sz = 2.5;
A = [sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1];
tform = affinetform3d(A);
brian2 = volshow(brian,renderingStyle="SlicePlanes",Transformation=tform);

%%
% Volume rendering
volshow(brian, RenderingStyle="VolumeRendering")

%%
% Volume viewer
volumeViewer(brian)



