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

%% VOLUME VISUALIZATION

% orthosliceViewer(brian)
volumeViewer(vol)

%% SEGMENTATION OF SINGLE SAGITTAL SLICE k

k = 135; % Slice number
slice_k = imrotate(squeeze(vol(k, :, :)), 90); % Slice extraction

% Filtering
% TODO: Filtering that could help (stuff w/ histograms?)

% Preprocessing
% TODO: Calling broken down substep functions
%       (right now the whole pipeline is inside seg())

[seg_k, area_k] = seg(slice_k); % Segmentation

% Displaying the slice with its cross-sectional area
figure, imshow(seg_k, [], 'InitialMagnification', 'fit')
title(['Sagittal Slice ' int2str(k) ' Cross-Sectional Area: ' int2str(area_k)])

%% SEGMENTATION OF WHOLE SAGITTAL VOLUME

% Rotating the volume such that montage outputs sagittal slices
% sag_brian = imrotate3(brian, 180, [0 1 1]);
% figure, montage(sag_brian)

figure
for k = 1:1:size(vol, 1)
    im_k = imrotate(squeeze(vol(k, :, :)), 90); % Extracting the sagittal slice
    [tumor, area] = seg(im_k); 

    % TODO: Find a better way to remove irrelevant tissue

    % Skipping slices without dense material
    if area < 100 % Empirical and incorrect!!
        continue
    else
        %TODO: Understand why imshow() parameters don't work
        imshow(tumor) % Displaying segmentation
        % disp(['Cross-sectional area of slice ' num2str(k) ' is ' num2str(area)]) % Displaying cross-sectional area
    end
end

% TODO: Create montage of sagittal tumor segmentations

%% SEGMENTATION OF WHOLE AXIAL VOLUME
% TODO: tutto

%% SEGMENTATION WITH NOISE
% TODO: tutto

%% SEGMENTATION FUNCTION

% TODO: Define thresholds as params if we decide on them being variable
% TODO: Break down into substep functions
function [im_out, area] = seg(im_in)

    bw_im = im_in > 150 & im_in < 220; % Thresholding
    label = bwlabel(bw_im); % Labeling continuous regions

    % Stats of continuous regions
    stats = regionprops(logical(bw_im), 'Solidity', 'Area');
    density = [stats.Solidity]; area = [stats.Area];

    % Supposing the tumor is the largest dense area
    denseArea = density > 0.6; % Empirical
    maxArea = max(area(denseArea));
    tumorLabel = find(area == maxArea);
    im_out = ismember(label, tumorLabel);

    area = maxArea; % Cross-sectional area of the supposed tumor

end
