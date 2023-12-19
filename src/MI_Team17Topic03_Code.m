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

wnd = [139.5 21.5 38 29];

%% VOLUME VISUALIZATION

orthosliceViewer(vol)
% volumeViewer(vol)

%% SEGMENTATION OF SINGLE SAGITTAL SLICE k

k = 135; % Slice number
slice_k = imrotate(squeeze(vol(k, :, :)), 90); % Slice extraction

% Filtering (Testing noise resistance)
% TODO: Create separate section
figure, 
subplot(2, 2, 1), imshow(imgaussfilt(slice_k)), title('Original')

noise = imnoise(slice_k, 'salt & pepper', 0.05);
subplot(2, 2, 2), imshow(noise), title('Noise')

filtered = medfilt2(noise);
filtered = imgaussfilt(filtered);
subplot(2, 2, 3), imshow(filtered), title('Filtered')

% Preprocessing
% TODO: Calling broken down substep functions
%       (right now the whole pipeline is inside seg())
%       gamma point operator (clean around the tumor)

[seg_k, area_k] = seg(slice_k, wnd); % Segmentation
[seg_filt, area_filt] = seg(filtered, wnd);

% Displaying the slice with its cross-sectional area
figure, 
subplot(2, 1, 1), imshow(seg_k, [], 'InitialMagnification', 'fit')
title(['Sagittal Slice ' int2str(k) ' Cross-Sectional Area: ' int2str(area_k)])

subplot(2, 1, 2), imshow(seg_filt, [], 'InitialMagnification', 'fit')
title(['Sagittal Slice ' int2str(k) ' Filtered Area: ' int2str(area_filt)])

%% SEGMENTATION OF WHOLE SAGITTAL VOLUME

% Rotating the volume such that montage outputs sagittal slices
% sag_brian = imrotate3(brian, 180, [0 1 1]);
% figure, montage(sag_brian)
figure
for k = 106:1:145 % size(vol, 1) % Empirical
    im_k = imrotate(squeeze(vol(k, :, :)), 90); % Extracting the sagittal slice
    [tumor, area] = seg(im_k, wnd); 

    % TODO: Find a better way to remove irrelevant tissue

    % Skipping slices without dense material
    if area < 0 % Empirical and incorrect!!
        continue
    else
        %TODO: Understand why imshow() parameters don't work
        imshow(tumor) % Displaying segmentation
        % disp(['Cross-sectional area of slice ' num2str(k) ' is ' num2str(area)]) % Displaying cross-sectional area
    end
end

% TODO: Create montage of sagittal tumor segmentations

%% SEGMENTATION OF WHOLE AXIAL VOLUME
% TODO: Check boundaries from orthosliceviewer

%% SEGMENTATION WITH NOISE
% TODO: Move here

%% SEGMENTATION FUNCTION

% TODO: Define thresholds as params if we decide on them being variable
% TODO: Break down into substep functions
function [im_out, area] = seg(im_in, rect)
    crp_im = imcrop(im_in, [], rect); % Cropping out the ROI

    % TODO: Different meaning between images
    bw_im = im2double(crp_im); % Standardizing the gray scale
    bw_im = bw_im > 0.588 & bw_im < 0.862; % Thresholding

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
