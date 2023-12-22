%% DATA EXTRACTION
clear variables; close all; clc

% Setting the path to the resources relatively to the running OS
if isunix
    path = '/Users/mattiapezzano/Documents/GitHub/proj-mi-2023/src';
else
    path = '';
end

load('MRIdata.mat');
% imshow(imread('memref.jpg'))

%% VOLUME VISUALIZATION

% Visualizing the volume to approximately assess:
% 2. The ranges of the slices of interest for the segmentation
% 3. The slices containing the biggest tumor extention
orthosliceViewer(vol);
% volumeViewer(vol)
sagittal_range = [105 145]; axial_range = [63 91]; % Empirical
big_sagittal = 123; big_axial = 78; % Empirical

% Deifning the dimensions of the volumes
sagittal_dim = [size(vol, 2) size(vol, 3) range(sagittal_range)];
axial_dim = [size(vol, 1) size(vol, 2) range(axial_range)];

% Using imcrop to manually extract the ROI coordinates
% imcrop(sagittal(vol, big_sagittal));
% imcrop(axial(vol, big_axial));
sagittal_window = [140 23 39 29]; axial_window = [138 107 42 39]; % Empirical

%% SEGMENTATION OF SINGLE SAGITTAL SLICE k

k = 120; % Slice number
slice_k = sagittal(vol, k); % Extracting the slice
[tumor_k, area_k] = segmentation2(slice_k, sagittal_window); % Segmenting the tumor

% Overlaying the segmented tumor on the original image
overlay_k = overlay(slice_k, tumor_k, sagittal_window);

% Displaying the segmentation
imshow(overlay_k, [], 'InitialMagnification', 'fit')
title(['Sagittal Slice: ' int2str(k) ' - Cross-Sectional Area: ' int2str(area_k)])
drawnow

%% SEGMENTATION OF WHOLE SAGITTAL VOLUME

% Rotating the volume such that montage outputs sagittal slices
% sag_brian = imrotate3(brian, 180, [0 1 1]);
% figure, montage(sag_brian)

sagittal_volume = ones(sagittal_dim);
for k = sagittal_range(1):1:sagittal_range(2)
    slice = sagittal(vol, k);
    [tumor, area] = segmentation2(slice, sagittal_window); 
    segmentation = overlay(slice, tumor, sagittal_window);
    imshow(segmentation, [], 'InitialMagnification', 'fit'), drawnow

    % TODO: Create a montage instead of an animation (RGB problem)
    % TODO: Understand how to display the cross-sectional area values
    % sagittal_volume(:, :, k) = segmentation;
end

%% SEGMENTATION OF WHOLE AXIAL VOLUME

axial_volume = ones(axial_dim);
for k = axial_range(1):1:axial_range(2)
    slice = axial(vol, k);
    [tumor, area] = segmentation2(slice, axial_window); 
    segmentation = overlay(slice, tumor, axial_window);
    % imshow(segmentation, [], 'InitialMagnification', 'fit'), drawnow

    % TODO: Create a montage instead of an animation (RGB problem)
    % TODO: Understand how to display the cross-sectional area values
    % axial_volume(:, :, k) = segmentation;
end

%% SEGMENTATION WITH NOISE

% TODO: Clean

% Filtering (Testing noise resistance)
% TODO: Create separate section
figure, 
subplot(2, 2, 1), imshow(imgaussfilt(sagittal_k)), title('Original')

noise = imnoise(sagittal_k, 'salt & pepper', 0.05);
subplot(2, 2, 2), imshow(noise), title('Noise')

filtered = medfilt2(noise);
filtered = imgaussfilt(filtered);
subplot(2, 2, 3), imshow(filtered), title('Filtered')

% Preprocessing
% TODO: Calling broken down substep functions
%       (right now the whole pipeline is inside seg())
%       gamma point operator (clean around the tumor)

[seg_k, area_k] = segmentation2(sagittal_k, sagittal_window); % Segmentation
[seg_filt, area_filt] = segmentation2(filtered, sagittal_window);

% Displaying the slice with its cross-sectional area
figure, 
subplot(2, 1, 1), imshow(seg_k, [], 'InitialMagnification', 'fit')
title(['Sagittal Slice ' int2str(k) ' Cross-Sectional Area: ' int2str(area_k)])

subplot(2, 1, 2), imshow(seg_filt, [], 'InitialMagnification', 'fit')
title(['Sagittal Slice ' int2str(k) ' Filtered Area: ' int2str(area_filt)])

%% FUNCTIONS

% Extracts a specific sagittal slice from a volume
function [im_out] = sagittal(volume, slice_num)
    im_out = imrotate(squeeze(volume(slice_num, :, :)), 90);
end

% Extracts a specific axial slice from a volume
function [im_out] = axial(volume, slice_num)
    im_out = squeeze(volume(:, :, slice_num));
end

% TODO: Define thresholds as params if we decide on them being variable
% TODO: Break down into substep functions
function [im_out, area] = segmentation2(im_in, rect)
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
    im_out = imfill(ismember(label, tumorLabel), 'holes');

    area = maxArea; % Cross-sectional area of the supposed tumor

end

% Overlays the second image on top of the first one in the given rectangle
function [im_out] = overlay(im_in1, im_in2, rect)
    pad_im_2 = padarray(im_in2, [rect(2) - 1, rect(1) - 1], 0, 'pre');
    im_out = imfuse(im_in1, pad_im_2);
end
