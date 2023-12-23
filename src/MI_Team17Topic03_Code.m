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
% 1. The ranges of the slices of interest for the segmentation
% 2. The slices containing the biggest tumor extention
orthosliceViewer(vol);
% volumeViewer(vol)
sagittal_range = [105 145]; axial_range = [63 91]; % coronal_range = []; % Empirical
big_sagittal = 123; big_axial = 78; % big_coronal = 158; % Empirical

% % Deifning the dimensions of the volumes
% sagittal_dim = [size(vol, 2) size(vol, 3) range(sagittal_range)];
% axial_dim = [size(vol, 1) size(vol, 2) range(axial_range)];

% Using imcrop to manually extract the ROI coordinates
% imcrop(slice(selectvol(vol, 'sagittal'), big_sagittal));
% imcrop(slice(selectvol(vol, 'axial'), big_axial));
% imcrop(slice(selectvol(vol, 'coronal'), big_coronal));
sagittal_window = [140 23 39 29]; axial_window = [138 107 42 39]; % Empirical

%% SEGMENTATION OF SINGLE SAGITTAL SLICE k

% TODO: Review segmentation parameters. Missing pieces for slice 135. Not
% good for the first thing the professor sees.

% Segmenting the tumor
k = 135; % Slice number
slice_k = slice(selectvol(vol, 'sagittal'), k); % Extracting the slice
[tumor_k, area_k] = segmentation2(slice_k, sagittal_window);

% Displaying the segmentation
imshow(tumor_k, [], 'InitialMagnification', 'fit')
title(['Sagittal Slice: ' int2str(k) ' - Cross-Sectional Area: ' int2str(area_k)])
drawnow

%% SEGMENTATION OF THE WHOLE VOLUME

clc, valid_type = {'axial', 'sagittal', 'coronal'};
while 1
    % Requesting user input (choice of type of slices)
    volume_type = input('Write the plane of interest: ', 's');
    if strcmp(volume_type, 'quit')
        clc
        break
    else
        while ~any(strcmp(volume_type, valid_type))
            if strcmp(volume_type, 'quit')
                clc
                break
            else
                volume_type = input('Write the plane of interest: ', 's');
            end
        end 
    end
    % TODO: Noise option

    % TODO: Improve pipeline
    % Segmenting the volume
    [sel_volume, sel_range, sel_window] = selectvol(vol, volume_type);
    segmentation3(sel_volume, sel_range, sel_window)
    close all, clc
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

%%
vol_t = selectvol(vol, 'axial');
subplot(1, 2, 1), montage(vol_t)
vol_n = addnoise(vol_t, 'salt & pepper', 0.05);
subplot(1, 2, 2), montage(vol_n)


%%
h = 100;
c_vol = selectvol(vol, 'sagittal');
vol_noise = zeros(size(c_vol));
vol_noise(:, :, h) = imnoise(slice(c_vol, h), 'salt & pepper', 0.10);
subplot(1, 2, 1), imshow(slice(c_vol, h))
subplot(1, 2, 2), imshow(squeeze(vol_noise(:, :, h)))
%% FUNCTIONS

% TODO: Understand if one can pass ranges and windows
% Takes as input the axial volume and returns a volume whose third dimension 
% are the slices in the plane of interest
function [vol_out, range_out, window_out] = selectvol(volume, type)
    switch type
        case 'axial'
            vol_out = volume;
            range_out = [63 91]; window_out = [138 107 42 39];
        case 'sagittal'
            vol_out = flip(flip(imrotate3(volume, 180, [0 -1 1]), 2), 3);
            range_out = [105 155]; window_out = [140 23 39 29];
        % case 'coronal'
        %     vol_out = imrotate3(volume, 180, [1 0 1]);
        %     range_out = [135 180]; window_out = [];
    end
end

% Extracts a specific slice from a volume (always the third dimension)
function [im_out] = slice(volume, slice_num)
    im_out = squeeze(volume(:, :, slice_num));
end

% Overlays the second image on top of the first one in the given rectangle
function [im_out] = overlay(im_in1, im_in2, rect)
    pad_im_2 = padarray(im_in2, [rect(2) - 1, rect(1) - 1], 0, 'pre');
    im_out = imfuse(im_in1, pad_im_2);
end

% Performs preprocessing and segmentation on a slice
function [im_out, area] = segmentation2(im_in, rect)
    crp_im = imcrop(im_in, [], rect); % Cropping out the ROI

    % TODO: Different meaning between images
    bw_im = im2double(crp_im); % Standardizing the gray scale
    bw_im = bw_im > 0.588 & bw_im < 0.862; % Thresholding % Empirical

    label = bwlabel(bw_im); % Labeling continuous regions

    % Stats of continuous regions
    stats = regionprops(logical(bw_im), 'Solidity', 'Area');
    density = [stats.Solidity]; area = [stats.Area];

    % Supposing the tumor is the largest dense area
    denseArea = density > 0.6; % Areas with empirical density
    maxArea = max(area(denseArea)); % Largest area with empirical density
    tumorLabel = find(area == maxArea); % Label of the supposed tumor
    tumor = ismember(label, tumorLabel); % Extracting the image from the label
    filled_tumor = imfill(tumor, 'holes'); % Correcting the image
    
    % Overlaying the segmented tumor on the original image
    im_out = overlay(im_in, filled_tumor, rect);
    
    area = maxArea; % Cross-sectional area of the supposed tumor

end

% Performs preprocessing and segmentation on a volume
function [] = segmentation3(volume, range, window)
    for k = range(1):1:range(2)
        slice_k = slice(volume, k);
        [tumor, area] = segmentation2(slice_k, window); 

        % TODO: Decide function behavior
        imshow(tumor, [], 'InitialMagnification', 'fit'), drawnow
    
        % TODO: Create a montage instead of an animation (RGB problem)
        % TODO: Understand how to display the cross-sectional area values
        % axial_volume(:, :, k) = tumor;
    end
end

% Adds noise to a volume
function [vol_out] = addnoise(vol_in, noise, degree)
    vol_out = zeros(size(vol_in));
    for k = 1:1:size(vol_in, 3)
        vol_out(:, :, k) = imnoise(slice(vol_in, k), noise, degree);
    end
end
