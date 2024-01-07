%% DATA EXTRACTIO
% N
clear variables; close all; clc

% Setting the path to the resources relatively to the running OS
if isunix
    path = '/Users/mattiapezzano/Documents/GitHub/proj-mi-2023/src';
else
    path = '';
end

load('MRIdata.mat');
% imshow(imread('memref.jpg'))

%% VOLUME VISUALIZATION (DATA PREPARATION)

% Visualizing the volume to approximately assess the ranges of the slices 
% of interest for the segmentation for each plane
sagittal_range = [100 145];
% axial_range = [60 90];
% coronal_range = [140 180];
orthosliceViewer(vol);

% Extracting the orthogonal planes according to the ranges (see function
% select_vol()) and permuting the volume such that the third dimension
% iterates through the slices of interest
% sg_planes = rot90(permute(squeeze(vol(100:145,:,:)),[2,3,1]));
% ax_planes = rot90(squeeze(vol(:,:,60:90)));
% cr_planes = rot90(permute(squeeze(vol(:,140:180,:)),[1,3,2]));

% For each plane identifying the slice with the biggest lesion extention and
% using it to empirically define the ROI window using imcrop()
% sagittal_window = [140 25 40 30]; % Original sagittal slice 125
% axial_window = [105 80 45 45]; % Original axial slice 79
% coronal_window = [100 20 45 30]; % Original coronal slice 158

%% 3. SEGMENTATION OF SAGITTAL SLICE 135

og_num = 135; % From reference
actual_num = 135 - sagittal_range(1); % Slice number
[sg_vol, sg_wnd] = selectvol(vol, 'sagittal'); % Selecting sagittal volume and ROI
slice_k = slice(sg_vol, actual_num); % Extracting the slice
[tumor_k, area_k] = segmentation2(slice_k, sg_wnd); % Segmenting the lesion

% Displaying the segmentation
imshow(tumor_k, [], 'InitialMagnification', 'fit')
title(['Original Sagittal Slice: ' int2str(og_num) ' - Cross-Sectional Area: ' int2str(area_k)])
drawnow

%% 4. SEGMENTATION OF THE WHOLE VOLUME (5. WITH NOISE INTRODUCTION)

% Command-line interface (CLI)
valid_volume = {'axial', 'sagittal', 'coronal', ''};
valid_noise = {'salt & pepper', 'gaussian', ''};
while 1
    % Requesting user to choose type of slice
    volume_type = input('Insert the plane of interest (enter to quit): ', 's');
    while ~any(strcmp(volume_type, valid_volume))
        volume_type = input('Insert the plane of interest: ', 's');
    end
    if strcmp(volume_type, '')
        clc, break
    end

    % Requesting user to add noise
    noise_type = input('Insert noise to add (enter if none): ', 's');
    while ~any(strcmp(noise_type, valid_noise))
        noise_type = input('Insert noise to add (enter if none): ', 's');
    end 
    if strcmp(noise_type, '')
        chosen_vol = vol;
    else
        % Requesting user noise level
        noise_level = input('Insert noise level: ', 's');
        while str2double(noise_level) < 0 | str2double(noise_level) > 1 | strcmp(noise_level, '')
            noise_level = input('Insert noise level: ', 's');
        end 

        % Supposing adding the noise is independent from slice type
        chosen_vol = addnoise(vol, noise_type, str2double(noise_level));
    end

    % Segmenting the volume
    [sel_volume, sel_window] = selectvol(chosen_vol, volume_type);
    segmentation3(sel_volume, sel_window);
    clc, close all
end

%% 5. 6. NOISE SENSITIVITY ANALYSIS
% Assessing sensitivity to different noise levels and whether the application 
% of a preemptive median filter enhances general performance

% Calculating the detected area relatively to different levels of noise
filtered = zeros(1, 100); noised = zeros(1, 100); i = 1;
for noise_level = 0.01:0.01:1
    im_noise = imnoise(slice_k, 'salt & pepper', noise_level); % Adding noise to the sample slice
    tmp_im = imcrop(im_noise, [], sg_wnd); % Cropping out the ROI
    tmp_im = im2double(tmp_im); % Standardizing the gray scale
    tmp_im = tmp_im > 0.50 & tmp_im < 0.85; % Empirical thresholding
    label = bwlabel(tmp_im); % Labeling continuous regions
    stats = regionprops(logical(tmp_im), 'Solidity', 'Area'); % Statistics of continuous regions
    density = [stats.Solidity]; area = [stats.Area];
    denseArea = density > 0.6; % Areas with empirical density
    area_noise = max(area(denseArea)); % Largest area with empirical density
    if isempty(area_noise)
        area_noise = 0;
    end
    noised(i) = area_noise; % Saving the cross-sectional area

    % Segmenting with the application of the median filter
    [~, area_filt] = segmentation2(im_noise, sg_wnd);
    if isempty(area_filt)
        area_filt = 0;
    end
    filtered(i) = area_filt; % Saving the cross-sectional area
    i = i + 1;
end

figure, plot(filtered), hold on, plot(noised), hold on
refline(0, area_k); % Cross-sectional area without noise
legend({'Filtered', 'Noise', 'Clean'}, 'Location', 'southwest')
ylabel('Cross-Sectional Area'), xlabel('Noise Percentage')

% Therefore we decide to always preemptively apply the median filter to 
% reduce eventual salt & pepper noise degradation

%%

% File 'true_segmentation.mat' was created using the volume segmenter
% application, applying Otsu thresholding and smoothing edges + removing
% some excess manually
true_seg = load('true_segmentation.mat').labels; % Ground truth

%% Sagittal Test with Noise
for noise_level = 0.1:0.1:0.8
    sg = addnoise(sg,'gaussian', noise_level);
    ground_truth = true_seg;
    assess_gamma(sg, ground_truth)
end
%no significant difference at increasing noise level

%% Axial Test
ground_truth = permute(true_seg, [3 1 2]);
assess_gamma(vol, 'axial',  ground_truth)

%% Coronal Test
ground_truth = permute(true_seg, [2 1 3]);
assess_gamma(cr, ground_truth)

%% FUNCTIONS

% Takes as input the axial volume and returns a volume whose third dimension 
% are the slices in the plane of interest
function [vol_out, window_out] = selectvol(volume, type)
    switch type
        case 'axial'
            vol_out = rot90(squeeze(volume(:, :, 60:90)));
            window_out = [105 80 45 45];
        case 'sagittal'
            vol_out = rot90(permute(squeeze(volume(100:145, :, :)), [2, 3, 1]));
            window_out = [140 25 40 30];
        case 'coronal'
            vol_out = rot90(permute(squeeze(volume(:, 140:180, :)),[1, 3, 2]));
            window_out = [105 20 45 30];
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

% 2. IMPLEMENTED WORKFLOW
% Performs preprocessing and segmentation on a slice
function [im_out, area] = segmentation2(im_in, rect)
    tmp_im = imcrop(im_in, [], rect); % Cropping out the ROI
    tmp_im = medfilt2(tmp_im); % Applying median filtering
    tmp_im = im2double(tmp_im); % Standardizing the gray scale
    tmp_im = tmp_im > 0.50 & tmp_im < 0.85; % Empirical thresholding
    label = bwlabel(tmp_im); % Labeling continuous regions

    % Statistics of continuous regions
    stats = regionprops(logical(tmp_im), 'Solidity', 'Area');
    density = [stats.Solidity]; area = [stats.Area];

    % % Supposing the tumor is the largest dense area
    denseArea = density > 0.6; % Areas with empirical density
    maxArea = max(area(denseArea)); % Largest area with empirical density
    lesionLabel = find(area == maxArea); % Label of the supposed tumor
    lesion = ismember(label, lesionLabel); % Extracting the image from the label
    filled_lesion = imfill(lesion, 'holes'); % Correcting the image

    % Overlaying the segmented tumor on the original image
    im_out = overlay(im_in, filled_lesion, rect);

    area = maxArea; % Cross-sectional area of the supposed tumor

end

% Performs preprocessing and segmentation on a volume
function [] = segmentation3(volume, window)
    for k = 1:size(volume, 3)
        slice_k = slice(volume, k);

        % TODO: Can feed to seg2 alreayd cropped
        [tumor, area] = segmentation2(slice_k, window);

        % TODO: Decide function behavior
        imshow(tumor, [], 'InitialMagnification', 'fit')
        title(['Cross-Sectional Area: ' num2str(area)])
        drawnow, pause(0.05)
        % TODO: Create a montage instead of an animation (RGB problem)
        % TODO: Understand how to display the cross-sectional area values
        % axial_volume(:, :, k) = tumor;
    end
end

% Adds noise to a volume
function [vol_out] = addnoise(vol_in, noise, degree)
    % tmp = zeros(size(vol_in));
    for k = 1:1:size(vol_in, 3)
        im = slice(vol_in, k);
        im_n = imnoise(im, noise, degree);
        vol_out(:, :, k) = im_n;
    end
    % vol_out = tmp;
end

% Assesses if a power point transformation before the implementation of our
% workflow enhances or worsens the performance by calculating the Dice
% coefficient relatively to the true segmentation (ground truth)
function [] = assess_gamma(volume, type, ground_truth)
    [sel_vol, sel_wnd] = selectvol(volume, type); 
    D = zeros(size(sel_vol, 3), 20); % 20 defined by gamma range
    for i = 1:size(sel_vol, 3)
        for gamma = 0.1:0.1:2
            tmp = imadjust(sel_vol(:,:,i), [0 1], [0 1], gamma);

            tmp = segmentation2(tmp, sel_wnd);
            % TODO: substitute w/ our seg2
            % tmp = imbinarize(tmp, graythresh(tmp));

            D(i, int8(gamma * 10)) = dice(rot90(squeeze(ground_truth(i, :, :))), tmp);
        end
    end
    figure
    for i = 1:size(sel_vol, 3)
        plot(D(i, :)), ylim([0 1]), xline(10,'LineWidth',3)
        pause(0.05), drawnow
    end
end
