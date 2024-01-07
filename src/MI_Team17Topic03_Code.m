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

% Extracting the orthogonal planes according to the ranges and permuting 
% the volume such that the third dimension iterates through the slices of 
% interest (see selectvol())
% sg_planes = rot90(permute(squeeze(vol(100:145,:,:)),[2,3,1]));
% ax_planes = rot90(squeeze(vol(:,:,60:90)));
% cr_planes = rot90(permute(squeeze(vol(:,140:180,:)),[1,3,2]));

% For each plane identifying the slice with the biggest lesion extention and
% using it to empirically define the ROI window using imcrop()
% sagittal_window = [140 25 41 31]; % Original sagittal slice 125
% axial_window = [105 80 41 46]; % Original axial slice 79
% coronal_window = [100 20 46 31]; % Original coronal slice 158

%% 3. SEGMENTATION OF SAGITTAL SLICE 135

og_num = 135; % From reference
actual_num = 135 - sagittal_range(1); % Slice number
[sg_vol, sg_wnd] = selectvol(vol, 'sagittal'); % Selecting sagittal volume and ROI
slice_k = slice(sg_vol, actual_num); % Extracting the slice

diocsn = imcrop(slice_k, [], sg_wnd); % TODO OOOOH

[tumor_k1, area_k] = segment(diocsn); % Segmenting the lesion

tumor1 = overlay(slice_k, tumor_k1, sg_wnd); % TODO DIOCAN

% Displaying the segmentation
imshow(tumor1, [], 'InitialMagnification', 'fit')
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
    show_volume_segmentation(sel_volume, sel_window);
    clc, close all
end

%% NOISE SENSITIVITY ANALYSIS
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
    [~, area_filt] = segment(im_noise, sg_wnd);
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

% The same test performed with gaussian noise showed that, regardless of
% filtering, the segmentation performance was significantly worsened.
% Therefore we decided upon not applying a preemptive averaging filter.

%% 6. PERFORMANCE ANALSYSIS WITH DICE COEFFICIENT

% The file 'true_segmentation.mat' was created using the volume segmenter
% application, applying Otsu thresholding and smoothing edges + removing
% some excesses manually
true_seg = load('true_segmentation.mat').labels; % Ground truth

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

    switch volume_type
        case 'axial'
            ground_truth = true_seg;
        case 'sagittal'
            ground_truth = permute(true_seg, [2 3 1]);
        case 'coronal'
            ground_truth = permute(true_seg, [3 1 2]);
    end
    
    % Assessing Dice coefficient
    assess_gamma(vol, volume_type, ground_truth)
    clc, close all
end

%% TODO: Degub
ground_truth = true_seg;
assess_gamma(vol, 'axial',  ground_truth)

%% FUNCTIONS

% Extracts the orthogonal planes according to the ranges and permuting 
% the volume such that the third dimension iterates through the slices of 
% interest
function [vol_out, window_out] = selectvol(volume, type)
    switch type
        case 'axial'
            vol_out = rot90(squeeze(volume(:, :, 60:90)));
            window_out = [105 80 45 40];
        case 'sagittal'
            vol_out = rot90(permute(squeeze(volume(100:145, :, :)), [2, 3, 1]));
            window_out = [140 25 40 30];
        case 'coronal'
            vol_out = rot90(permute(squeeze(volume(:, 140:180, :)),[1, 3, 2]));
            window_out = [105 20 30 45];
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
% Performs preprocessing and segmentation on a slice by assuming that the
% lesion is the largest relatively dense area in the ROI
function [im_out, area] = segment(im_in)
    tmp_im = medfilt2(im_in); % Applying preemptive median filtering
    tmp_im = im2double(tmp_im); % Standardizing the gray scale
    tmp_im = tmp_im > 0.50 & tmp_im < 0.85; % Empirical thresholding
    label = bwlabel(tmp_im); % Labeling continuous regions
    stats = regionprops(logical(tmp_im), 'Solidity', 'Area'); % Stats of continuous regions
    density = [stats.Solidity]; area = [stats.Area]; disp(size(area))
    denseArea = density > 0.6; % Areas with empirical density
    maxArea = max(area(denseArea)); disp(size(maxArea)), disp('----') % Largest area with empirical density
    lesionLabel = find(area == maxArea); % Label of the supposed tumor
    lesion = ismember(label, lesionLabel); % Extracting the image from the label
    im_out = imfill(lesion, 'holes'); % Correcting the image
    area = maxArea; % Cross-sectional area of the supposed tumor
end

% Performs preprocessing and segmentation on a volume
function [] = show_volume_segmentation(volume, window)
    for k = 1:size(volume, 3)
        slice_k = slice(volume, k); % Loading the slice
        cropped_slice_k = imcrop(slice_k, [], window); % Cropping out the ROI
        [lesion, area] = segment(cropped_slice_k); % Isolating the lesion
        segmented_lesion = overlay(slice_k, lesion, window); % Segmenting the lesion
        imshow(segmented_lesion, [], 'InitialMagnification', 'fit') % Displaying the segmentation
        title(['Cross-Sectional Area: ' num2str(area)])
        drawnow, pause(0.05)
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

% Assesses whether a gamma correction enhances or worsens the performance
% of our segmentation method by calculating the Dice coefficient and
% comparing it with Otsu's method
function [] = assess_gamma(volume, type, ground_truth)
    figure, [sel_vol, sel_wnd] = selectvol(volume, type); % Selecting

    D1 = zeros(size(sel_vol, 3), 20); % 20 defined by gamma range
    D2 = zeros(size(sel_vol, 3), 20); % 20 defined by gamma range
    for i = 1:size(sel_vol, 3)
        for gamma = 0.1:0.1:2
            tmp = imadjust(sel_vol(:, :, i), [0 1], [0 1], gamma);
            tmp = imcrop(tmp, [], sel_wnd); % Cropping out the ROI
            tmp1 = segment(tmp); % Isolating the lesion with our workflow
            tmp2 = imbinarize(tmp, graythresh(tmp)); % Isolating the lesion with Otsu's method

            % Saving the Dice coefficients
            D1(i, int8(gamma * 10)) = dice(rot90(squeeze(ground_truth(:, :, i))), tmp1);
            D2(i, int8(gamma * 10)) = dice(rot90(squeeze(ground_truth(:, :, i))), tmp2);
        end

        plot(D1(i, :)), hold on, plot(D2(i, :))
        xline(10,'LineWidth',3) % i-th slice with gamma = 1
        ylim([0 1]), xlim([0 20])
        ylabel('Dice Coefficient'), xlabel('Brighter <-- 10 γ --> Darker')
        legend({'Our Workflow', 'Otsu', ' γ = 1'}, 'Location', 'southwest')
        title(['Analysis over ' type ' planes'])
        pause(0.05), drawnow, hold off
    end
end
