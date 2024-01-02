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
sagittal_range = [105 145]; axial_range = [60 90]; % coronal_range = []; % Empirical
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

% Segmenting the tumor
k = 135; % Slice number
slice_k = slice(selectvol(vol, 'sagittal'), k); % Extracting the slice
[tumor_k, area_k] = segmentation2(slice_k, sagittal_window);

% Displaying the segmentation
imshow(tumor_k, [], 'InitialMagnification', 'fit')
title(['Sagittal Slice: ' int2str(k) ' - Cross-Sectional Area: ' int2str(area_k)])
drawnow

%% SEGMENTATION OF THE WHOLE VOLUME (WITH NOISE)

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
        while str2double(noise_level) < 0 | str2double(noise_level) > 1
            noise_level = input('Insert noise level: ', 's');
        end 

        % Supposing adding the noise is independent from slice type
        chosen_vol = addnoise(vol, noise_type, str2double(noise_level));

    end

    % TODO: Prompted slice choice but problem in range selection domain

    % Segmenting the volume
    [sel_volume, sel_range, sel_window] = selectvol(chosen_vol, volume_type);
    segmentation3(sel_volume, sel_range, sel_window);
    clc, close all
end

%% NOISE SENSITIVITY ANALYSIS

% Selecting a random slice as reference (supposing valid for all)
slice_k = slice(selectvol(vol, 'sagittal'), 135);
[seg_norm, area_norm] = segmentation2(slice_k, sagittal_window);

% Calculating the detected area relatively to different levels of noise
filtered = zeros(1, 100); noised = zeros(1, 100); i = 1;
for k = 0.01:0.01:1
    % Noise
    im_noise = imnoise(slice_k, 'salt & pepper', k);
    [seg_noise, area_noise] = segmentation2(im_noise, sagittal_window);
    if isempty(area_noise)
        area_noise = 0;
    end
    noised(i) = area_noise;

    % Filter
    filt = medfilt2(im_noise);
    [seg_filt, area_filt] = segmentation2(filt, sagittal_window);
    if isempty(area_filt)
        area_filt = 0;
    end
    filtered(i) = area_filt;

    i = i + 1;
end

figure, plot(filtered), hold on, plot(noised), hold on, refline(0, area_norm);
legend({'Filtered','Noise', 'Clean'}, 'Location', 'southwest')

% Deciding to always preventively apply the median filter to reduce
% eventual salt & pepper noise degradation

%% GAUSSIAN NOISE SENSITIVITY

h = 1/9.*ones(3,3);
tmp_im = imnoise(slice_k, "gaussian", 0.2);
subplot(2, 2, 1), imshow(tmp_im)
title('Gaussian Noise')

tmp_im = imfilter(tmp_im, h, 'conv');
subplot(2, 2, 2), imshow(tmp_im)
title('Filtered')

thr = 0.5;
db_im = im2double(tmp_im);
bw_im = db_im > thr; %& bw_im < 0.86;
subplot(2, 2, 3), imshow(bw_im)
title(['Threshold: ' num2str(thr)])

gamma = 3;
im_gamma = imadjust(db_im, [0 1], [0 1], gamma);
im_gamma = im_gamma > thr;
subplot(2, 2, 4), imshow(im_gamma, [], 'InitialMagnification', 'fit')
title(['Gamma' num2str(gamma)]), drawnow

%% experiments and notes

% xy 'axial' (iterate over third dimension of vol) plane intersting slices: from 60 to 90
% yz 'coronal' (iterate over second dimension of vol) plane interesting slices: from 140 to 180
% xz 'sagittal' (iterate over first dimension of vol) plane interesting slices: from 100 to 145
sg=rot90(permute(squeeze(vol(100:145,:,:)),[2,3,1])); %46 slices
ax=rot90(squeeze(vol(:,:,60:90))); %31 slices
cr=rot90(permute(squeeze(vol(:,140:180,:)),[1,3,2])); %41 slices

t_vol=vol(100:145,140:180,60:90);

% file 'true_segmentation.mat' was created using volume segmenter
% application, applying otsu thresholding and smoothing edges + removing
% some excess manually

true_seg=load('true_segmentation.mat').labels;

% estimation for gaussian noise may be done with std2 function, maybe
% comparing filtered and unfiltered image

%% sagittal test
for j=0.1:0.1:0.8
    sg=addnoise(sg,'gaussian',j);
    D=zeros(size(sg,3),20);
    for i=1:size(sg,3)
        for gamma=0.1:0.1:2
        tmp=imadjust(squeeze(sg(60:90,140:180,i)),[0 1], [0 1], gamma);
        tmp=imbinarize(tmp,graythresh(tmp));
        D(i,int8(gamma*10))=dice(rot90(squeeze(true_seg(i,:,:))),tmp);
        end
    end
    figure
    for i=1:size(sg,3)
        plot(D(i,:))
        hold on
    end
    xline(10,'LineWidth',3)
    %this means gamma<1 improves segmentation performance with otsu method
end
%testing with different noise levels: higher noise improve performance
%regardless of gamma, only for gaussian
%salt & pepper makes it worse but gamma doesn't have an effect for higher
%noise levels
sg=rot90(permute(squeeze(vol(100:145,:,:)),[2,3,1]));
%% axial test
D=zeros(size(ax,3),20);
for i=1:size(ax,3)
    for gamma=0.1:0.1:2
    tmp=imadjust(ax(100:145,140:180,i),[0 1], [0 1], gamma);
    tmp=imbinarize(tmp,graythresh(tmp));
    D(i,int8(gamma*10))=dice(squeeze(true_seg(:,:,i)),tmp);
    end
end
figure
for i=1:size(ax,3)
    plot(D(i,:))
    hold on
end
xline(10,'LineWidth',3)

%% coronal test
D=zeros(size(cr,3),20);
for i=1:size(cr,3)
    for gamma=0.1:0.1:2
    tmp=imadjust(cr(60:90,100:145,i),[0 1], [0 1], gamma);
    tmp=imbinarize(tmp,graythresh(tmp));
    D(i,int8(gamma*10))=dice(rot90(squeeze(true_seg(:,i,:))),tmp);
    end
end
figure
for i=1:size(cr,3)
    plot(D(i,:))
    hold on
end
xline(10,'LineWidth',3)

%% compare otsuthresh with our method
for j=0.1:0.1:0.8
    sg=addnoise(sg,'gaussian',j);
    D=zeros(size(sg,3),20);
    for i=1:size(sg,3)
        for gamma=0.1:0.1:2
        tmp=imadjust(squeeze(sg(60:90,140:180,i)),[0 1], [0 1], gamma);
        % segment with our method
        D(i,int8(gamma*10))=dice(rot90(squeeze(true_seg(i,:,:))),tmp);
        end
    end
    figure
    for i=1:size(sg,3)
        plot(D(i,:))
        hold on
    end
    xline(10,'LineWidth',3)
    %this means gamma<1 improves segmentation performance with otsu method
end
%% FUNCTIONS

% Takes as input the axial volume and returns a volume whose third dimension 
% are the slices in the plane of interest
function [vol_out, range_out, window_out] = selectvol(volume, type)
    switch type
        case 'axial'
            vol_out = volume;
            range_out = [60 90]; window_out = [138 107 42 39];
        case 'sagittal'
            vol_out = flip(flip(imrotate3(volume, 180, [0 -1 1]), 2), 3);
            range_out = [105 150]; window_out = [140 23 39 29];
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
    tmp_im = imcrop(im_in, [], rect); % Cropping out the ROI
    tmp_im = medfilt2(tmp_im); % Applying median filtering
    
    % h = 1/9.*ones(3,3); tmp_im = imfilter(tmp_im, h, 'conv');

    % tmp_im = imsharpen(tmp_im); % Applying image sharpening
    tmp_im = im2double(tmp_im); % Standardizing the gray scale

    % gamma = 2; tmp_im = imadjust(tmp_im, [0 1], [0 1], gamma);

    % TODO: Different meaning between images
    tmp_im = tmp_im > 0.50 & tmp_im < 0.85; % Thresholding % Empirical
    % tmp_im = edge(tmp_im, 'canny');
   
    label = bwlabel(tmp_im); % Labeling continuous regions

    % Stats of continuous regions
    stats = regionprops(logical(tmp_im), 'Solidity', 'Area');
    density = [stats.Solidity]; area = [stats.Area];

    % % Supposing the tumor is the largest dense area
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
        pause(0.05)
        % TODO: Create a montage instead of an animation (RGB problem)
        % TODO: Understand how to display the cross-sectional area values
        % axial_volume(:, :, k) = tumor;
    end
end

% TODO: Understnd why preallocating fucks things up
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
