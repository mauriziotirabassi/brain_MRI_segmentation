# Brain Lesion Segmentation in MRI Images

## Overview

This project implements a **robust pipeline for automatic brain lesion segmentation** from MRI images. The workflow extracts the lesion region, highlights its surface, and calculates cross-sectional area. The method is tested for sensitivity to image noise, demonstrating robustness against “salt & pepper” noise.

## Features

- **Automatic Lesion Detection:** Segment lesions using filtering, thresholding, and region property analysis.  
- **ROI-Based Processing:** Lesion identified in empirically defined regions of interest.  
- **Noise Robustness:** Tested against additive Gaussian and salt & pepper noise.  
- **Accuracy Assessment:** Dice coefficient used to compare automated segmentation with manual ground truth.  
- **Visualization:** Sagittal, coronal, and axial views with optional noise addition via MATLAB CLI.  
- **Comparison:** Outperforms or matches Otsu thresholding in baseline images; handles noise effectively.  

## Methodology

1. **Volume Visualization:** Used `orthosliceViewer()` and `imcrop()` to select relevant slices and define ROI.  
2. **Lesion Segmentation:**  
   - Threshold standardized image (0.5–0.85 range).  
   - Label continuous regions with `bwlabel()`.  
   - Compute region properties (`solidity`, `area`) using `bwregionprops()`.  
   - Fill contiguous lesion areas using `imfill()`.  
3. **Noise Sensitivity Analysis:** Added varying levels of Gaussian and salt & pepper noise; applied median filter for robustness.  
4. **Manual Segmentation:** Used MATLAB Volume Segmenter as ground truth for Dice coefficient evaluation.  
5. **Dice Coefficient:** Evaluates segmentation accuracy (0 = no overlap, 1 = perfect overlap).  

## Results

- **Sagittal Slice 135:** Lesion cross-sectional area detected: 332 pixels.  
- **Volume Segmentation:** CLI allows visualization across planes, with optional noise simulation.  
- **Noise Performance:** Preemptive median filtering improves accuracy against salt & pepper noise.  
- **Comparison:** The proposed method is resilient to gamma corrections and matches or outperforms Otsu thresholding.  
