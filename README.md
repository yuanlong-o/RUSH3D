# RUSH3D — Real-time, Ultra-large-Scale, High-resolution 3D imaging platform

The RUSH3D achieves uniform resolutions of 2.6×2.6×6 µm3 across a volume of 8,000 × 6,000 × 400 µm3, at 20 volumes per second with low phototoxicity.

Please submit any issues or questions to this repository's [issue tracker](https://github.com/yuanlong-o/RUSH3D/issues).

## System requirements
* Software:
  * Matlab-compatible version of Windows or Linux (see https://www.mathworks.com/support/requirements/matlab-system-requirements.html)
  * Matlab R2020a or newer
  * Toolboxes: signal processing, image processing, statistics and machine learning, curve fitting, parallel computing
  * Matlab-compatible CUDA driver (>= 10.1 for Matlab R2020a)
  * Tested on Ubuntu Linux 20.04 with Matlab R2020a, CUDA driver v11, CUDA toolkit v10.1
* Hardware:
  * Matlab system requirements (see https://www.mathworks.com/support/requirements/matlab-system-requirements.html)
  * RAM should be about twice the size of the raw data in one patch. 
  * Matlab-compatible Nvidia GPU (see https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html], >~ 10 GB GPU-RAM recommended
  * Tested on a workstation with two Intel Xeon Gold 6136 CPUs (12 cores each), 256 GB RAM, 4 TB NVMe flash disk, three Nvidia Titan V GPUs with 12 GB GPU-RAM each

## Installation
1. Add the src folder and all subfolders to Matlab path
2. Ensure Matlab's `mex` compiler is configured with a compatible compiler, as listed here: https://www.mathworks.com/support/requirements/supported-compilers.html
3. Change directory to subfolder `src/cuda` and run `build.m` to run the mexcuda compiler to build CUDA and C++/MEX binaries. This should result in a file `mex_lfm_convolution_vec.mexa64` (file extension may differ depending on platform)

Typical installation time: 5 minutes. Download times depend on connection speed; should be <1 hour on a reasonably fast connection

## Usage
* See comments in parameter setting script for documentation of required and optional arguments
* In general, to run the RUSH3D pipeline, pass at least the required arguments to the main function. 
  
  * Parameters `frames_x_offset`, `frames_y_offset`, `frames_dx` (central microlens offets in x and y and microlens pitch, all in units of pixels) can be conveniently determined using the LFDisplay software published by the Levoy lab at Stanford University: http://graphics.stanford.edu/software/LFDisplay/
  * Parameter `indir` should point to a folder containing a series of .tif files containing one frame of raw data each. Files will be read in in alphabetic order.
  * Replace `<psfdir>` with the path to the directory containing the PSF file.
