# RUSH3D — Real-time, Ultra-large-Scale, High-resolution 3D imaging platform

The RUSH3D achieves uniform resolutions of 2.6×2.6×6 µm3 across a volume of 8,000 × 6,000 × 400 µm3, at 20 volumes per second with low phototoxicity.

Please submit any issues or questions to this repository's [issue tracker](https://github.com/yuanlong-o/RUSH3D/issues).

## System requirements
* Software:
  * Matlab-compatible version of Windows or Linux (see https://www.mathworks.com/support/requirements/matlab-system-requirements.html)
  * Matlab R2021a or newer
  * Toolboxes: signal processing, image processing, statistics and machine learning, curve fitting, parallel computing
  * Matlab-compatible CUDA driver (>= 10.1 for Matlab R2021a)
  * Tested on Windows 10 with Matlab R2021a
* Hardware:
  * Matlab system requirements (see https://www.mathworks.com/support/requirements/matlab-system-requirements.html)
  * RAM should be about twice the size of the raw data in one patch. 
  * Matlab-compatible Nvidia GPU (see https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html], >~ 10 GB GPU-RAM recommended
  * Tested on a workstation with two Intel Xeon Gold 6136 CPUs (12 cores each), 512 GB RAM, 4 TB NVMe flash disk, three Nvidia Titan V GPUs with 12 GB GPU-RAM each
  * Mechanical design files of scanning light-field sensor (see `camera mechanical design-upload`)


## Usage
* You can download the raw data and system PSF from: https://drive.google.com/drive/folders/1Vu0NbfKpJaAzkTZOIb1kqtKEy2WZWVjL (120 frames and PSF from 12 ROIs).
* First, run the file `main_global.m` This file is used to convert the light-field rawdata to different view images (pixel realign) and divided the whole FOV into several blocks.
  * Make sure the `centerX` and `centerY` (coordinate of center pixel of center microlens of whole FOV) and other parameters correct.
  * Put the raw data and PSF in correct dir. 
* Second, run the `main_patch.bat`. This batch file will open several matlab scripts one by one and run automatically, processing corresponding patches to extract neuron signal. Three questions need to be answered:
  * dir of cod file. eg: `D:\RUSH3D\va.31`.
  * patch(es) need(s) to be processed. You have to follow the form, eg: `(1,2,3)`. Your input must start with `(` and end with `)`, and use `,` as delimiters.
  * delay time (s). A `.m` process will start to run after an excat delay time to avoid confict of memory or other computational resources. eg: `200`. Over 200 is recommended.
* Finally, run the `collect_all_trace.m` after set the result file directory above as `filepath`. This batch is used to concatenate all neurons information from different patches, providing a global result.


