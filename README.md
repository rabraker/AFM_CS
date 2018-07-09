# AFM_CS

This branch contains all of the code used in our paper
R. A. Braker, Y. Luo, L. Y. Pao and S. B. Andersson  
"Hardware Demonstration of Atomic Force Microscopy Imaging via
Compressive Sensing and $\mu$-path Scans", American Control Conf., July 2018. pp. 6037-6042.


The relevant directory structure looks like:

├── afm_imaging
│   ├── acc2018-data
│   │   ├── cs-data
│   │   ├── force-map
│   │   ├── frf_data
│   │   ├── raster
│   │   └── test
│   ├── controls
│   ├── data
│   │   ├── cs-data
│   │   └── sys_id
│   ├── fpga_VIs
│   │   ├── floating_point_toolkit
│   │   └── fpga_subVIs
│   ├── matlab-code
│   │   ├── classes
│   │   ├── dev
│   │   ├── examples
│   │   ├── figures
│   │   ├── force-modeling
│   │   ├── functions
│   │   ├── models
│   │   ├── system_id
│   │   └── tests
│   ├── sandbox
│   ├── subVIs
│   ├── sysID
│   └── UnitTests
├── latex
│   ├── data
│   └── figures
├── latex-rev1
│   ├── data
│   └── figures
├── libraries_host
├── presentation
│   └── figures
│       ├── blocks
├── reconstruction
│   ├── BP
│   │   └── l1magic
│   ├── SMP
│   └── testfile
├── SMP_1D


# Directory description.
## afm_imaging
All of the analysis and experimental code is in this directory. This includes both Labview code for both the host side and the FPGA. the sub-directory matlab-code contains code to generate CS-patterns, raster patterns as well as all the analysis code used to post-process the data. 

The actual reconstruction algorithm is from an external library, l1-magic, and is contained in the higher leval directory 
`recustruction/BP/l1magic.`

The actual raw data files are not included here due to their size (over half a gigabyte). 

## latex an latex-rev1 and presentation
These directories contain the manuscript latex code and figures for the initial submisssion (`latex`) and the revision (`latex-rev1`) while `presentation` contains all the latex and figures for the ACC conference presentation.

## 
 
