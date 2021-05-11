# MATLAB implementation of B-COSFIRE filters
B-COSFIRE filters are non-linear trainable filters for detection of elongated patterns in images.  
This is the code of the trainable non-linear B-COSFIRE filters for delineation of elongated patterns
in images.  
The B-COSFIRE filters are proposed in the paper:

[Azzopardi, G., Strisciuglio, N., Vento, M., Petkov, N.: Trainable COSFIRE  filters for vessel delineation with application to retinal images. Medical Image Analysis 19(1), 46 - 57, 2015](http://www.cs.rug.nl/~george/articles/MEDIA2015.pdf)


## Applications
For applications of B-COSFIRE filters to different kinds of images and problems, please refer to the following codes.

### ExampleBloodVesselSegmentation.m - 
The example code for the configuration of a line detector and a line-ending detector and their 
applications to the segmentation of blood vessels in retinal images. 
The final response is the summation of the responses of the two filters. 


### PavementCrackDelineation.m (CAIP 2017)
Application of the B-COSFIRE filters in an image processing pipeline for the detection of cracks in pavement images presented at CAIP17.
The code provides the benchmark results reported in the paper
_Strisciuglio, N. Azzopardi, G. Petkov, N._ "Detection of Curved Lines with B-COSFIRE Filters: A Case Study on Crack Delineation", IWOBI 2017.

If you use this code, please cite the paper:

    @Inbook{Strisciuglio2017,
    author="Strisciuglio, Nicola and Azzopardi, George and Petkov, Nicolai",
    title="Detection of Curved Lines with B-COSFIRE Filters: A Case Study on Crack Delineation",
    bookTitle="Computer Analysis of Images and Patterns: 17th International Conference, CAIP 2017, Ystad, Sweden, August 22-24, 2017, Proceedings, Part I",
    year="2017",
    pages="108--120",
    isbn="978-3-319-64689-3",
    doi="10.1007/978-3-319-64689-3_9",
    }

If you use the data, please cite:
_Zou, Q., Li, Q., Zhang, F., Wang, Z.X.Q., Wang, Q.: Path voting based pavement crack detection from laser range images. In: IEEE ICDSP, pp. 432?436 (2016)_

The images used in this example are available from the 
[website of Dr. Zou](https://sites.google.com/site/qinzoucn/).

### INRIAImages.m (IWOBI 2017)
Application of the B-COSFIRE filters for detection of elongated structures in images.  
This code provides the benchmark results on the images of the INRIA data
set used in the paper  
_Strisciuglio, N. Petkov, N._ "Delineation of line patterns in images using B-COSFIRE filters", IWOBI 2017.

If you use this code, please cite the paper:

	@INPROCEEDINGS{7985538,
	author={N. Strisciuglio and N. Petkov},
	booktitle={2017 International Conference and Workshop on Bioinspired Intelligence (IWOBI)},
	title={Delineation of line patterns in images using B-COSFIRE filters},
	year={2017},
	pages={1-6},
	doi={10.1109/IWOBI.2017.7985538},
	month={July},}


The images used in this example are available at [this website](http://www-sop.inria.fr/members/Florent.Lafarge/benchmark/line-network_extraction/line-networks.html).


## Reference publications
If you use the code of B-COSFIRE filters for your application, please cite the following articles. 

__Original paper:__  

	@article{BCOSFIRE-MedIA2015,
	title = "Trainable {COSFIRE} filters for vessel delineation with application to retinal images ",
	journal = "Medical Image Analysis ",
	volume = "19",
	number = "1",
	pages = "46 - 57",
	year = "2015",
	note = "",
	issn = "1361-8415",
	doi = "http://dx.doi.org/10.1016/j.media.2014.08.002",
	author = "George Azzopardi and Nicola Strisciuglio and Mario Vento and Nicolai Petkov",
	} 
 
__Supervised learning of B-COSFIRE filters:__  

	@article{BCOSFIRE-selection2016,
	author={Strisciuglio, Nicola and Azzopardi, George and Vento, Mario and Petkov, Nicolai},
	title={Supervised vessel delineation in retinal fundus images with the automatic selection of {B-COSFIRE} filters},
	journal={Machine Vision and Applications},
	year={2016},
	pages={1?13},
	issn={1432-1769},
	} 


## Changelog
__20 Aug 2017__  
CrackDetectionCluster.m - Experimental code (and data) to replicate results in the CAIP17 paper.

__3 Jul 2017__  
applyCOSFIRE_inhib.m:132/137 - Approximatation of the shifting amount corrected  
The results published in the paper "_Strisciuglio, N. Petkov, N._ "Delineation of line patterns in images using B-COSFIRE filters", IWOBI 2017." are slightly different due to this bug fixing.

