# PTOA-MRI-T1rho-and-T2star
Matlab files for processing PTOA knee MRI cartilage segmentations and calculating T1rho and T2* 

A collection of Matlab M-files for reading Philips DICOM MRI T1rho and T2star knee image data, registering the spin lock or echo times using Melastix and elastix, reading OsiriX segmentation CSV files and calculating T1rho or T2star for regions of interest (ROIs).

The ROIs are based on the lateral and medial compartments of the femur and tibia.  The compartments are divided into anterior (trochlea in the femur), central, and posterior regions.  These cartilage ROIs are further divided into deep and superficial layers.  Thus, there are 24 regions in each knee:

 1. femur lateral trochlea deep ROI,
 2. femur lateral trochlea superficial ROI,
 3. femur lateral central deep ROI,
 4. femur lateral central superficial ROI,
 5. femur lateral posterior deep ROI,
 6. femur lateral posterior superficial ROI,
 7. femur medial trochlea deep ROI,
 8. femur medial trochlea superficial ROI,
 9. femur medial central deep ROI,
10. femur medial central superficial ROI,
11. femur medial posterior deep ROI,
12. femur medial posterior superficial ROI,
13. tibia lateral anterior deep ROI,
14. tibia lateral anterior superficial ROI,
15. tibia lateral central deep ROI,
16. tibia lateral central superficial ROI,
17. tibia lateral posterior deep ROI,
18. tibia lateral posterior superficial ROI,
19. tibia medial trochlea deep ROI,
20. tibia medial trochlea superficial ROI,
21. tibia medial central deep ROI,
22. tibia medial central superficial ROI,
23. tibia medial posterior deep ROI,
24. tibia medial posterior superficial ROI,

Please see PTOA_ImageAnalysisPipeline2.pdf for a flowchart for the sequence of execution of Matlab M-files.

Please see the header comments in the M-files.  Dependent functions are listed in the comments.  There are a few helper M-files and a few unused M-files that are not required for the T1rho and T2* analyzes.

Note that the data must be in subject directories under visit subdirectories (BASELINE, YEAR1, or YEAR2).
