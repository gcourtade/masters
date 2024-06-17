# A Computational Pipeline to Study Donor Substrate Binding to the Wild-Type and Mutants of a Xanthan Gum Glycosyltransferase

## Tova Alenfalk
This repository contains the scripts and results from my master thesis "A Computational Pipeline to Study Donor Substrate Binding to the Wild-Type and Mutants of a Xanthan Gum Glycosyltransferase". 

## Content:
The coordinate files for the different wild type conformations and mutants used in the thesis, as well as the different ligands, is found in *coordinate_files*. This folder also contains the autobox files used to the define the search box for the docking. 

*Docking* contains the script for generating an ensemble of docked poses using the docking program GNINA. It outputs an sdf file containing the docked poses and an excel sheet with the docking scores and the calculated attributes for the poses (i.e. the uridine RMSD, distance between the C6 atom of the donor substrate sugar and residue 307, and dihedral angles).

The result files (sdf files with poses and corresponding excel files) from the docking done in my thesis is found in *results*. The representative UDP poses for the different WT conformations selected from the docked UDP poses (used to calculate the uridine RMSD) is also found in this folder. The results are further analysed in the script *analysis_docking_results*.
