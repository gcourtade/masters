# Computational analysis of _de novo_ enzymes

This GitHub repository contains all code used in the computational analysis in Ingrid Moa Bolstad's specialisation project on generating _de novo_ PET-hydrolyzing enzymes.



## Contents
The folder _Data_ contains all input files used in the computational analysis, and output text files are saved there as well. All output figures are stored in _Figures_.

_rmsd\_from\_files.py_ finds the catalytic and binding site root-mean-square deviation (RMSD) between the novel proteins and LCC-ICCG as a reference using [MDAnalysis](https://docs.mdanalysis.org/1.1.0/documentation_pages/analysis/rms.html). The results are sorted from low to high and plotted as barplots.

_molprobity\_pipeline.py_ is the script used to conduct the batch MolProbity analysis of all _de novo_ proteins and outputs a MolProbity summary table for each protein (_Analysis\_results.txt_). It requires local installation of the [MolProbity software](https://github.com/rlabduke/MolProbity).

_protein\_selection.py_ sorts proteins by their MolProbity score and by their Ramachandran outliers and plots the results as barplots.

_AF\_parameters.py_ sorts the proteins by their AlphaFold prediction parameters pLDDT and prediction-to-backbone RMSD and the results are plotted as barplots.

_ramachandran\_plots.py_ plots Ramachandran distribution plots using [MDAnalysis](https://docs.mdanalysis.org/stable/documentation_pages/analysis/dihedrals.html#MDAnalysis.analysis.dihedrals.Ramachandran).

_avg\_residue\_numbers.py_ finds the average residue numbers for each input design.


