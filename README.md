# RefineAssist
This script will calculate the solvent-accessible surface (SASA) for each amino acid within a protein and point out which amino acids are in the wrong layer.



## DESCRIPTION:
This script does the following:

1. Prints out PYMOL commands to select spesific amino acids in a protein depending on their layer calculated by their SASA (solvent-accessible surface area). These commands allow the user to colour different parts of a protein and better visualise it in order to take better decisions on refining it.
2. Prints out a sequence alignment which shows which amino acid belongs to which layers in the protein accodring to its SASA value, in which secondary structure does each amino acid belong to, the protein's amino acid sequence, and points out with a * which amino acids belong to the wrong layer and should be mutated.

This script was written to run on GNU/Linux, it was not tested in Windows or MacOS.
This script will mostly be useful to refine proteins after a Rosetta FFL (Fold From Loop) computation, but can still be used to refine any protein.
Contact the author at sari.sabban@gmail.com for any questions regarding this script.



## HOW TO USE:
To use follow these steps:

1. Install biopython by running the following command in terminal (pip install biopython) you need pip to be installed, if you do not have it you can install it in linux by running the following command in terminal (sudo apt-get install python-pip).
2. Install DSSP in linux by running the following command in terminal (sudo apt-get install dssp).
3. Enter file name when prompted. File must be in the same directory as this script.



## REFERENCES:

1. Refere to the paper by (Koga et.al., 2012 - [PMID: 23135467](https://www.ncbi.nlm.nih.gov/pubmed/23135467)) Methods section's Sequence Design Protocol for more explanation on protein refinement and each layer's SASA parameters.
2. Refere to the paper by (Correia et.al., 2014 - [PMID: 24499818](https://www.ncbi.nlm.nih.gov/pubmed/24499818)) for details about the Rosetta Fold From Loop (FFL) protocol.
3. Refere to this [webpage](goo.gl/NsQubf) for details about Rosetta LayerDesign Protocol and which amino acids should be in which layer.
4. Refere to this [webpage](http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ) for how to use BioPython.
5. Refere to the references (Cock et.al., 2009 - [PMID: 19304878](https://www.ncbi.nlm.nih.gov/pubmed/19304878)) and (Hamelryck and Manderick, 2003 - [PMID: 14630660](https://www.ncbi.nlm.nih.gov/pubmed/14630660)) regarding the applications used here from biopython.
6. Refere to the references (Kabsch W. and Sander C., 1983 - [PMID: 6667333](https://www.ncbi.nlm.nih.gov/pubmed/6667333)) regading DSSP.
7. Refere to the reference (Tien et.al., 2013 - [PMID: 24278298](https://www.ncbi.nlm.nih.gov/pubmed/24278298)) regarding Wilke SASA parameters.
