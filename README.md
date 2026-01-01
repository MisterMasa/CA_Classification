Carbonic Anhydrase Motif Classification
This document describes a Python-based tool for classifying protein sequences of carbonic anhydrases (CA) into Alpha, Beta, Gamma, or None families based on motif detection using regular expressions. It is designed for reproducible bioinformatics workflows and outputs results in standard FASTA, JSON, and text formats.
Features
•	Reads input sequences from a FASTA file (CA.fasta).
•	Detects motifs corresponding to Alpha, Beta, and Gamma carbonic anhydrases.
•	Classifies each sequence into one of the families: Alpha, Beta, Gamma, None, or Ambiguous.
•	Outputs:
•	- alpha.txt, beta.txt, gamma.txt, none.txt: classified sequences in FASTA format.
•	- detailed_results.json: full classification details for each sequence.
•	- statistics.txt: overall statistics of classification results.
Usage
Clone the repository:
•	git clone https://github.com/USERNAME/CA_Classification.git
•	cd CA_Classification
Place your input FASTA file as CA.fasta in the project directory.
Run the script:
•	python CA_Class.py
Check the output files:
•	- Classified FASTA files (alpha.txt, beta.txt, gamma.txt, none.txt)
•	- statistics.txt for summary statistics
•	- detailed_results.json for detailed results
Example Output
Statistics:
Beta: 183 (91.5%)
Alpha: 10 (5.0%)
None: 7 (3.5%)
Citation
If you use this code in your research or project, please cite the following article:

Gheibzadeh MS, Manyumwa CV, Tastan Bishop Ö, Shahbani Zahiri H, Parkkila S, Zolfaghari Emameh R.
Genome Study of α-, β-, and γ-Carbonic Anhydrases from the Thermophilic Microbiome of Marine Hydrothermal Vent Ecosystems.
Biology (Basel). 2023 May 25;12(6):770. doi: 10.3390/biology12060770. PMID: 37372055; PMCID: PMC10295459.
License
You are free to use and modify this code. Please ensure proper citation of the above article when using this repository in academic or professional work.
