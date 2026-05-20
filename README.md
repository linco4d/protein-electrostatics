# protein-electrostatics

This project explores how structural protein families can be sampled, cleaned, and prepared for electrostatics analysis at scale. I used Python to pull classification and structure data from CATH and the RCSB Protein Data Bank, organize representative proteins by fold family, clean raw structures with tools such as PDBFixer and OpenMM, and prepare fragments for downstream electrostatic potential calculations with `pdb2pqr` and APBS. The work gave me hands-on experience building small research pipelines, handling imperfect scientific data, and combining scripting, file processing, visualization, and external bioinformatics tools into a reproducible workflow.

`src/cath.py` is an early parsing script for the raw CATH domain files. It was useful for validating how domain identifiers map onto the CATH hierarchy and for shaping the first version of the data model.

`src/cath_pipeline.py` is the main dataset construction script. It downloads current CATH classification tables, parses domain metadata into TSV outputs, groups entries by fold category, and selects representative structures for each superfamily using either deterministic spacing or random sampling.

`src/get_cath_domains.py` downloads domain-level structures directly from the CATH API. It uses parallel requests, retry logic, and resume-friendly checks to make large pulls more reliable.

`src/cath_rossmann_filter.py` focuses on a Rossmann-like subset of the dataset. It filters domains by CATH code, annotates higher-confidence superfamily matches, and can generate FASTA output from local mmCIF files with `gemmi`.

`src/fetch_rep_structures.py` takes representative PDB IDs from a TSV and downloads full PDB or mmCIF files from RCSB. This is the bridge between the classification tables and the structure-processing stage of the workflow.

`src/failed_mmcifs.py` retries only the structures that failed during an earlier mmCIF download pass. It makes the pipeline easier to recover without rerunning a full batch job.

`src/preproc.py` handles structure cleaning and fragment generation. It uses PDBFixer and OpenMM when available, falls back to a lighter ATOM-only cleanup path when needed, and then slices cleaned proteins into residue k-mers for localized downstream analysis.

`src/calc_electrostatics.py` prepares cleaned fragments for electrostatics calculations. It converts PDB fragments to PQR with `pdb2pqr`, patches APBS input settings, and runs APBS to generate electrostatic potential maps.

`src/frequency_analysis.py` summarizes how often different CATH architectures appear in the processed dataset. It produces tabular outputs that help quantify which fold classes are common, rare, or overrepresented.

`src/fold_type_distribution.py` turns those distributions into plots across multiple CATH levels, including class, architecture, topology, and homologous superfamily. It was useful for quick exploratory analysis and for spotting broad structural trends.

`src/identify_motifs.py` visualizes common versus rare structural motifs from the architecture frequency tables. It produces figures that make the dataset easier to interpret at a glance.
