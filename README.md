# NIRFP_design


This repository contains scripts for designing covalent small molecule binders, with a focus on targeting fluorogenic small molecules to create near-infrared fluorescent proteins (NIRFPs). Specifically, we used endogenous biliverdin to develop NIR-I fluorescent proteins (650–1000 nm) and the cyanine dye FD-1080 to develop NIR-II fluorescent proteins (1000–1700 nm). For further background information, please refer to our manuscript.


Although our primary application is fluorescent proteins, all scripts in this repository are generalizable and can be used to design proteins that covalently bind other types of small molecule ligands.


The pipeline features a physics-based generative approach for constructing covalent small molecule binders using only the configuration and chemistry of the target small molecule. The design process is divided into four main steps:


---


## Pipeline Overview


### 1. Generating the Rotamer Interaction Field (RIF)


The rotamer interaction field (RIF) for each small molecule ligand is generated using the `rifgen` app from the **rifdock** package. The resulting RIF is saved in phmap (parallel hashmap) format, allowing for efficient look-up during structure generation. Example scripts can be found in the `1_rif` subfolder.


### 2. Building the Protein Backbone Structure


We use the `covalent_ligand_binder` app in our **ProBuilder** package to construct covalent small molecule binders in a bottom-up fashion. An extended protein chain is sampled with randomly assigned secondary structure elements. The covalently bound small molecule ligand is modeled as a noncanonical amino acid (NCAA), and a motif containing this NCAA (or a random rotamer of it) is randomly placed onto the extended chain. Monte Carlo fragment assembly is then used to fold the chain into a well-structured protein, guided by the RIF and RPX score. Final designs feature the ligand tightly packed and covalently tethered to the protein.


### 3. Sequence Optimization


The generated structures undergo sequence design using **Rosetta** and **ProteinMPNN** to optimize both protein folding and ligand interaction. The covalently bound ligand is modeled as a NCAA, which is kept fixed throughout the sequence optimization process.


### 4. Design Filtering and Selection


Final designs are filtered based on ligand-binding energy, shape complementarity, hydrogen-bonding networks, and structural agreement with **AlphaFold2** predictions. The best designs are then selected for experimental characterization.


---


## Getting Started


All scripts and example command lines are organized in the respective subfolders. This repository will be continuously updated to improve user-friendliness and further facilitate the design of covalent small molecule binders for next-generation tool and therapeutic development.


---


If you use this repository or the associated scripts, please cite our manuscript.