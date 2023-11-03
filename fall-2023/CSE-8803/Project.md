# Project notes

## Project proposal

https://www.overleaf.com/project/651d9b2fea99146e108c545d

## Project proposal feedback

- I suggest you make sure about the data availability sooner rather than later, e.g., if you can find the data you need, including the patient-level data (the exposed pathogens for a given patient). It is unclear from the proposal if you have obtained such data.
- Also make sure you are able to find the structure data you need. I assume you need to map the sequence in IDEB to PDB? This process can take some time, so better to start early and make a plan B.
- Are there any existing studies solving the same problem that you formulated in the “Problem Description” section? If so, check if those publications already provided processed data that you can use.
- The proposal mentioned that the entire protein structure will be represented as a graph -- please start testing your GNN early, especially to see if you can fit your GNN on a GPU. If not, may need to consider only using the binding surface of the protein rather than the entire structure in the GNN.

## Data

Antigenic distances:
https://figshare.com/collections/A_benchmark_dataset_of_protein_antigens_for_antigenicity_measurement/4961501/2

## 3D protein representations

- Links from "Geometric Latent Diffusion Models for 3D Molecule Generation" paper:
  - Represent molecules as atomic graphs in 3D:
    - "Schnet: A continuous-filter convolutional neural network for modeling quantum interactions.""
  - **Represent proteins as proximity spatial graphs over amino acids:**
    - "Learning from protein structure with geometric vector perceptrons."
  - See "Molecule Generation in 3D" section.
- ScanNet: an interpretable geometric deep learning model for structure-based protein binding site prediction
  - Here, we introduce ScanNet, an end-to-end, interpretable geometric deep learning model that learns features directly from 3D structures. ScanNet builds representations of atoms and amino acids based on the spatio-chemical arrangement of their neighbors. We train ScanNet for detecting protein–protein and protein–antibody binding site

## Antibody fit

https://vina.scripps.edu/ - "Molecular docking" software.

## Architecture

Also from "Geometric Latent Diffusion Models for 3D Molecule Generation" paper:
- DMs (diffusion models) have shown promising results, i.e., target drug design (Lin et al., 2022) and antigen-specific antibody generation.


## Backup

Instead of focusing on patient-specific prediction, since this data may be hard to find, we could just predict antigenic distance.

"antigen-specific antibody generation" (note from "Geometric Latent Diffusion Models for 3D Molecule Generation" above) may be a good backup idea.
