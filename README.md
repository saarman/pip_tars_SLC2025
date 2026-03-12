# Ecological partitioning of Culex pipiens and tarsalis in SLC, Utah

**Brief overview:**  
  This repository supports analyses of *Culex pipiens* and *Culex tarsalis* ecological partitioning in Salt Lake City (SLC), Utah.  

Core objectives include mapping spatial and temporal partitioning, and modeling occurance and abundance with GLMM. The repository contains all code and small metadata files for these analyses Large and/or geospatial inputs, as well as final outputs, are stored externally (see below).

## Folder structure

- **`scripts/`**  (see `scripts/README.md` for table of contents)   
– R, Python, and Bash scripts in numbered order (e.g., `01_prepare_data.R`, `02_map_weekly.R`).  

- **`help/`**  (Long-form help and documentation; see `help/README.md` for details)  
1. Connecting to RStudio on CHPC  
2. Submitting Slurm (interactive & batch) jobs  
3. Linking GitHub with RStudio (creating a PAT, cloning, committing, pushing)  

- **`input/`** (e.g., sample metadata, lookup tables)      
– Small reference files (e.g., sample metadata, lookup tables).      
– *NOTE*: Does not include large raw genomic and geospatial inputs, these are stored externally (see below).

- **`data/`** (ignored by GitHub; see `data/README.md` for external locations)

- **`results/`** (ignored by GitHub; see `results/README.md` for external locations)

- **`figures/`**  (see `figures/README.md` for list of figures; larger files ignored by GitHub)      
– Small output figures.      
– Large high-resolution figures ignored, see `figures/fignored`


## Getting started

1. Clone the repo:
  ```bash
git clone https://github.com/saarman/pip_tars_SLC2025.git
```

2. Install required R packages:
  ```bash
Rscript scripts/00_install_packages.R
```
3. Follow the pipeline:
  - Run `01_prepare_data.R` to load and format data.
- Then run subsequent scripts in numeric order (e.g., `02_run_PCA.R`, `03_build_RF_model.R`, etc.).

4. For HPC setup (RStudio + Slurm), and GitHub ⇄ RStudio integration, see:
  `docs/README.md`


## Citation & License

If you use these scripts or refer to our data, please cite:  
  > This repository for now. 

Licensed under MIT. See `LICENSE` for full terms.
