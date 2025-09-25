
This repository contains the data and codes used to generate the figures and results in the manuscript: **'El Nino amplified food security in Early Modern Europe'** (Esmaili et al., 2025, *in review*).

The repository is structured as follows:
- `data/`: contains the raw and data used in the analysis, except for the ModERA and OWDA datasets that can be publicly found online.
- `processed_data/`: contains the processed data used in the analysis
- `preprocessing codes/`: contains the code used to preprocess the raw data
- `analysis notebooks/`: contains the code used to perform the analysis and generate the figures
- `figures/`: contains the figures generated from the analysis

To create the correct conda environment, please use the provided `environment.yml` file. You can create the environment by running the following command in your terminal:

```bash
conda env create -f environment.yml
```