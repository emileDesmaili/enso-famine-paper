# El Niño Amplified Food Security in Early Modern Europe

This repository contains the data and code used to generate the figures and results in the manuscript: **"El Niño amplified food security in Early Modern Europe"** (Esmaili et al., 2025, in review).


## Repository Structure

```
enso-famine-paper/
├── data/                    # Raw and processed data used in the analysis
├── processed_data/          # Final processed datasets used for analysis
├── preprocessing codes/     # Scripts for data preprocessing and cleaning
├── analysis notebooks/      # Jupyter notebooks for analysis and figure generation
├── figures/                 # Generated figures and visualizations
├── environment.yml          # Conda environment specification
└── README.md               # This file
```

### Directory Details

- **`data/`**: Contains the raw and intermediate data files used in the analysis. Note that ModERA and OWDA datasets are not included here as they can be publicly accessed online (see Data Sources section).

- **`processed_data/`**: Contains the final processed datasets that are directly used in the analysis notebooks.

- **`preprocessing codes/`**: Python scripts and notebooks used to clean, process, and prepare the raw data for analysis.

- **`analysis notebooks/`**: Jupyter notebooks containing the main analysis code, statistical tests, and figure generation scripts.

- **`figures/`**: Output directory containing all generated figures, plots, and visualizations used in the manuscript.

## Setup Instructions

### Prerequisites
- Python 3.10+
- Conda or Miniconda

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/emileDesmaili/enso-famine-paper.git
   cd enso-famine-paper
   ```

2. **Create and activate the conda environment:**
   ```bash
   conda env create -f environment.yml
   conda activate famine-enso
   ```




*Note: External datasets (ModERA and OWDA) are not included in this repository but can be accessed from their respective public sources.*


## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{esmaili2025enso,
  title        = {El Niño amplified food security in Early Modern Europe},
  author       = {Esmaili, Emile and Puma, Michael J. and Ludlow, Francis and Jobbová, Eva and Kumar, Janavi and Holm, Poul and Ljungqvist, Fredrik Charpentier and Matthews, John Alphonsus and Dahl, Johannes Rom and Seim, Andrea},
  year         = {2025},
  note         = {in review}

```

## Contributing

This repository is maintained for research reproducibility. For questions or issues, please open an issue or contact the corresponding author.

## License

MIT License

## Contact

For questions about this research or code, please contact:
- Emile Esmaili: ee4561@princeton.edu
- Repository maintainer: Emile Esmaili

---

**Note**: This repository is associated with a manuscript currently under review. Some details may be updated following peer review and publication.