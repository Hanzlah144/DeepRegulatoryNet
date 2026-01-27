# DeepRegulatoryNet

A machine learning-based pipeline for reconstructing circRNA–miRNA–mRNA regulatory networks and performing integrative functional analyses from transcriptomic datasets.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![Platform: Linux](https://img.shields.io/badge/Platform-Linux-FF6600?logo=linux&logoColor=white)](https://github.com/<USERNAME>/<REPO>)
[![Platform: Windows](https://img.shields.io/badge/Platform-Windows-0078D6?logo=windows&logoColor=white)](https://github.com/<USERNAME>/<REPO>)
[![Platform: macOS](https://img.shields.io/badge/Platform-macOS-4CAF50?logo=apple&logoColor=white)](https://github.com/<USERNAME>/<REPO>)


## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Parameters](#parameters)
- [Input File Structure](#input-file-structure)
- [Usage](#usage)
- [Running with Test Data](#running-with-test-data)
- [License](#license)
- [Contact](#contact)

## Introduction

DeepRegulatoryNet is a comprehensive bioinformatics pipeline designed to investigate the complex regulatory relationships among circRNA-miRNA-mRNA networks. It employs machine learning models to predict circRNA-miRNA binding sites, constructs regulatory networks, and performs downstream analyses, including functional enrichment, protein-protein interaction mapping, and drug-gene interaction analysis to identify potential biomarkers and therapeutic targets.



## Features

- **Machine Learning Prediction**: Uses machine learning for accurate circRNA-miRNA binding site prediction
- **Regulatory Network Construction**: Builds comprehensive circRNA → miRNA → mRNA interaction networks
- **Multi-omics Integration**: Combines circRNA, miRNA, and mRNA expression data
- **Enrichment Analysis**: identifies pathways and biological processes enriched among overlapping target genes
- **PPI Network Analysis**: Constructs and analyzes protein-protein interaction networks
- **Drug-Gene Interactions**: Identifies potential drug targets through gene-drug association analysis
- **Comprehensive Reporting**: Generates detailed Excel reports and CSV outputs
- **Cross-platform compatibility**: works on ```Linux```, ```Windows```, and ```macOS``` systems
- **Debug and logging support**: Debug mode to track workflow execution and troubleshoot issues.

## Installation

### Dependencies

- Support Python 3.8+
- [pandas](https://pandas.pydata.org/)
- [NumPy](https://numpy.org/)
- [CatBoost](https://catboost.ai/)
- [scikit-learn](https://scikit-learn.org/stable/)
- [joblib](https://joblib.readthedocs.io/)
- [NetworkX](https://networkx.org/)
- [PyVis](https://pyvis.readthedocs.io/)
- [Matplotlib](https://matplotlib.org/)
- [Seaborn](https://seaborn.pydata.org/)
- [matplotlib-venn](https://pypi.org/project/matplotlib-venn/)
- [GSEApy](https://gseapy.readthedocs.io/)
- [Requests](https://requests.readthedocs.io/)
- [Beautiful Soup 4](https://www.crummy.com/software/BeautifulSoup/bs4/doc/)
- [urllib3](https://urllib3.readthedocs.io/)


### Using ```pip```

```bash
$ git clone https://github.com/Hanzlah144/DeepRagulatory-Net.git
```
Navigate into the ```DeepRegulatoryNet``` repository
```
$ cd DeepRegulatoryNet
```

Install all Python dependencies required to run ```DeepRegulatoryNet```
```
$ pip install -r requirements.txt
```

### Using ```conda```
Clone the Git Repository
```bash
$ git clone https://github.com/Hanzlah144/DeepRagulatory-Net.git
```
Navigate into the repository
```
$ cd DeepRegulatoryNet
```

Set up your environment and install dependencies
```
$ conda env create -f environment.yml
```

Activate the environment

```
$ Conda activate deepregulatorynet
```


Ensure all required model files are present in the `Model_Files/` directory:
   - `calibrated_catboost_site_type_model.pkl`
   - `label_encoder.pkl`
   - `robust_scaler.pkl`

## Parameters

The pipeline accepts the following command-line arguments:
| Parameter | Description | Requirement |
|-----------|-------------|------------|
| `--circ` | Path to a text file containing circRNA IDs, one per line | *Mandatory* |
| `--mirna` | Path to a text file containing miRNA IDs, one per line | *Mandatory* |
| `--deg` | Path to a text file containing differentially expressed gene (DEG) symbols, one per line | *Mandatory* |
| `--debug` | Enable debug logging to display detailed tracebacks and workflow messages | *Optional* |

## Input Files Structure

To run **DeepRegulatoryNet**, you need three input files:

1. circRNA IDs data  
2. miRNA IDs data 
3. DEG (Differentially Expressed Genes) data 

All input files should be in plain text format (`.txt`) with **one identifier per line**.

### circRNA Input Data File
Must contain human circRNA IDs, one per line. IDs should start with `hsa_circ_` following standard human circRNA nomenclature.

```
hsa_circ_000001
hsa_circ_000002
hsa_circ_000003
```

### miRNA Input Data File
Must contain human miRNA IDs, one per line. IDs should start with `hsa-miR` following miRBase human nomenclature.

```
hsa-miR-1
hsa-miR-21
hsa-miR-155
```

### DEG Input Data File
Must contain gene symbols, one per line. Gene symbols should follow ```HGNC``` standard nomenclature with the following structure; no specific prefix is required.

```
TP53
EGFR
MYC
BRCA1
```


## Usage
Runs the DeepRegulatoryNet pipeline on circRNA, miRNA, and DEG files to analyze regulatory interactions, with an optional `--debug` mode for detailed logging.


```bash
python DeepRegulatoryNet.py --circ <circRNA_file> --mirna <miRNA_file> --deg <DEG_file> [--debug]
```

## Running with Test Data

Executes the ```DeepRegulatoryNet``` using test data in the `Test_Data/` directory. use the following command:


```bash
python DeepRegulatoryNet.py --circ Test_Data/DEcircRNA.txt --mirna Test_Data/DEmiRNA.txt --deg Test_Data/DEG.txt
```

This will process the test datasets and generate output files in the `output/` directory. The test data includes:
- `DEcircRNA.txt`: Sample circRNA IDs
- `DEmiRNA.txt`: Sample miRNA IDs
- `DEG.txt`: Sample differentially expressed genes

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Reference
In process.


## Contact

For questions, issues, or contributions, please contact:

- **Maintainer**: [Your Name]
- **Email**: [your.email@example.com]

---
