# DeepRegulatoryNet

A machine learning-based pipeline for reconstructing circRNA–miRNA–mRNA regulatory networks and performing integrative functional analyses from transcriptomic datasets.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![Platform: Linux](https://img.shields.io/badge/Platform-Linux-FF6600?logo=linux&logoColor=white)](https://github.com/<USERNAME>/<REPO>)
[![Platform: Windows](https://img.shields.io/badge/Platform-Windows-0078D6?logo=windows&logoColor=white)](https://github.com/<USERNAME>/<REPO>)
[![Platform: macOS](https://img.shields.io/badge/Platform-macOS-4CAF50?logo=apple&logoColor=white)](https://github.com/<USERNAME>/<REPO>)


## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [DeepRegulatoryNet Installation](#installation)
- [Parameters](#parameters)
- [Input File Structure](#input-files-structure)
- [Usage](#usage)
- [Test DeepRegulatoryNet with Example Data](#test-deepregulatorynet-with-example-data)
- [License](#license)
- [Authors and Contacts](#authors-and-contacts)

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

- Support Python 3.9+
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


### Installation Using ```pip```

```bash
$ git clone https://github.com/Hanzlah144/DeepRegulatoryNet.git
```
Navigate into the ```DeepRegulatoryNet``` repository
```
$ cd DeepRegulatoryNet
```

Install all Python dependencies required to run ```DeepRegulatoryNet```
```
$ pip install -r requirements.txt
```

### Installation Using ```conda```
Clone the Git Repository
```bash
$ git clone https://github.com/Hanzlah144/DeepRegulatoryNet.git
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

### Handling Trained Model Files

If you encounter issues with model files in the `trained models/` directory (such as missing or corrupted files), follow these steps to retrieve them using Git LFS (Large File Storage):

```bash

$ cd trained\ models

$ git lfs install

$ git lfs pull
```

**Note**: This project uses Git LFS to manage large model files efficiently. If you cloned the repository before Git LFS was initialized, ensure you have Git LFS installed on your system before pulling the model files. For more information, visit [Git LFS Documentation](https://git-lfs.github.com/).

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

## Test DeepRegulatoryNet with Example Data

Executes the ```DeepRegulatoryNet``` using test data in the `examples/` directory. use the following command:


```bash
python DeepRegulatoryNet.py --circ examples/DEcircRNA.txt --mirna examples/DEmiRNA.txt --deg examples/DEG.txt --debug
```

This will process the test datasets and generate output files in the `output/` directory. The test data includes:
- `DEcircRNA.txt`: Sample circRNA IDs
- `DEmiRNA.txt`: Sample miRNA IDs
- `DEG.txt`: Sample differentially expressed genes

## License

This DeepRegulatoryNet is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Reference
In process.


## Authors and Contacts

**Muhammad Hanzlah Shah Khalid**  
*Integrative Omics and Molecular Modeling Laboratory, Department of Bioinformatics and Biotechnology, Government College University Faisalabad (GCUF), Faisalabad, 38000, Pakistan*  
Email: [nk6692381@gmail.com](mailto:nk6692381@gmail.com)

**Rana Sheraz Ahmad**  
*Integrative Omics and Molecular Modeling Laboratory, Department of Bioinformatics and Biotechnology, Government College University Faisalabad (GCUF), Faisalabad, 38000, Pakistan*  
Email: [ranasheraz.202101902@gcuf.edu.pk](mailto:ranasheraz.202101902@gcuf.edu.pk)

**Fatima Noor**  
*Institute of Molecular Biology and Biotechnology (IMBB), The University of Lahore, Pakistan*  
Email: [fatimanoor1122@yahoo.com](mailto:fatimanoor1122@yahoo.com)

**Faheem Ahmed Khan**  
*Stem Cell and Cancer Research Indonesia; Department of Transfusion Medicine and Clinical Microbiology, Faculty of Allied Health Sciences, Chulalongkorn University, Bangkok, 10330, Thailand*  
Email: [faheemgenetics@yahoo.com](mailto:faheemgenetics@yahoo.com)

**Muhammad Tahir ul Qamar**  
*Integrative Omics and Molecular Modeling Laboratory, Department of Bioinformatics and Biotechnology, Government College University Faisalabad (GCUF), Faisalabad, 38000, Pakistan*  
Email: [m.tahirulqamar@hotmail.com](mailto:m.tahirulqamar@hotmail.com)
