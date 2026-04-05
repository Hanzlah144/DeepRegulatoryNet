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
- [License and Issues](#license-and-issues)
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

### Prerequisites

- Python 3.9+
- Internet access for external APIs used by the pipeline (NIH CircInteractome, miRDB, STRING-DB, ChEMBL)
- The `trained models/` directory must contain the pre-trained model artifacts:
  - `catboost_model.pkl`
  - `label_encoder.pkl`
  - `scaler.pkl`

### Dependencies

The pipeline depends on the following Python packages:

- pandas>=1.5.0
- numpy>=1.23.0
- catboost>=1.2.0
- scikit-learn==1.5.1
- joblib>=1.2.0
- networkx>=3.0
- pyvis>=0.3.0
- matplotlib>=3.6.0
- seaborn>=0.12.0
- matplotlib-venn>=0.11.9
- gseapy>=1.0.0
- requests>=2.28.0
- requests-cache>=1.0.0
- beautifulsoup4>=4.11.0
- urllib3>=1.26.0
- lxml>=4.9.0
- openpyxl>=3.0.0
- chembl-webresource-client
- ipython>=8.0.0

### Installation Using `pip`

```bash
$ git clone https://github.com/Hanzlah144/DeepRegulatoryNet.git
$ cd DeepRegulatoryNet
$ pip install -r requirements.txt
```

### Installation Using `conda`

```bash
$ git clone https://github.com/Hanzlah144/DeepRegulatoryNet.git
$ cd DeepRegulatoryNet
$ conda env create -f environment.yml
$ conda activate deepregulatorynet
```

### Handling Trained Model Files

If model files are missing or corrupted in `trained models/`, retrieve them with Git LFS:

```bash
$ cd 'trained models'
$ git lfs install
$ git lfs pull
```

**Note**: The project uses Git LFS for large model files. If you cloned the repository before Git LFS was configured, install Git LFS first and then pull the model files.

## Parameters

The pipeline accepts the following command-line arguments:
| Parameter | Description | Requirement |
|-----------|-------------|------------|
| `--circ` | Path to a text file containing circRNA IDs, one per line | *Mandatory* |
| `--mirna` | Path to a text file containing miRNA IDs, one per line | *Mandatory* |
| `--deg` | Path to a text file containing differentially expressed gene (DEG) symbols, one per line | *Mandatory* |
| `--debug` | Enable debug logging and tracebacks | *Optional* |
| `--mode` | Pipeline mode: `quick` or `full` | *Optional* |
| `--max_genes` | Maximum number of hub genes for drug-gene analysis when running in quick mode | *Optional* |

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

DeepRegulatoryNet processes circRNA, miRNA, and DEG input files to reconstruct regulatory networks and perform functional analyses. The pipeline runs in two modes: `full` (default) and `quick`.

### Full Mode (Default)
Performs complete analysis including all steps: binding site prediction, network construction, enrichment, PPI, and drug-gene interactions.

```bash
python DeepRegulatoryNet.py --circ <circRNA_file> --mirna <miRNA_file> --deg <DEG_file> [--debug]
```

### Quick Mode
Limits drug-gene analysis to a subset of hub genes for faster execution.

```bash
python DeepRegulatoryNet.py --circ <circRNA_file> --mirna <miRNA_file> --deg <DEG_file> --mode quick --max_genes 50 [--debug]
```

### Pipeline Workflow
1. **Binding Site Prediction**: Uses ML to predict circRNA-miRNA interactions from NIH CircInteractome data.
2. **Network Construction**: Builds tripartite circRNA→miRNA→mRNA networks.
3. **Enrichment Analysis**: Identifies enriched pathways and GO terms for overlapping genes.
4. **PPI Analysis**: Constructs and analyzes protein-protein interaction networks.
5. **Drug-Gene Interactions**: Maps hub genes to ChEMBL drug targets.

Outputs are saved in the `output/` directory, including CSV files, Excel reports, GraphML networks, and visualizations. A `pipeline.log` file records execution details.

## Test DeepRegulatoryNet with Example Data

Executes the ```DeepRegulatoryNet``` using test data in the `examples/` directory. use the following command:


```bash
python DeepRegulatoryNet.py --circ examples/DEcircRNA.txt --mirna examples/DEmiRNA.txt --deg examples/DEG.txt --debug
```

This will process the test datasets and generate output files in the `output/` directory. The test data includes:
- `DEcircRNA.txt`: Sample circRNA IDs
- `DEmiRNA.txt`: Sample miRNA IDs
- `DEG.txt`: Sample differentially expressed genes

## License and Issues

This DeepRegulatoryNet is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Submit issues or contributions via [GitHub Issues](https://github.com/Hanzlah144/DeepRegulatoryNet/issues).

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
