[README.md](https://github.com/user-attachments/files/33481240/README.md)
# Association Between Tobramycin MIC and BPC in *Pseudomonas aeruginosa*

## Project description

This project investigates whether *Pseudomonas aeruginosa* isolates with higher tobramycin minimum inhibitory concentration (MIC) values also tend to have higher biofilm prevention concentration (BPC) values.

MIC represents the lowest antibiotic concentration that prevents visible growth of planktonic bacterial cells, whereas BPC represents the concentration that prevents biofilm formation.

## Research question

Do *Pseudomonas aeruginosa* isolates with higher tobramycin MIC values also tend to have higher BPC values?

## Dataset

The publicly available data were obtained from the GitHub repository associated with the original study:

https://github.com/gdewael/biofilm-amr/tree/main/data/isolates

The analysis uses the following files:

- `MIC_IMC.csv`
- `BPC_IMC.csv`

The original files contained 60 rows representing technical replicates of 30 bacterial isolates. The technical replicates were checked and combined to produce one representative MIC value and one representative BPC value for each isolate. The two datasets were then matched using isolate identifiers, resulting in 30 matched MIC–BPC pairs.

## Analysis

Data preparation and statistical analysis were performed in R. The analysis included:

- Checking the technical replicates
- Combining repeated measurements
- Matching MIC and BPC measurements by isolate ID
- Checking for missing values and duplicate isolate IDs
- Calculating descriptive statistics
- Calculating Spearman’s rank correlation
- Performing 10,000 paired bootstrap resamples
- Calculating the bootstrap standard error
- Calculating a 95% BCa confidence interval

The final figures were generated using Python and Matplotlib.

## Requirements

### R

The R analysis requires the `boot` package. Install it once by running:

```r
install.packages("boot")
```

### Python

The Python visualisation requires `pandas` and `matplotlib`. Install them by running:

```bash
python3 -m pip install pandas matplotlib
```

## Project structure

```text
PROJECT/
├── README.md
├── 01_prepare_data.r
├── 02_analysis.r
├── 03_correlation.r
├── 04_visualisation.py
├── MIC_IMC.csv
├── BPC_IMC.csv
├── figure1_MIC_BPC_association.png
├── figure2_bootstrap_uncertainty.png
└── report.pdf
```

## How to run the analysis

1. Download or clone this repository.
2. Open the project folder in Visual Studio Code.
3. Make sure that the scripts and input CSV files are in the same project folder.
4. Open an R terminal and run the scripts in the following order:

```r
source("01_prepare_data.r")
source("02_analysis.r")
source("03_correlation.r")
```

## How to generate the figures

After completing the R analysis, run the Python visualisation script from the VS Code terminal:

```bash
python3 04_visualisation.py
```

This script generates:

- `figure1_MIC_BPC_association.png`
- `figure2_bootstrap_uncertainty.png`

## Main results

The cleaned dataset contained 30 matched isolates.

- Median MIC: 2 µg/mL
- Median BPC: 16 µg/mL
- Spearman correlation: 0.899
- Bootstrap standard error: 0.040
- 95% BCa confidence interval: 0.777–0.948

The results indicate a strong positive association between tobramycin MIC and BPC. Isolates with relatively high MIC values generally also had relatively high BPC values. This association does not demonstrate causation.

