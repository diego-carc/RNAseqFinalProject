# Differential expression analysis

## Introduction 
This project is part of the course *Intro RNA-seq LCG-UNAM 2024*, which is an introduction to RANseq data differential expression analysis tools available at *Bioconductor*.
In the following report, I present the results of a differential expression analysis using the dataset *Transcriptomic profiling reveals distinct modes of aging in the kidney* (SRP165875).
This dataset was generated as part of a study amied at revealing the molecular mechanisms under the aging in the kidneys, obtaining 188 samples of mouse kidney at different life stages, that is, 6, 12 and 18 months of birth.

## Data retrieval
To retreive the data, I created a script calling the function `create_rse_manual` available at   `recount3`.

## Data preprocessing
Since filtering by expression and normalization are crucial steps in a differential expression analysis, I developed a script applying this processes to the rae data. Then I compared the data ditribution before and after the abovementioned transformations.
Data pre filtering:
![Raw data distribution](figures/rawDataDist.png)
Data post filtering:
![Normalized data](figures/filtDataDist.png)

## Variable analysis
Before making the differential expression analysis, I verified the relations between the variables in the dataset.
First, I used a boxplot to visualize the effect of the age of the mouses in the overall alignment:
!(figures/ageBoxPlot.png)
## Differential expression analysis
Once I filtered the data, I started to making the differential expression analysis.
To this end, I created the script *diffExp.R*, in which I 
## Conclusion
