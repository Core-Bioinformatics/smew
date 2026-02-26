# smew: Spatial metabolomics enhanced workflow 

## Usage

### Preprocessing step

To create a shiny app using *smew*, you need a processed (e.g. normalised) table of per-pixel peak intensities and a corresponding metadata table with spatial coordinates loaded in your workspace. 

The intensity matrix is expected to be a dataframe where rows correspond to pixels (named by pixel ID) and columns correspond to peaks (named by m/z value).

For example, your intensity matrix and metadata table may look something like:

|    |m/z 101.2345   |  m/z 102.85754 | m/z 303.35855   | m/z 344.48575  | m/z 321.38583  | m/z 112.28485 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---: |
| pixel_1| 5774.812  | 675.361  | 23.555  | 8444.958  | 777.234  | 20.332  |
| pixel_2| 8794.013  | 444.523  | 81.294  | 6775.393  | 899.284  | 10.275  |
| pixel_3| 6777.358  | 857.585  | 14.326  | 9468.367  | 747.385  | 24.521  |
| | ...  | ...  | ...  | ...  | ... | ... |

pixel_id    | Sample   | Treatment | Fibrotic |
| :---:        | :---:       | :---:   | :---: |
| pixel_1 | sample_1        | treatment_1 | normal |
| pixel_2 | sample_1        | treatment_1 | normal |
| pixel_3 | sample_1       | treatment_1 | fibrotic |
| pixel_4 | sample_2       | treatment_1 | fibrotic |
| pixel_5 | sample_2       | treatment_1 | normal |
| pixel_6 | sample_2       | treatment_1 | fibrotic |
| | ...  | ...  | ...  | 



The first column of the metadata table must match the row names of the intensity matrix, named *pixel_id*, followed by *x* and *y* columns containing spatial coordiantes. There must be another column Sample which describes the sample a pixel corresponds to. Other columns can contain sample-wide information (e.g. treatment group) and other metadata information containing individual pixel information. 

### Creating an app

Once your data is in the format above, you can create an app in just 1 line of code. In the case, suppose you have data from mouse in negative ion mode where the sample name is denoted by 'Sample_Name', you have 2 extra sample-wide columns called 'Timepoint' and 'Treatment' and 1 pixel-wise metadata column called 'Fibrotic'. In this case we only want to include peaks which can be annotated as known metabolites.

```{r}
library(smew) 

generateShinyApp(
  shiny.dir = 'MyShinyApp',
  intensity.matrix = my.intensity.matrix,
  metadata = my.metadata,
  sample.id.column = 'Sample_Name',
  sample.wide.columns = c('Timepoint','Treatment'),
  only.annotated = TRUE,
  organism = 'Mouse'
)

#run shiny app
shiny::runApp(shiny.dir)
```

## Installation guide

To use *smew*, you need R >= 4.0. Currently, *smew* can only be installed from GitHub, either by cloning the repository and using *devtools::install()* or using *devtools::install_github("Core-Bioinformatics/MSIToolkit")*. You need to make sure all dependencies are installed using the following:

### Required CRAN packages (use *install.packages()*) ###

* bsicons
* bslib
* ClustAssess
* data.table
* dendextend
* dplyr
* dunn.test
* DT
* ggplot2
* ggplotify
* ggnewscale
* ggrepel
* ggVennDiagram
* ggrastr
* grDevices
* harmony
* heatmaply
* htmlwidgets
* igraph
* jpeg
* Matrix
* matrixStats
* methods
* parallel
* patchwork
* pbapply
* plotly
* RColorBrewer
* RcppML
* reshape2
* scales
* shiny
* shinyjqui
* shinyjs
* shinyWidgets
* stats
* stringr
* tibble
* tidyr
* UpSetR
* utils
* viridis
* visNetwork

### Required Bioconductor packages (use *BiocManager::install()*) ###

* BayesSpace
* BiocSingular
* GENIE3
* mixOmics
* S4Vectors
* scater
* scran
* SingleCellExperiment

### Required packages to install from GitHub

To download plotly outputs to file, you may also need to run *webshot::install_phantomjs()*