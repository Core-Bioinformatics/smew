library(shiny)
library(bulkAnalyseR)
library(ggplot2)
library(mixOmics)
library(tidyverse)

new.expression.matrix <- data.table::fread('../prep_datasets/IPF_new_expression_matrix.csv',nThread = 10)
rownames(new.expression.matrix)=new.expression.matrix$V1
new.expression.matrix = new.expression.matrix[,2:ncol(new.expression.matrix)]
new.expression.matrix = as.data.frame(new.expression.matrix)
new.metadata <- read.csv('../prep_datasets/IPF_new_metadata.csv',row.names = 1)

bulked <- create_bulk_exp(new.expression.matrix,new.metadata,sample.id.column = "Group",sample.wide.columns = c('Treatment','Timepoint','Rep'))
# show the user the assigned identities and ask if they want to include everything or just the annotated ones
bulk.expression.matrix = bulked$expression_matrix
bulk.metadata = bulked$metadata
anno = bulked$annotation_table
only.annotated=T
if (only.annotated){
  anno = anno[anno$adduct!='NA',]
  bulk.expression.matrix = bulk.expression.matrix[anno$m_z,]
}

library(shiny)
library(bulkAnalyseR)
library(ggplot2)
library(mixOmics)
organism='Mouse'
ui <- fluidPage(
  navbarPage(
    'MSI toolkit',
    theme = shinythemes::shinytheme('flatly'),
    header = tags$head(tags$style('body {overflow-y: scroll;}')),
    footer = bookmarkButton(),

    tabPanel(title = 'Pseudobulk Analysis',
             tabsetPanel(
               # map against all known metabolites and create anno table with KEGG ID, Name and m/z
               BulkQCpanelUI(id='BulkQC', bulk.metadata = bulk.metadata),
               BulkDEpanelUI(id='BulkDE', bulk.metadata = bulk.metadata),
               BulkDESummaryPanelUI(id='BulkSummaryDE', bulk.metadata = bulk.metadata),
               BulkORAPanelUI(id='BulkORA', bulk.metadata = bulk.metadata)
               )
             )
  )
)

server <- function(input, output, session) {
  BulkQCpanelServer(id='BulkQC', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno)
  bulkDEres <- BulkDEpanelServer(id='BulkDE', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno)
  BulkDESummaryPanelServer(id='BulkSummaryDE', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = bulkDEres)
  BulkORAPanelServer(id='BulkORA', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = bulkDEres)
}

shinyApp(ui,server)
