# library(shiny)
# library(bulkAnalyseR)
# library(ggplot2)
# library(mixOmics)
# library(tidyverse)
# library(MSIToolKit)
# library(data.table)
# 
# new.expression.matrix <- data.table::fread('../../prep_datasets/IPF_new_expression_matrix.csv',nThread = 10)
# colnames(new.expression.matrix)=gsub('X','',colnames(new.expression.matrix))
# rownames(new.expression.matrix)=new.expression.matrix$V1
# new.expression.matrix = new.expression.matrix[,2:ncol(new.expression.matrix)]
# 
# new.metadata <- read.csv('../../prep_datasets/IPF_new_metadata.csv',row.names = 1)
# blue = read.csv('../../../IPF/New/Clustering/blue_spots.csv',header=F)
# red = read.csv('../../../IPF/New/Clustering/red_spots.csv',header=F)
# yellow = read.csv('../../../IPF/New/Clustering/yellow_spots.csv',header=F)
# new.metadata$bisecting_kmeans_clusters = ifelse(new.metadata$spot_id%in%blue$V1,'red',ifelse(new.metadata$spot_id%in%red$V1,'blue','yellow'))
# bulked <- create_bulk_exp(new.expression.matrix,new.metadata,sample.id.column = "Group",sample.wide.columns = c('Treatment','Timepoint','Rep'))
# # show the user the assigned identities and ask if they want to include everything or just the annotated ones
# bulk.expression.matrix = bulked$expression_matrix
# bulk.metadata = bulked$metadata
# bulk.metadata = bulk.metadata[,c(1,4,3,2)]
# anno = bulked$annotation_table
# only.annotated=T
# if (only.annotated){
#   anno = anno[anno$adduct!='NA',]
#   bulk.expression.matrix = bulk.expression.matrix[anno$m_z,]
#   new.expression.matrix = as.data.frame(new.expression.matrix)[,anno$m_z]
# }
# 
# library(shiny)
# library(bulkAnalyseR)
# library(ggplot2)
# library(mixOmics)
# organism='Mouse'
# ui <- fluidPage(
#   navbarPage(
#     'MSI toolkit',
#     theme = shinythemes::shinytheme('flatly'),
#     header = tags$head(tags$style('body {overflow-y: scroll;}')),
#     footer = bookmarkButton(),
#     tabPanel(title = 'Introduction',
#              tabsetPanel(
#                 IntroAnnopanelUI(id='Anno', bulk.metadata = bulk.metadata),
#                 IntroSpatialVisPanelUI(id='spatialVis', bulk.metadata = bulk.metadata, full.metadata = new.metadata),
#              )
#     ),
#     tabPanel(title = 'Pseudobulk Analysis',
#              tabsetPanel(
#                BulkQCpanelUI(id='BulkQC', bulk.metadata = bulk.metadata),
#                BulkDEpanelUI(id='BulkDE', bulk.metadata = bulk.metadata),
#                BulkDESummaryPanelUI(id='BulkSummaryDE', bulk.metadata = bulk.metadata),
#                BulkORAPanelUI(id='BulkORA', bulk.metadata = bulk.metadata)
#              )
#     ),
#     tabPanel(title = 'Region Analysis',
#              tabsetPanel(
#                RegionClusterPanelUI(id='RegionCluster', bulk.metadata = bulk.metadata, full.metadata = new.metadata),
#                RegionDEpanelUI(id='RegionDE', full.metadata = new.metadata, bulk.metadata),
#                BulkORAPanelUI(id='RegionORA', bulk.metadata = bulk.metadata)
#              )
#   )
# )
# )
# 
# server <- function(input, output, session) {
#   IntroAnnopanelServer(id='Anno', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno)
#   IntroSpatialVisPanelServer(id='spatialVis', bulk.metadata = bulk.metadata, full.expression.matrix = new.expression.matrix, full.metadata = new.metadata, anno = anno)
#   BulkQCpanelServer(id='BulkQC', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno)
#   bulkDEres <- BulkDEpanelServer(id='BulkDE', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno)
#   BulkDESummaryPanelServer(id='BulkSummaryDE', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = bulkDEres)
#   BulkORAPanelServer(id='BulkORA', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = bulkDEres)
#   clusters <- RegionClusterPanelServer(id='RegionCluster', full.expression.matrix = as.data.frame(t(new.expression.matrix)), full.metadata = new.metadata, bulk.metadata = bulk.metadata, anno = anno)
#   regionDEres <- RegionDEpanelServer(id='RegionDE', full.expression.matrix = new.expression.matrix, full.metadata = new.metadata, bulk.metadata = bulk.metadata, region.clusters = clusters, anno = anno)
#   BulkORAPanelServer(id='RegionORA', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = regionDEres)
# 
#   }
# 
# shinyApp(ui,server)
# 
