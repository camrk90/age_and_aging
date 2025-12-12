#!/usr/bin/env Rscript

library(flowCore)
library(umap)
library(ggplot2)
library(monocle3)
load("/home/ckelsey4/cayo_flow_data/flow_map.RData")

flow_path="/data/CEM/smacklab/cayo_flow_cytometry/" 

fcs = read.FCS(paste(flow_path, "Flow_Data_2019/Cayo Santiago Week 3 (Oct-Nov 28-1)/2019-11-01/04Z2019-11-01.0001.fcs", 
                     sep = ""),emptyValue=FALSE,transformation=FALSE)

# Call it "e.matrix" i.e., expression matrix so that I can copy-paste code, but dimensions are flow gates instead of expressed genes
e.matrix = fcs@exprs

# Instead of cell IDs, the columns are "events"
rownames(e.matrix) = paste0('event',formatC(seq_len(nrow(e.matrix)),width=8,flag=0))
e.matrix<- e.matrix[, 5:14]

# Drop forward scatter < 225 and and values for FSC or SSC that are 1000 (the maximum)
e.gated = e.matrix[e.matrix[,'FSC-A'] >= 225 & e.matrix[,'FSC-A'] < 1000 & e.matrix[,'SSC-A'] < 1000,]

# Scale data
e.scaled = scale(e.gated)[,grep('^FSC-A$',colnames(e.gated)):ncol(e.gated)]

# Treat this as mock single-cell RNAseq data
expression.data = t(e.scaled)

options(stringsAsFactors=FALSE)

cell.metadata = data.frame(cell_name = colnames(expression.data)[1:10000])
rownames(cell.metadata) = colnames(expression.data)[1:10000]

gene.metadata = data.frame(gene_short_name = rownames(expression.data))
rownames(gene.metadata) = rownames(expression.data)

cds = new_cell_data_set(exp(expression.data[,1:10000]),cell.metadata,gene.metadata)

cds = preprocess_cds(cds)
cds = reduce_dimension(cds)
cds = cluster_cells(cds)
plot_cells(cds,color_cells_by='cluster')
plot_cells(cds,color_cells_by='partition')

channels<- colnames(e.scaled)
plot_cells(cds, genes=channels[1:10],scale_to_range = F)


save.image("/home/ckelsey4/cayo_flow_data/flow_map.RData")



