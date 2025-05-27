library(data.table)
library(Matrix)
library(CARD)
library(SingleCellExperiment) 
library(MuSiC) 
library(scatterpie)

var.deconvolu <- "celltype"
obj@meta.data <- obj@meta.data
################## deconvolution for each sample
results.card <- data.frame()
samples <- c("A1", "B1", "C1", "D1")
for(sample in samples){ 
    spot.temp <- obj@meta.data %>% filter(orig.ident == sample) %>% pull(barcode)
    spatialcount <- obj@assays$Spatial@counts[,spot.temp]
    spatialInfo <- data.frame(x = obj@images[[sample]]@coordinates$row, 
                              y = -obj@images[[sample]]@coordinates$col, 
                              row.names = rownames(obj@images[[sample]]@coordinates))[spot.temp,]

    CARD_obj = createCARDObject(
      sc_count = mat.count[, meta$cells],
      sc_meta = meta,
      spatial_count = spatialcount,
      spatial_location = spatialInfo,
      ct.varname = var.deconvolu,
      ct.select = unique(meta[,var.deconvolu]),
      sample.varname = NULL,
      minCountGene = 100,
      minCountSpot = 1) 
    CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

    #### save result
    proportion <- CARD_obj@Proportion_CARD
    spatial_location <- CARD_obj@spatial_location
    results.card <- rbind(results.card, cbind(spatial_location, proportion))
}
