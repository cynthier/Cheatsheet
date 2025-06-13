library(SeuratData)
library(SeuratDisk)
library(Seurat)

objs.sub$pseudotime <- pse[colnames(objs.sub), "pseudotime"]

#################### remove cells
meta <- objs.sub@meta.data
ggplot(meta, aes(x = pseudotime, fill = group)) + 
    geom_density(alpha = 0.5) + 
    theme_minimal(base_size = 16) +
    # scale_fill_manual(values = colors) + 
    geom_vline(xintercept = 65) + 
    geom_vline(xintercept = 2) 

objs.sub <- subset(objs.sub, pseudotime < 65 & pseudotime > 2)

####################
#### function
process_G2G <- function(obj, pseudotime, h5d.name, exp.filename){ 
     obj$time <- (obj@meta.data[, pseudotime]-min(obj@meta.data[, pseudotime]))/(max(obj@meta.data[, pseudotime])-min(obj@meta.data[, pseudotime]))
    
    # print(max(obj@meta.data[, "time"]))
    # print(min(obj@meta.data[, "time"]))
    
    exp.data <- GetAssayData(object = obj, assay = "RNA", layer = "data")
    print(dim(exp.data))
    print(dim(obj))
    
    write.table(exp.data, file = paste0(exp.filename, ".csv"), quote = F, sep = ",", row.names = T)
    SaveH5Seurat(obj, file = paste0(h5d.name, ".h5Seurat"), overwrite = TRUE)
    Convert(paste0(h5d.name, ".h5Seurat"), dest = "h5ad") 
}


#### need to create obj again
objs.sub$cell <- colnames(objs.sub)
obj.temp <- CreateSeuratObject(count = GetAssayData(objs.sub, layer = "counts"), meta.data = objs.sub@meta.data, min.cells = 0, min.features = 0)

#### all one cell to max pseudotime, ensure the aligment in control and mutant
objs.sub.con <- subset(obj.temp, group == "Control")
objs.sub.ko <- subset(obj.temp, group == "Vcl cKO")


#### convert control data
objs.sub.con <- NormalizeData(objs.sub.con)
objs.sub.con <- ScaleData(objs.sub.con)
objs.sub.con[["RNA3"]] <- as(object = objs.sub.con[["RNA"]], Class = "Assay")
DefaultAssay(objs.sub.con) <- "RNA3"
objs.sub.con[["RNA"]] <- NULL
objs.sub.con <- RenameAssays(object = objs.sub.con, RNA3 = 'RNA')
process_G2G(obj = objs.sub.con, pseudotime = "pseudotime", h5d.name = "1_con_BranchA", exp.filename = "1_exp_con_BranchA") 

#### convert mutant data
objs.sub.ko <- NormalizeData(objs.sub.ko)
objs.sub.ko <- ScaleData(objs.sub.ko)
objs.sub.ko[["RNA3"]] <- as(object = objs.sub.ko[["RNA"]], Class = "Assay")
DefaultAssay(objs.sub.ko) <- "RNA3"
objs.sub.ko[["RNA"]] <- NULL
objs.sub.ko <- RenameAssays(object = objs.sub.ko, RNA3 = 'RNA')
process_G2G(obj = objs.sub.ko, pseudotime = "pseudotime", h5d.name = "1_ko_BranchA", exp.filename = "1_exp_ko_BranchA")
