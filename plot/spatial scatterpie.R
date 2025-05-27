############# define function
plot.spatial.pie <- function(data = plot.data, 
                             image_path = image_path, 
                             sample = sample, 
                             pie_scale = 0.8, 
                             columns){
    spatial_coord <- data.frame(obj@images[[sample]]@coordinates) %>%
            tibble::rownames_to_column("barcodeID") %>% dplyr::mutate(imagerow_scaled = imagerow *
            obj@images[[sample]]@scale.factors$lowres, imagecol_scaled = imagecol *
            obj@images[[sample]]@scale.factors$lowres) %>% dplyr::inner_join(data %>%
            tibble::rownames_to_column("barcodeID"), by = "barcodeID")
            
    img <- png::readPNG(image_path)
     
    img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1,
            "npc"), height = grid::unit(1, "npc"))
     
    #### barplot for percentage
    p1 <- scatterpie_pie <- ggplot2::ggplot() + 
        ggplot2::annotation_custom(grob = img_grob,
        xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
        # scatterpie::geom_scatterpie(data = spatial_coord, ggplot2::aes(x = imagecol_scaled,y = imagerow_scaled), cols = c("LI", "SI", "ST"), color = NA,alpha = 1, pie_scale = pie_scale) +
        scatterpie::geom_scatterpie(data = spatial_coord, ggplot2::aes(x = imagecol_scaled,y = imagerow_scaled), cols = columns, color = NA, alpha = 1, pie_scale = pie_scale) +
        ggplot2::scale_y_reverse() + ggplot2::ylim(nrow(img),
        0) + ggplot2::xlim(0, ncol(img)) + cowplot::theme_half_open(11,
        rel_small = 1) + ggplot2::theme_void() + ggplot2::coord_fixed(ratio = 1,
         xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") & 
    # scale_fill_manual(values = c("LI" = "#fa58d0", "SI" = "#3ADF00", "ST" = "#6666FF")) 
    # &
        scale_fill_manual(values = c("MES" = "#FFE4B5", "EPI" = "#E41A1C", "NC" = "#4DAF4A"))
    
    return(p1)
}

propo <- obj@meta.data[,c("x", "y", "MES", "EPI", "NC", "orig.ident")]

p1 <- plot.spatial.pie(data = plot.data, ### metadata containing the x and y coordinates
                       image_path = "/outs/spatial/tissue_lowres_image.png", 
                       sample = "A1", 
                       pie_scale = 1.7, 
                       columns = c("MES", "EPI", "NC"))

options(repr.plot.width = 10, repr.plot.height = 5)
pdf("./spatial.region.final.proportion.pdf", width = 10, height = 5)
p1 
dev.off()
