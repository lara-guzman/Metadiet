#Readme

    #Index
      #Libraries and references
      #1.Databases and list
      #2.PCA (Exploratory)
      #3.PLSDA-MULTIBLOCK (MixOmics-DIABLO MODEL)
      #4. Univariate analysis - Correlations


#Libraries and references
#mixOmics
#Guide:  https://mixomicsteam.github.io/mixOmics-Vignette/id_01.html
# Chapter N-Integration
#General conditions
#Type of data:e.g. microarray, mass spectrometry-based proteomics 
        #and metabolomics) or sequenced-based count data (RNA-seq, 16S, shotgun metagenomics)
#Normalization. The package does not handle normalisation. But it is necessary
        #Normalization in this case is CLRs for all database
#Prefiltering. It recommend pre-filtering the data to less than 10K predictor variables per data set
#Data format. Use matrix decomposition techniques
        # data frames have n observations or samples in rows 
        # P predictors or variables (e.g. genes, proteins, OTUs) in columns.

#Installation.
## install BiocManager if not installed
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      BiocManager::valid()
      library(BiocManager)

      ## install mixOmics
      BiocManager::install('mixOmics')
      ABiocManager::install("mixOmicsTeam/mixOmics")
      install_github("mixOmicsTeam/mixOmics")
      # https://www.bioconductor.org/packages/release/bioc/html/mixOmics.html
            #Libraries     
            library(mixOmics)
            library(BiocParallel)
            library(mixOmics)                                           
            library(devtools)
            library(ggplot2)
            library(lattice)
            library(MASS)
            library(rgl)
            library(igraph)
            library(ellipse)
            library(corpcor)
            library(RColorBrewer)
            library(plyr)
            library(parallel)
            library(dplyr)
            library(reshape2)
            library(methods)
            library(matrixStats)
            library(rARPACK)
            library(gridExtra)
            library(tidyr)
            library(vegan)
            library(factoextra) #Libreria principal
            library(ggplot2) #para graficar
            library(FactoMineR) #para el PCA
            library(readr)
            library(openxlsx)
            library(tibble)
            library(ggplot2)
            library(tidyr)
            library(dplyr)
            library(tidyverse)
            library(ggbreak)
            library(yulab.utils)
            library(aplot)
            library(ggfun)
            library(ggplotify)
            library(gridGraphics)
            library(ggplot2)
            library(openxlsx)
            library(RColorBrewer)
            library(viridis) 
            library(ggplot2)
            library(dplyr)
            library(ggbreak)
            library(viridis)
            library(caret)
      
      
      
      
#Directory
      setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet")
 
      
#1.Databases and list 
      
    #1.1. Databases
      setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/RAW")
      #Import data set from RAW in the global environment of R (Encoding "Automatic" and Heading "YES")
        #B_OTU_Plants_CLR
        #E_OTU_Microbiota_CLR
        #G_Metabolome_IsoMS
    
    #1.2 List
        Plants <- list (B_OTU_Plants_CLR[c(-1,-2)])
        Metabolites <- list (scale (G_Metabolome_IsoMS[c(-1,-2)])) #Metabolome database no compositional;beeter correaltions with metadiet
        Microbiota <- list (E_OTU_Microbiota_CLR[c(-1,-2)])
        Groups <- list (Group)

    #1.3 Plants classified by nutrition characteristics
        setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/RAW")
        nutri_class <- read.xlsx("D_tRNL_Nutritional_classification.xlsx")  # Nutritional classification
        # Read the file with CLR values
        otu_clr_data <- read_delim("B_OTU_Plants_CLR.csv", delim = ";", locale = locale(decimal_mark = ",", grouping_mark = "."))
        otu_clr_data <- read_delim("B_OTU_Plants_CLR.csv", delim = ";", 
                                   locale = locale(decimal_mark = ",", grouping_mark = "."), 
                                   col_types = cols(.default = "c"))  # Read all as text
        # Save the original order
        otu_clr_data <- otu_clr_data %>%
          mutate(Original_Order = row_number())  # New column to store the original order
        # Convert ID column to row names
        otu_clr_data <- otu_clr_data %>%
          column_to_rownames(var = "ID")
        # Convert values to numeric
        otu_clr_data_numeric <- otu_clr_data %>%
          mutate(across(-Original_Order, ~ as.numeric(gsub(",", ".", .))))  # Convert all except Original_Order
        # Filter taxa names that are present in the nutritional classification
        taxa_names <- colnames(otu_clr_data_numeric)
        nutri_class <- nutri_class %>%
          filter(Taxa %in% taxa_names)
        # Long format and join with classification
        otu_clr_long <- otu_clr_data_numeric %>%
          rownames_to_column(var = "Sample") %>%
          pivot_longer(cols = -c(Sample, Original_Order), names_to = "Taxa", values_to = "CLR_Value") %>%
          left_join(nutri_class, by = "Taxa")
        # Average CLR by nutritional classification
        otu_clr_summary <- otu_clr_long %>%
          group_by(Sample, Original_Order, Nutritional_Classification) %>%
          summarise(Mean_CLR = mean(CLR_Value, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = Nutritional_Classification, values_from = Mean_CLR)
        # Sort according to the original order from the CSV file
        otu_clr_summary <- otu_clr_summary %>%
          arrange(Original_Order) %>%
          select(-Original_Order)  # Remove auxiliary column
        # Write to CSV
        write_csv(otu_clr_summary, "B_OTU_Plants_CLR_nutritional_class.csv")
        
        
              #1.3.1 PLOT OF Plants classified by nutrition characteristics
              
              # Remove NA and filter out negative values
              otu_clr_summary_clean <- otu_clr_summary %>%
                select(-`NA`)
              otu_clr_long <- otu_clr_summary_clean %>%
                pivot_longer(cols = -Sample, names_to = "Nutritional_Classification", values_to = "Mean_CLR") %>%
                filter(Mean_CLR >= 0)
              # Order samples
              otu_clr_long <- otu_clr_long %>%
                mutate(Sample = factor(Sample, levels = unique(otu_clr_summary_clean$Sample)))
              # Calculate the number of observations per group
              label_data <- otu_clr_long %>%
                group_by(Nutritional_Classification) %>%
                summarise(n = n(), .groups = "drop")
              taxa_count <- nutri_class %>%
                group_by(Nutritional_Classification) %>%
                summarise(n_taxa = n(), .groups = "drop")
              ggplot(otu_clr_long, aes(x = Nutritional_Classification, y = Mean_CLR, fill = Nutritional_Classification)) +
                # Violin plot without fill or border
                geom_violin(trim = FALSE, color = "gray", fill = NA) +
                # Individual points
                geom_jitter(width = 0.2, alpha = 0.7, size = 1, color = "black") +
                # Group mean point
                stat_summary(fun = mean, geom = "crossbar", color = "black", size = 0.5) +
                # Labels for number of observations (above)
                geom_text(data = label_data, aes(x = Nutritional_Classification, 
                                                 y = max(otu_clr_long$Mean_CLR) + 0.1, 
                                                 label = paste0("n = ", n)),
                          inherit.aes = FALSE, size = 4, angle = 90, fontface = "italic") +
                # Labels for number of taxa at the bottom
                geom_text(data = taxa_count, aes(x = Nutritional_Classification, 
                                                 y = -0.1, 
                                                 label = paste0(n_taxa, " taxa")),
                          inherit.aes = FALSE, size = 3, fontface = "italic") +
                scale_fill_brewer(palette = "Set3") +
                labs(
                  title = "",
                  x = "Nutritional Classification",
                  y = "Mean CLR"
                ) +
                theme_minimal(base_size = 14) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                guides(fill = "none")
              
      #1.4 Microbitoa classified by family and plot

              
            #1.4.1 Top 5 Families with Highest Mean CLR        
                  # 1. Read CSV files
                  otu_microbiota <- read_csv2("E_OTU_Microbiota_CLR.csv")
                  taxonomy <- read_csv2("F_Taxonomy_Microbiota.csv")
                  
                  # 2. Ensure that all "OtuXXXX" columns are numeric
                  otu_microbiota_clean <- otu_microbiota %>%
                    mutate(across(starts_with("Otu"), ~ as.numeric(.)))
                  
                  # 3. Convert the matrix to long format
                  otu_long <- otu_microbiota_clean %>%
                    pivot_longer(cols = starts_with("Otu"), names_to = "OTU", values_to = "CLR_Value") %>%
                    rename(Sample = ID)
                  
                  # 4. Join with taxonomy table to get Family
                  otu_long_tax <- otu_long %>%
                    left_join(taxonomy %>% select(OTU, Family), by = "OTU")
                  
                  write.xlsx(otu_long_tax, file = "otu_long_tax.xlsx")
                  
                  
                  # 5. Calculate mean CLR per Family
                  top5_familias <- otu_long_tax %>%
                    group_by(Family) %>%
                    summarise(mean_CLR = mean(CLR_Value, na.rm = TRUE)) %>%
                    arrange(desc(mean_CLR)) %>%
                    slice_head(n = 5) %>%
                    pull(Family)
                  
                  # 6. Filter only top 5 families
                  otu_top5 <- otu_long_tax %>%
                    filter(Family %in% top5_familias)
                  
                  # 7. Count number of OTUs per family
                  taxa_count <- taxonomy %>%
                    filter(Family %in% top5_familias) %>%
                    group_by(Family) %>%
                    summarise(n_taxa = n(), .groups = "drop")
                  
                  # 8. Create plot
                  ggplot(otu_top5, aes(x = Family, y = CLR_Value, fill = Family)) +
                    geom_violin(trim = FALSE, color = "gray", fill = NA) +
                    geom_jitter(width = 0.2, alpha = 0.7, size = 1, color = "black") +
                    stat_summary(fun = mean, geom = "crossbar", color = "black", size = 0.5) +
                    geom_text(data = taxa_count, aes(x = Family, y = max(otu_top5$CLR_Value, na.rm = TRUE) + 0.2,
                                                     label = paste0("n = ", n_taxa)),
                              inherit.aes = FALSE, size = 3, angle = 90, fontface = "italic") +
                    scale_fill_brewer(palette = "Set3") +
                    labs(
                      x = "Family",
                      y = "CLR Value",
                      title = "Top 5 Families with Highest Mean CLR"
                    ) +
                    theme_minimal(base_size = 14) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    guides(fill = "none")
              
              #1.4.2 Top 5 Families with most quantities of OTUs 

                    otu_microbiota <- read_csv2("E_OTU_Microbiota_CLR.csv")
                    taxonomy <- read_csv2("F_Taxonomy_Microbiota.csv")
                    
                    # 2. Convertir valores CLR a numérico (excepto columnas no OTU)
                    otu_microbiota_clean <- otu_microbiota %>%
                      mutate(across(starts_with("Otu"), ~ as.numeric(.)))
                    
                    # 3. Convertir a formato largo
                    otu_long <- otu_microbiota_clean %>%
                      pivot_longer(cols = starts_with("Otu"), names_to = "OTU", values_to = "CLR_Value") %>%
                      rename(Sample = ID)
                    
                    # 4. Unir con taxonomía
                    otu_tax <- otu_long %>%
                      left_join(taxonomy %>% select(OTU, Family), by = "OTU")
                    
                    # 5. Filtrar OTUs con CLR > 0 y calcular número de OTUs únicos por familia
                    top5_familias <- otu_tax %>%
                      filter(CLR_Value > 0, !is.na(Family)) %>%
                      group_by(Family, OTU) %>%
                      summarise(has_value = TRUE, .groups = "drop") %>%
                      group_by(Family) %>%
                      summarise(n_OTUs = n(), .groups = "drop") %>%
                      arrange(desc(n_OTUs)) %>%
                      slice_head(n = 5) %>%
                      pull(Family)
                    
                    # 6. Filtrar datos CLR solo para esas familias
                    otu_top5 <- otu_tax %>%
                      filter(Family %in% top5_familias)
                    
                    # 7. Etiquetas: número de muestras por familia
                    label_data <- otu_top5 %>%
                      group_by(Family, Sample) %>%
                      summarise(Mean_CLR = mean(CLR_Value, na.rm = TRUE), .groups = "drop") %>%
                      group_by(Family) %>%
                      summarise(n = n(), .groups = "drop")
                    
                    # 8. Número de OTUs por familia (para etiqueta inferior)
                    taxa_count <- otu_tax %>%
                      filter(Family %in% top5_familias) %>%
                      distinct(Family, OTU) %>%
                      group_by(Family) %>%
                      summarise(n_OTUs = n(), .groups = "drop")
                    
                    # 9. Gráfico
                    ggplot(otu_top5, aes(x = Family, y = CLR_Value, fill = Family)) +
                      geom_violin(trim = FALSE, color = "gray", fill = "lightblue", scale = "width") +
                      geom_jitter(width = 0.2, alpha = 0.7, size = 0.1, color = "black") +
                      stat_summary(fun = mean, geom = "crossbar", color = "black", size = 0.3) +
                      scale_fill_brewer(palette = "Set3") +
                      labs(
                        title = "",
                        x = "",
                        y = "CLR"
                      ) +
                      theme_minimal(base_size = 14) +
                      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                      guides(fill = "none")             

      #    1.5 Metabolome classification an taxonomy  and plot
                    
                    setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/RAW")
                    
                    # Read the Excel file (replace with correct path)
                    MetCla <- read.xlsx("H_Metabolomics_Taxonomy_&_Classification.xlsx", sheet = 1)  # You can change the sheet number
                    
                    # Count the number of compounds per category
                    MetCla_summary <- as.data.frame(table(MetCla$classhmdb))
                    colnames(MetCla_summary) <- c("classhmdb", "count")
                    
                    # Custom color palette
                    num_colors <- length(unique(MetCla_summary$classhmdb))
                    palette_colors <- colorRampPalette(viridis::viridis(num_colors))(num_colors)
                    
                    # Plot with segmented and flipped X axis
                    ggplot(MetCla_summary, aes(y = reorder(classhmdb, -count), x = count, fill = classhmdb)) +
                      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
                      scale_fill_manual(values = palette_colors) +
                      scale_x_break(c(30, 70), scales = c(0.8, 0.2)) +  # 📌 Segment: 80% - 20%
                      coord_flip() +  # 🔄 Flip the plot
                      theme_minimal(base_size = 14) +
                      labs(
                        title = "",
                        x = "Number of Metabolites",
                        y = "Class by HMDB"
                      ) +
                      theme(
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                        axis.text.y = element_text(size = 12, hjust = 1),  # 🔥 Vertical text
                        axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor = element_blank()
                      )
                    
                    # Create frequency table for superclass
                    
                    MetCla1 <- read.xlsx("H_Metabolomics_Taxonomy_&_Classification.xlsx", sheet = 1)  # You can change the sheet number
                    MetCla2 <- read.xlsx("H_Metabolomics_Taxonomy_&_Classification.xlsx", sheet = 2)  # You can change the sheet number
                    
                    MetCla_summary_S1 <- as.data.frame(table(MetCla1$superclasshmdb))
                    MetCla_summary_S2 <- as.data.frame(table(MetCla2$Tentative_superclass_classification_no_HMDB))
                    
                    # Keep "Not available" only from S2
                    na_S2 <- MetCla_summary_S2 %>% filter(Var1 == "Not available")
                    
                    # Filter valid data (excluding "Not available" from S1)
                    MetCla_S1_filtered <- MetCla_summary_S1 %>% filter(Var1 != "Not available")
                    
                    # Filter valid data from S2 (excluding nothing)
                    MetCla_S2_filtered <- MetCla_summary_S2 %>% filter(Var1 != "Not available")
                    
                    # Combine all: S1 (without NA), S2 (without NA), and NA once (from S2)
                    MetCla_combined <- bind_rows(MetCla_S1_filtered, MetCla_S2_filtered, na_S2)
                    
                    # Group and sum frequencies
                    MetCla_summary_total <- MetCla_combined %>%
                      group_by(Var1) %>%
                      summarise(Freq = sum(Freq)) %>%
                      arrange(desc(Freq))
                    
                    # View result
                    print(MetCla_summary_total)
                    sum(MetCla_summary_total$Freq)
                    
                    # Rename columns
                    colnames(MetCla_summary_total) <- c("superclasshmdb", "count")
                    
                    # Order data by count
                    MetCla_summary_total <- MetCla_summary_total %>%
                      arrange(desc(count)) %>%
                      mutate(label_pos = count + max(count) * 0.05)  # Label position
                    
                    # Create plot
                    ggplot(MetCla_summary_total, aes(x = reorder(superclasshmdb, -count), y = count)) +
                      geom_segment(aes(xend = superclasshmdb, yend = 0), color = "black", size = 1) +  # Black lines
                      geom_point(aes(y = count), size = 9, color = "black", fill = "#ffa500", shape = 21) +  # Orange dots
                      geom_text(aes(y = count, label = count), color = "black", size = 4, fontface = "bold") +  # Count labels
                      coord_flip() +  # Flip the plot
                      theme_minimal(base_size = 14) +
                      labs(
                        title = "",
                        x = "Superclass by HMDB*",
                        y = "Number of Metabolites"
                      ) +
                      theme(
                        axis.text.y = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor = element_blank()
                      )
                    
              #1.5.1 Category
                    
                    setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/RAW")
                    
                    # Read the Excel file (replace with correct path)
                    MetCla1 <- read.xlsx("H_Metabolomics_Taxonomy_&_Classification.xlsx", sheet = 1)  # You can change the sheet number
                    MetCla1$Category
                    # Count the number of compounds per category
                    MetCla_category <- as.data.frame(table(MetCla1$Category))
                    colnames(MetCla_category) <- c("Category", "count")
                    
                    # Rename columns
                    colnames(MetCla_category) <- c("superclasshmdb", "count")
                    
                    # Order data by count
                    MetCla_category <- MetCla_category %>%
                      arrange(desc(count)) %>%
                      mutate(label_pos = count + max(count) * 0.05)  # Label position
                    
                    # Create plot
                    ggplot(MetCla_category, aes(x = reorder(superclasshmdb, -count), y = count)) +
                      geom_segment(aes(xend = superclasshmdb, yend = 0), color = "black", size = 1) +  # Black lines
                      geom_point(aes(y = count), size = 9, color = "black", fill = "#6897bb", shape = 21) +  # Orange dots
                      geom_text(aes(y = count, label = count), color = "black", size = 4, fontface = "bold") +  # Count labels
                      coord_flip() +  # Flip the plot
                      theme_minimal(base_size = 14) +
                      labs(
                        title = "",
                        x = "Category",
                        y = "Number of Metabolites"
                      ) +
                      theme(
                        axis.text.y = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor = element_blank()
                      )

                    
                    

#2.PCA (Exploratory)
    res.pca_plants_I <- PCA(B_OTU_Plants_CLR[c(-1,-2)], scale.unit = FALSE, ncp = 10, graph = TRUE)
    res.pca_metabolome <- PCA(Metabolites, scale.unit = FALSE, ncp = 10, graph = TRUE)
    res.pca_microbiota <- PCA(E_OTU_Microbiota_CLR[c(-1,-2)], scale.unit = FALSE, ncp = 10, graph = TRUE)
 

    
    #coordinates
     coords_Diet <- res.pca_plants_I$ind$coord
          dist_matrix_Diet <- vegdist(coords_Diet[, c(1, 2)], method = "euclidean")
     coords_metabolome <- res.pca_metabolome$ind$coord
          dist_matrix_metabolome <- vegdist(coords_metabolome[, c(1, 2)], method = "euclidean")
     coords_microbiota <- res.pca_microbiota$ind$coord
          dist_matrix_Microbiota <- vegdist(coords_microbiota[, c(1, 2)], method = "euclidean")
 
          
    #2.1criteria to define groups: Diet, metabolome and microbiota
       setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/PCA")
       
       #2.1.Baseline(Time1) vs Chronic (Time2) 
           #Diet
           Biplot_Baseline_Post_ByDiet = fviz_pca_biplot(res.pca_plants_I, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Group, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 2, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Time"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Group)
           ggpubr::ggpar(Biplot_Baseline_Post_ByDiet, title = "Plants intake", subtitle = "Figure.1_PCA_Time1_vs_Time2", caption = FALSE, xlab = "PC1(10.7%)", ylab = "PC2(9.7%)", legend.position = "top", legendsize = 8)
                 #permanova  
                 adonis_Biplot_Baseline_Post_ByDiet <- adonis2(dist_matrix_Diet ~ Group$Group)
                 write.table( adonis_Biplot_Baseline_Post_ByDiet,file = "Adonis_PCA/adonis_Figure.1_PCA_BaselineVS_Post_Diet.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
           #Metabolome
            Biplot_Baseline_Post_Metabolome = fviz_pca_biplot(res.pca_metabolome, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Group, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 2, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Time"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Group)
            ggpubr::ggpar(Biplot_Baseline_Post_Metabolome, title = "Metabolome profile", subtitle = "Figure.2_PCA_Time1_vs_Time2", caption = FALSE, xlab = "PC1(15.8%)", ylab = "PC2(10.5%)", legend.position = "top", legendsize = 8)
                  #permanova  
                  adonis_Biplot_Baseline_Post_Metabolome <- adonis2(dist_matrix_metabolome ~ Group$Group)
                  write.table( adonis_Biplot_Baseline_Post_Metabolome,file = "Adonis_PCA/adonis_Figure.2_PCA_BaselineVS_Post_Metabolome.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
           #Microbiota
              Biplot_Baseline_Post_Microbiota = fviz_pca_biplot(res.pca_microbiota, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Group, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 2, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Time"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Group)
              ggpubr::ggpar(Biplot_Baseline_Post_Microbiota, title = "Microbiota profile", subtitle = "Figure.3__PCA_Time1_vs_Time2", caption = FALSE, xlab = "PC1(12.2%)", ylab = "PC2(8.5%)", legend.position = "top", legendsize = 8)
                  #permanova  
                  adonis_Biplot_Baseline_Post_Microbiota <- adonis2(dist_matrix_Microbiota ~ Group$Group)
                  write.table( adonis_Biplot_Baseline_Post_Microbiota,file = "Adonis_PCA/adonis_Figure.3_PCA_BaselineVS_Post_Microbiota.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                  
            #Conclusion: The comparison of individuals between the baseline and 19 days did not show changes in diet or microbiota. It is observed that there is a difference between the two times with respect to the metabolome.
            
            
        #2.2.Diet groups 
                    
          #Diet groups identified  
          fviz_eig(res.pca_plants_I, addlabels = TRUE, ylim = c(0, 50))
                
          #Diet Groups defined by Hierarchical Clustering and Cut the Tree 
          Dendograma_Diet <- hcut(B_OTU_Plants_CLR[c(-1,-2)],hc_func = "hclust", hc_method = "ward.D", hc_metric =  "euclidean", k = 3, stand = TRUE)
          fviz_dend(Dendograma_Diet, rect = TRUE, cex = 1, labelsize = 10, k_colors = c("#2E9FDF", "#E7B800", "#FC4E07"), main = "Figure.4_Dendogram Diet - hclust_ward.D_Euclidean")
        
                      
                    #This classification was used to establishing two types of diet -> Group$Two_Diet.Type_Dendrogram 
                    Group$Two_Diet.Type_Dendrogram 
                    
                          
                    #Two diet by distance matrix
                    Biplot_Diet_Dendogram = fviz_pca_biplot(res.pca_plants_I, select.var = list(cos2 = 0.25), geom.var = c("point", "text"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Two_Diet.Type_Dendrogram, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 3, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Diet"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Two_Metabolome.type)
                    ggpubr::ggpar(Biplot_Diet_Dendogram, title = "By Diet", subtitle = "Figure.5_PCA_Dendogram_Diet", caption = FALSE, xlab = "PC1(10.7%)", ylab = "PC2(9.7%)", legend.position = "top", legendsize = 8)
                    
                        # Description of dimension to define plants associated to the three components
                        res.desc_Plants_I <- dimdesc(res.pca_plants_I, axes = c(1,2,3), proba = 0.05) 
                        res.desc_Plants_I$Dim.1
                        res.desc_Plants_I$Dim.2
                        res.desc_Plants_I$Dim.3
                        
                        #Export excel
                        write.table(res.desc_Plants_I$Dim.1,file = "Correlations_PCA/Figure_5_CorrelationPC1_Plants.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                        write.table(res.desc_Plants_I$Dim.2,file = "Correlations_PCA/Figure_5_CorrelationPC2_Plants.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                    
                      #permanova  
                       adonis_result_Diet_dendogram <- adonis2(dist_matrix_Diet ~ Group$Two_Diet.Type_Dendrogram)
                       write.table(adonis_result_Diet_dendogram,file = "Adonis_PCA/adonis_Figure.5_PCA_Dendogram_Diet.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                       
                            #By batch
                                  Biplot_batch_trnl = fviz_pca_biplot(res.pca_plants_I, geom.var = FALSE, axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Batch_tRNL, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 3, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Group"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "euclid", addEllipses = TRUE, group = Group$Batch_tRNL)
                                  ggpubr::ggpar(Biplot_batch_trnl, title = "By batch trnl", subtitle = "PC1vsPC2", caption = FALSE, xlab = "PC1(15.9%)", ylab = "PC2(10.5%)", legend.position = "top", legendsize = 8)
                                  #permanova  
                                  adonis_result_metabolome_batch_trnl <- adonis2(dist_matrix_Diet ~ Group$Batch_tRNL)
                                  write.table( adonis_result_metabolome_batch_trnl,file = "Adonis_PCA/adonis_result_metabolome_batch_trnl.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                  
         #2.3.Metabolome groups 
           #Full database uchuva proyet. TENTATIVE 
              eig.val_pca_METABOLOME<- get_eigenvalue(res.pca_metabolome) #Use this value for visualization of PCA biplot final
              fviz_eig(res.pca_metabolome, addlabels = TRUE, ylim = c(0, 50))
        
              #Two metabolomes
              Biplot_Metabolome_PC1_PC2 = fviz_pca_biplot(res.pca_metabolome, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Two_Metabolome.type, palette = c("darkblue", "#fc4e07"), pointshape = 21, pointsize = 2, labelsize = 3, center = c(0, 1), fill.var = "gray", col.var = "black", legend.title = list(fill = "Metabolome profile"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Two_Metabolome.type)
              ggpubr::ggpar(Biplot_Metabolome_PC1_PC2, title = "By Metabolome", subtitle = "Figure.6_PCA_Metabolome_PC1_PC2", caption = FALSE, xlab = "PC1(15.8%)", ylab = "PC2(10.5%)", legend.position = "top", legendsize = 8)
              c("darkblue", "#fc4e07")

                    # Description of dimension 
                    res.desc_Metabolome <- dimdesc(res.pca_metabolome, axes = c(1,2,3), proba = 0.05) #res.pca se optiene del codigo anterior en PCA
                    res.desc_Metabolome$Dim.1
                    res.desc_Metabolome$Dim.2
                    res.desc_Metabolome$Dim.3
                    
                    #Export excel
                    write.table(res.desc_Metabolome$Dim.1,file = "Correlations_PCA/Figure_6_CorrelationPC1_Metabolome.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                    write.table(res.desc_Metabolome$Dim.2,file = "Correlations_PCA/Figure_6_CorrelationPC2_Metabolome.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 

              #permanova  
              adonis_result_metabolome <- adonis2(dist_matrix_metabolome ~ Group$Two_Metabolome.type)
              write.table(adonis_result_metabolome,file = "Adonis_PCA/adonis_Figure.6_PCA_Metabolome_PC1_PC2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
              
              

          #2.4.MICROBIOTA groups 
              #Full database uchuva proyet. 
              eig.val_pca_Microbiota<- get_eigenvalue(res.pca_microbiota) #Use this value for visualization of PCA biplot final
              fviz_eig(res.pca_microbiota, addlabels = TRUE, ylim = c(0, 50))
        
              #Three microbiota profile 
              Biplot_Microbiota__PC1_PC2 = fviz_pca_biplot(res.pca_microbiota, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Microbiota.profile, palette = c("#00AFBB","#ecd817", "#b8360a"), pointshape = 21, pointsize = 3, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Microbiota"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Microbiota.profile)
              ggpubr::ggpar(Biplot_Microbiota__PC1_PC2, title = "By Microbiota", subtitle = "Figure.7_PCA_Microbiota__PC1_PC2", caption = FALSE, xlab = "PC1(12.2%)", ylab = "PC2(8.5%)", legend.position = "top", legendsize = 8)
              
              # Description of dimension 
              res.desc_Microbiota <- dimdesc(res.pca_microbiota, axes = c(1,2,3), proba = 0.05) #res.pca se optiene del codigo anterior en PCA
              res.desc_Microbiota$Dim.1
              res.desc_Microbiota$Dim.2
              res.desc_Microbiota$Dim.3
              
                      #permanova  
                      adonis_result_microbiota <- adonis2(dist_matrix_Microbiota ~ Group$Microbiota.profile)
                      write.table(adonis_result_metabolome,file = "Adonis_PCA/adonis_Figure.7_PCA_Microbiota__PC1_PC2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      
        #2.1 Estratification by subject
                      
                      # --- INICIO DEL CÓDIGO ACTUALIZADO ---
                      
                      # 0. Preparar la estratificación por sujeto (CORREGIDO)
                      
                      # Extrae solo el número inicial (ej. convierte "2T1" y "2T2" en "2")
                      Subject_ID_Short <- sub("T.*", "", rownames(Group))
                      
                      # Verificar que ahora tienes solo 9 IDs únicos (o 18 si son 18 sujetos en total)
                      print(unique(Subject_ID_Short)) 
                      # Deberías ver: [1] "2" "3" "4" "5" "7" "8" "9" "10" "11" "12" ... etc.
                      
                      # Definir el control de permutaciones bloqueado por el nuevo ID de sujeto
                      library(permute)
                      ctrl <- how(plots = Plots(strata = Subject_ID_Short), nperm = 999)
                      
                      # Alinear la matriz de distancias y la variable Group antes del PERMANOVA
                      # Esto es CRUCIAL para evitar el error 'subíndice fuera de los límites'
                      dist_matrix_Diet_aligned <- as.matrix(dist_matrix_Diet)[rownames(Group), rownames(Group)]                      
                      
                      
                      #2.PCA (Exploratory)
                      res.pca_plants_I <- PCA(B_OTU_Plants_CLR[c(-1,-2)], scale.unit = FALSE, ncp = 10, graph = TRUE)
                      res.pca_metabolome <- PCA(Metabolites, scale.unit = FALSE, ncp = 10, graph = TRUE)
                      res.pca_microbiota <- PCA(E_OTU_Microbiota_CLR[c(-1,-2)], scale.unit = FALSE, ncp = 10, graph = TRUE)
                      
                      #coordinates
                      coords_Diet <- res.pca_plants_I$ind$coord
                      dist_matrix_Diet <- vegdist(coords_Diet[, c(1, 2)], method = "euclidean")
                      coords_metabolome <- res.pca_metabolome$ind$coord
                      dist_matrix_metabolome <- vegdist(coords_metabolome[, c(1, 2)], method = "euclidean")
                      coords_microbiota <- res.pca_microbiota$ind$coord
                      dist_matrix_Microbiota <- vegdist(coords_microbiota[, c(1, 2)], method = "euclidean")
                      
                      
                      #2.1criteria to define groups: Diet, metabolome and microbiota
                      setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/PCA")
                      
                      #2.1.Baseline(Time1) vs Chronic (Time2) 
                      #Diet
                      Biplot_Baseline_Post_ByDiet = fviz_pca_biplot(res.pca_plants_I, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Group, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 2, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Time"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Group)
                      ggpubr::ggpar(Biplot_Baseline_Post_ByDiet, title = "Plants intake", subtitle = "Figure.1_PCA_Time1_vs_Time2", caption = FALSE, xlab = "PC1(10.7%)", ylab = "PC2(9.7%)", legend.position = "top", legendsize = 8)
                      #permanova  
                      # MODIFICADO: Se añade 'permutations = ctrl'
                      adonis_Biplot_Baseline_Post_ByDiet <- adonis2(dist_matrix_Diet ~ Group$Group, permutations = ctrl)
                      write.table( adonis_Biplot_Baseline_Post_ByDiet,file = "Adonis_PCA/adonis_Figure.1_PCA_BaselineVS_Post_Diet.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      #Metabolome
                      Biplot_Baseline_Post_Metabolome = fviz_pca_biplot(res.pca_metabolome, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Group, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 2, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Time"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Group)
                      ggpubr::ggpar(Biplot_Baseline_Post_Metabolome, title = "Metabolome profile", subtitle = "Figure.2_PCA_Time1_vs_Time2", caption = FALSE, xlab = "PC1(15.8%)", ylab = "PC2(10.5%)", legend.position = "top", legendsize = 8)
                      #permanova  
                      # MODIFICADO: Se añade 'permutations = ctrl'
                      adonis_Biplot_Baseline_Post_Metabolome <- adonis2(dist_matrix_metabolome ~ Group$Group, permutations = ctrl)
                      write.table( adonis_Biplot_Baseline_Post_Metabolome,file = "Adonis_PCA/adonis_Figure.2_PCA_BaselineVS_Post_Metabolome.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      #Microbiota
                      Biplot_Baseline_Post_Microbiota = fviz_pca_biplot(res.pca_microbiota, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Group, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 2, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Time"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Group)
                      ggpubr::ggpar(Biplot_Baseline_Post_Microbiota, title = "Microbiota profile", subtitle = "Figure.3__PCA_Time1_vs_Time2", caption = FALSE, xlab = "PC1(12.2%)", ylab = "PC2(8.5%)", legend.position = "top", legendsize = 8)
                      #permanova  
                      # MODIFICADO: Se añade 'permutations = ctrl'
                      adonis_Biplot_Baseline_Post_Microbiota <- adonis2(dist_matrix_Microbiota ~ Group$Group, permutations = ctrl)
                      write.table( adonis_Biplot_Baseline_Post_Microbiota,file = "Adonis_PCA/adonis_Figure.3_PCA_BaselineVS_Post_Microbiota.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      
                      #Conclusion: The comparison of individuals between the baseline and 19 days did not show changes in diet or microbiota. It is observed that there is a difference between the two times with respect to the metabolome.
                      
                      
                      #2.2.Diet groups 
                      
                      #Diet groups identified  
                      fviz_eig(res.pca_plants_I, addlabels = TRUE, ylim = c(0, 50))
                      
                      #Diet Groups defined by Hierarchical Clustering and Cut the Tree 
                      Dendograma_Diet <- hcut(B_OTU_Plants_CLR[c(-1,-2)],hc_func = "hclust", hc_method = "ward.D", hc_metric =  "euclidean", k = 3, stand = TRUE)
                      fviz_dend(Dendograma_Diet, rect = TRUE, cex = 1, labelsize = 10, k_colors = c("#2E9FDF", "#E7B800", "#FC4E07"), main = "Figure.4_Dendogram Diet - hclust_ward.D_Euclidean")
                      
                      
                      #This classification was used to establishing two types of diet -> Group$Two_Diet.Type_Dendrogram 
                      Group$Two_Diet.Type_Dendrogram 
                      
                      
                      #Two diet by distance matrix
                      Biplot_Diet_Dendogram = fviz_pca_biplot(res.pca_plants_I, select.var = list(cos2 = 0.25), geom.var = c("point", "text"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Two_Diet.Type_Dendrogram, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 3, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Diet"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Two_Metabolome.type)
                      ggpubr::ggpar(Biplot_Diet_Dendogram, title = "By Diet", subtitle = "Figure.5_PCA_Dendogram_Diet", caption = FALSE, xlab = "PC1(10.7%)", ylab = "PC2(9.7%)", legend.position = "top", legendsize = 8)
                      
                      # Description of dimension to define plants associated to the three components
                      res.desc_Plants_I <- dimdesc(res.pca_plants_I, axes = c(1,2,3), proba = 0.05) 
                      res.desc_Plants_I$Dim.1
                      res.desc_Plants_I$Dim.2
                      res.desc_Plants_I$Dim.3
                      
                      #Export excel
                      write.table(res.desc_Plants_I$Dim.1,file = "Correlations_PCA/Figure_5_CorrelationPC1_Plants.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      write.table(res.desc_Plants_I$Dim.2,file = "Correlations_PCA/Figure_5_CorrelationPC2_Plants.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      
                      #permanova  
                      # MODIFICADO: Se añade 'permutations = ctrl'
                      adonis_result_Diet_dendogram <- adonis2(dist_matrix_Diet ~ Group$Two_Diet.Type_Dendrogram, permutations = ctrl)
                      write.table(adonis_result_Diet_dendogram,file = "Adonis_PCA/adonis_Figure.5_PCA_Dendogram_Diet.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      
                      #By batch
                      Biplot_batch_trnl = fviz_pca_biplot(res.pca_plants_I, geom.var = FALSE, axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Batch_tRNL, palette = c("#00AFBB", "#b8360a"), pointshape = 21, pointsize = 3, labelsize = 3, center = c(0, 1), fill.var = "black", col.var = "black", legend.title = list(fill = "Group"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "euclid", addEllipses = TRUE, group = Group$Batch_tRNL)
                      ggpubr::ggpar(Biplot_batch_trnl, title = "By batch trnl", subtitle = "PC1vsPC2", caption = FALSE, xlab = "PC1(15.9%)", ylab = "PC2(10.5%)", legend.position = "top", legendsize = 8)
                      #permanova
                      
                      
                      #Metaboloma
                      

                      #2.3.Metabolome groups 
                      #Full database uchuva proyet. TENTATIVE 
                      eig.val_pca_METABOLOME<- get_eigenvalue(res.pca_metabolome) #Use this value for visualization of PCA biplot final
                      fviz_eig(res.pca_metabolome, addlabels = TRUE, ylim = c(0, 50))
                      
                      #Two metabolomes
                      Biplot_Metabolome_PC1_PC2 = fviz_pca_biplot(res.pca_metabolome, geom.var = c("point"), axes = c(1,2), geom.ind = c("text","point"), fill.ind = Group$Two_Metabolome.type, palette = c("darkblue", "#fc4e07"), pointshape = 21, pointsize = 2, labelsize = 3, center = c(0, 1), fill.var = "gray", col.var = "black", legend.title = list(fill = "Metabolome profile"), repel = TRUE, mean.point = TRUE, mean.pointsize = 20, ellipse.type = "confidence", addEllipses = TRUE, group = Group$Two_Metabolome.type)
                      ggpubr::ggpar(Biplot_Metabolome_PC1_PC2, title = "By Metabolome", subtitle = "Figure.6_PCA_Metabolome_PC1_PC2", caption = FALSE, xlab = "PC1(15.8%)", ylab = "PC2(10.5%)", legend.position = "top", legendsize = 8)
                      c("darkblue", "#fc4e07")
                      
                      # Description of dimension 
                      res.desc_Metabolome <- dimdesc(res.pca_metabolome, axes = c(1,2,3), proba = 0.05) #res.pca se optiene del codigo anterior en PCA
                      res.desc_Metabolome$Dim.1
                      res.desc_Metabolome$Dim.2
                      res.desc_Metabolome$Dim.3
                      
                      #Export excel
                      write.table(res.desc_Metabolome$Dim.1,file = "Correlations_PCA/Figure_6_CorrelationPC1_Metabolome.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      write.table(res.desc_Metabolome$Dim.2,file = "Correlations_PCA/Figure_6_CorrelationPC2_Metabolome.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      
                      # Alinear la matriz de distancias y la variable Group antes del PERMANOVA (CRUCIAL)
                      dist_matrix_metabolome_aligned <- as.matrix(dist_matrix_metabolome)[rownames(Group), rownames(Group)]
                      
                      # permanova (MODIFICADO: Se añade 'permutations = ctrl' y se usa la matriz alineada)
                      adonis_result_metabolome <- adonis2(dist_matrix_metabolome_aligned ~ Group$Two_Metabolome.type, permutations = ctrl)
                      
                      write.table(adonis_result_metabolome,file = "Adonis_PCA/adonis_Figure.6_PCA_Metabolome_PC1_PC2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE)
                      
                      # Imprimir resultados para verificar
                      print(adonis_result_metabolome)
                      
                      
                           
#3.PLSDA-MULTIBLOCK (MixOmics-DIABLO MODEL)
              
              
          #3.1 Metadata for paper metadiet
                M <- list(Plants = Plants[[1]], 
                          Metabolites = Metabolites[[1]], 
                          Microbiota = Microbiota[[1]])
                
                            library(caret)
                            # Identifying plants with zero variance
                            zero_var_features <- nearZeroVar(M$Plants, saveMetrics = TRUE)
                            # View plants with zero variance
                            print(zero_var_features$nzv)
                            # Filter plants with zero variance and set it to M.
                            M$Plants <- M$Plants[, !zero_var_features$nzv]
          
                     
                            otu_clr_summary       
                 #Nutritional category
                            
                            M2 <- list(Plants_N = B_OTU_Plants_CLR_nutritional_class, 
                                      Metabolites = Metabolites[[1]], 
                                      Microbiota = Microbiota[[1]])
                            
                            library(caret)
                            # Identifying plants with zero variance
                            zero_var_features_N <- nearZeroVar(M2$Plants_N, saveMetrics = TRUE)
                            # View plants with zero variance
                            print(zero_var_features_N$nzv)
                            # Filter plants with zero variance and set it to M.
                            M2$Plants_N <- M2$Plants_N[, !zero_var_features_N$nzv]
                            
          
          #3.2.Correlation between matrices using PLS
                setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output")
          
              #Correlation and P-Value between matrices using PLS for Metadiet paper
                  
                #Plants vs Metabolites
                  res1.pls.METADIEt <- pls(Metabolites[[1]], Plants[[1]], ncomp = 1) 
                  cor(res1.pls.METADIEt$variates$X, res1.pls.METADIEt$variates$Y)
                  cor.test(res1.pls.METADIEt$variates$X, res1.pls.METADIEt$variates$Y)[["p.value"]]
                      write.table(cor(res1.pls.METADIEt$variates$X, res1.pls.METADIEt$variates$Y),file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_R2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      write.table(cor.test(res1.pls.METADIEt$variates$X, res1.pls.METADIEt$variates$Y)[["p.value"]],file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_Pvalue.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      
                            res1.pls.METADIEt2 <- pls(Metabolites[[1]], Plants[[1]], ncomp = 2) 
                            cor(res1.pls.METADIEt2$variates$X, res1.pls.METADIEt2$variates$Y)
                            cor.test(res1.pls.METADIEt2$variates$X, res1.pls.METADIEt2$variates$Y)[["p.value"]]
                            write.table(cor(res1.pls.METADIEt2$variates$X, res1.pls.METADIEt2$variates$Y),file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_R2_pc2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                            write.table(cor.test(res1.pls.METADIEt2$variates$X, res1.pls.METADIEt2$variates$Y)[["p.value"]],file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_Pvalue_pc2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 

                                  #PLOT ARROW Figure.0a  CORRELATION BY PLS Metabolites VS plants
                                  res1.pls.ARROW <- pls(Metabolites[[1]], Plants[[1]])
                                    plotArrow(res1.pls.ARROW, ind.names = FALSE, 
                                              group = Group$Group, 
                                              col.per.group = c("darkblue", "#fc4e07", "forestgreen"),
                                              comp = c(1,2), 
                                              ind.names.size = 5, legend = TRUE, arrow.size = 0.3,arrow.length = 0.3, pch.size = 2, arrow.lwd = 4,
                                              title = 'Metabolites vs Plants') 
                                              par(new = TRUE) # Permite superponer gráficos en la misma ventana
                                              plot(0, 0, type = "n", xlim = par("usr")[1:2], ylim = par("usr")[3:4], xlab = "", ylab = "", axes = FALSE)
                                              grid(col = "gray", lty = "dotted", lwd = 3)
                                              text(x = 1, y = max(par("usr")[3:4]) * 0.9, 
                                                   labels = "Metabolome vs Plants", 
                                                   col = "black", cex = 1.5, font = 2)
                                              

                                              
                                              
                                              
                                             #nutritional
                                              #Plants vs Metabolites
                                              res1.pls.METADIEt_N <- pls(Metabolites[[1]], B_OTU_Plants_CLR_nutritional_class, ncomp = 1) 
                                              cor(res1.pls.METADIEt_N$variates$X, res1.pls.METADIEt_N$variates$Y)
                                              cor.test(res1.pls.METADIEt_N$variates$X, res1.pls.METADIEt_N$variates$Y)[["p.value"]]
                                              write.table(cor(res1.pls.METADIEt_N$variates$X, res1.pls.METADIEt_N$variates$Y),file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_R2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                              write.table(cor.test(res1.pls.METADIEt_N$variates$X, res1.pls.METADIEt_N$variates$Y)[["p.value"]],file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_Pvalue.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                              
                                              res1.pls.METADIEt2 <- pls(Metabolites[[1]], Plants[[1]], ncomp = 2) 
                                              cor(res1.pls.METADIEt2$variates$X, res1.pls.METADIEt2$variates$Y)
                                              cor.test(res1.pls.METADIEt2$variates$X, res1.pls.METADIEt2$variates$Y)[["p.value"]]
                                              write.table(cor(res1.pls.METADIEt2$variates$X, res1.pls.METADIEt2$variates$Y),file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_R2_pc2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                              write.table(cor.test(res1.pls.METADIEt2$variates$X, res1.pls.METADIEt2$variates$Y)[["p.value"]],file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_Pvalue_pc2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                              
                                              #PLOT ARROW Figure.0a  CORRELATION BY PLS Metabolites VS plants
                                              res1.pls.ARROW_N <- pls(Metabolites[[1]], B_OTU_Plants_CLR_nutritional_class)
                                              plotArrow(res1.pls.ARROW_N, ind.names = FALSE, 
                                                        group = Group$Group, 
                                                        col.per.group = c("darkblue", "#fc4e07", "forestgreen"),
                                                        comp = c(1,2), 
                                                        ind.names.size = 5, legend = TRUE, arrow.size = 0.3,arrow.length = 0.3, pch.size = 2, arrow.lwd = 4,
                                                        title = 'Metabolites vs Plants') 
                                              par(new = TRUE) # Permite superponer gráficos en la misma ventana
                                              plot(0, 0, type = "n", xlim = par("usr")[1:2], ylim = par("usr")[3:4], xlab = "", ylab = "", axes = FALSE)
                                              grid(col = "gray", lty = "dotted", lwd = 3)
                                              text(x = 1, y = max(par("usr")[3:4]) * 0.9, 
                                                   labels = "Metabolome vs Plants", 
                                                   col = "black", cex = 1.5, font = 2)
                                              
                                              
                                              
                                              
                                              
                #Metabolites vs Microbiota
                  res2.pls.METADIEt <- pls(Metabolites[[1]], Microbiota[[1]], ncomp = 1)
                  cor(res2.pls.METADIEt$variates$X, res2.pls.METADIEt$variates$Y)
                  cor.test(res2.pls.METADIEt$variates$X, res2.pls.METADIEt$variates$Y)[["p.value"]]
                      write.table(cor(res2.pls.METADIEt$variates$X, res2.pls.METADIEt$variates$Y),file = "PLS_correlation_between_matrices/PLS_Correlation_Microbiota_vs_Metabolites_R2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      write.table(cor.test(res2.pls.METADIEt$variates$X, res2.pls.METADIEt$variates$Y)[["p.value"]],file = "PLS_correlation_between_matrices/PLS_Correlation_Microbiota_vs_Metabolites_Pvalue.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 

                                  res2.pls.METADIEt2 <- pls(Metabolites[[1]], Microbiota[[1]], ncomp = 2 , scale = TRUE) 
                                  cor(res2.pls.METADIEt2$variates$X, res2.pls.METADIEt2$variates$Y)
                                  cor.test(res2.pls.METADIEt2$variates$X, res2.pls.METADIEt2$variates$Y)[["p.value"]]
                                  write.table(cor(res2.pls.METADIEt2$variates$X, res2.pls.METADIEt2$variates$Y),file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_R2_pc2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                  write.table(cor.test(res2.pls.METADIEt2$variates$X, res2.pls.METADIEt2$variates$Y)[["p.value"]],file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_Pvalue_pc2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                  
                                  #PLOT ARROW Figure.0b  CORRELATION BY PLS Metabolites VS Microbiota_comp_1-2
                                  res2.pls.ARROW <- pls(Metabolites[[1]], Microbiota[[1]])
                                  plotArrow(res2.pls.ARROW, ind.names = FALSE, 
                                            group = Group$Group, 
                                            col.per.group = c("darkblue", "#fc4e07", "forestgreen"),
                                            comp = c(1,2), 
                                            ind.names.size = 5, legend = TRUE, arrow.size = 0.3,arrow.length = 0.3, pch.size = 2, arrow.lwd = 4,
                                            title = 'Metabolome vs Microbiota') 
                                  par(new = TRUE) # Permite superponer gráficos en la misma ventana
                                  plot(0, 0, type = "n", xlim = par("usr")[1:2], ylim = par("usr")[3:4], xlab = "", ylab = "", axes = FALSE)
                                  grid(col = "gray80", lty = "dotted", lwd = 3)
                                  text(x = 1, y = max(par("usr")[3:4]) * 0.9, 
                                       labels = "Metabolome vs Microbiota", 
                                       col = "black", cex = 1.5, font = 2)
          
                #Plants vs Microbiota
                  res3.pls.METADIEt <- pls(Plants[[1]], Microbiota[[1]], ncomp = 1)
                  cor(res3.pls.METADIEt$variates$X, res3.pls.METADIEt$variates$Y)
                  cor.test(res3.pls.METADIEt$variates$X, res3.pls.METADIEt$variates$Y)[["p.value"]]
                      write.table(cor(res3.pls.METADIEt$variates$X, res3.pls.METADIEt$variates$Y),file = "PLS_correlation_between_matrices/PLS_Correlation_Microbiota_vs_Plants_R2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      write.table(cor.test(res3.pls.METADIEt$variates$X, res3.pls.METADIEt$variates$Y)[["p.value"]],file = "PLS_correlation_between_matrices/PLS_Correlation_Microbiota_vs_Plants_Pvalue.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                  
                  
                      res3.pls.METADIEt2 <- pls(Plants[[1]], Microbiota[[1]], ncomp = 2 , scale = TRUE) 
                      cor(res3.pls.METADIEt2$variates$X, res3.pls.METADIEt2$variates$Y)
                      cor.test(res3.pls.METADIEt2$variates$X, res3.pls.METADIEt2$variates$Y)[["p.value"]]
                      write.table(cor(res3.pls.METADIEt2$variates$X, res3.pls.METADIEt2$variates$Y),file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_R2_pc2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      write.table(cor.test(res3.pls.METADIEt2$variates$X, res3.pls.METADIEt2$variates$Y)[["p.value"]],file = "PLS_correlation_between_matrices/PLS_Correlation_Plants_vs_Metabolites_Pvalue_pc2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      
                                      #PLOT ARROW
                                      res3.pls.arrow <- pls(Plants[[1]], Microbiota[[1]])
                                      plotArrow(res3.pls.arrow, ind.names = FALSE, 
                                                group = Group$Group, 
                                                col.per.group = c("darkblue", "#fc4e07", "forestgreen"),
                                                comp = c(1,2), 
                                                ind.names.size = 5, legend = TRUE, arrow.size = 0.3,arrow.length = 0.3, pch.size = 2, arrow.lwd = 4,
                                                title = 'Plants vs Microbiota') 
                                      par(new = TRUE) # Permite superponer gráficos en la misma ventana
                                      plot(0, 0, type = "n", xlim = par("usr")[1:2], ylim = par("usr")[3:4], xlab = "", ylab = "", axes = FALSE)
                                      grid(col = "gray80", lty = "dotted", lwd = 3)
                                      text(x = 1, y = max(par("usr")[3:4]) * 0.9, 
                                           labels = "Plants vs Microbiota", 
                                           col = "black", cex = 1.2, font = 2)
                                      
                      
          
                                      
            #3.2.Parameter choice
                                      
            #design for Metadiet (It is the closest to 1 since we want to demonstrate if the matrices correlate.)
                      
                                      #METABOLOME                
                    design1 <- matrix(0.9, ncol = length(M), nrow = length(M), 
                                dimnames = list(names(M), names(M))) #Inicialmente realice el diseño con 1 dada la alta correalción
                                diag(design1) <- 0
                                design1 
                          


          #3.4.MixOmics_DIABLO
             
            #Step 1 Generation of plsda blockS
                  # Blocks to metadiet:The first 5 components are chosen for the analysis since this is where the greatest amount of explained variability is found                diablo.metadiet_DietaPC1_PC2 <- block.plsda(M, Group$Two_Diet.Type_Dendrogram, ncomp = 10, design = design1)
                    diablo.metadiet_DietDendograma<- block.plsda(M, Group$Two_Diet.Type_Dendrogram, ncomp = 5, design = design1)
                    diablo.metadiet_Metabolome.typePC1_PC2 <- block.plsda(M, Group$Two_Metabolome.type, ncomp = 5, design = design1)
                    diablo.metadiet_MicrobiotaPC1_PC2 <- block.plsda(M, Group$Microbiota.profile, ncomp = 5, design = design1)
          
            #Step 2 Seed and crossvalidation
              set.seed(123) # For reproducibility, remove for your analyses
          
                    #Blocks of metadiet and crossvalidation 
                    perf.diablo.Diet_Dendogram = perf(diablo.metadiet_DietDendograma, validation = 'Mfold', folds = 10, nrepeat = 10)
                    perf.diablo.Metabolome.typePC1_PC2 = perf(diablo.metadiet_Metabolome.typePC1_PC2, validation = 'Mfold', folds = 10, nrepeat = 10)
                    perf.diablo.MicrobiotaPC1_PC2 = perf(diablo.metadiet_MicrobiotaPC1_PC2, validation = 'Mfold', folds = 10, nrepeat = 10)
              
          
            #Step 3. Number of components to Metadiet after crossvaldation
                      #perf.diablo.tcga$error.rate  # Lists the different types of error rates
                      # Plot of the error rates based on weighted vote
          
                      #Diet
                      plot(perf.diablo.Diet_Dendogram)
                      perf.diablo.Diet_Dendogram$choice.ncomp$WeightedVote #optimal number of components according to the prediction distance and type of error rate (overall or balanced)
                      ncomp_Diet <- perf.diablo.Diet_Dendogram$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] # final ncomp value
                           
                      #Metabolome
                      plot(perf.diablo.Metabolome.typePC1_PC2)
                      perf.diablo.Metabolome.typePC1_PC2$choice.ncomp$WeightedVote #optimal number of components according to the prediction distance and type of error rate (overall or balanced)
                      ncomp_Metabolome <- perf.diablo.Metabolome.typePC1_PC2$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] # final ncomp value
          
                      #Microbiota
                      plot(perf.diablo.MicrobiotaPC1_PC2)
                      perf.diablo.MicrobiotaPC1_PC2$choice.ncomp$WeightedVote #optimal number of components according to the prediction distance and type of error rate (overall or balanced)
                      ncomp_Microbiota <- perf.diablo.MicrobiotaPC1_PC2$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] # final ncomp value
          
            # Step 4.Number of variables to select and final model
            
                      set.seed(123) # for reproducibility
                            #Paper Metadiet
                      
                                #Blocks to cross validation
                                # test.keepX for all models
                                    test.keepX <- list(Plants =c(5:15, seq(20, 50, 5)),
                                                      Metabolites = c(5:20, seq(25, 100, 25)),
                                                      Microbiota = c(seq(5, 50, 5), seq(60, 200, 20)))
                                        
                                    
                                #Final Model BY Diet
                                        tune.diablo.Diet_PC1_PC2 <- tune.block.splsda(M, Group$Two_Diet.Type_Dendrogram, ncomp_Diet, 
                                                                    test.keepX = test.keepX, design = design1,
                                                                    validation = 'Mfold', folds = 10, nrepeat = 10, 
                                                                    BPPARAM = BiocParallel::SnowParam(workers = 2),
                                                                    dist = "centroids.dist")
                                        list.keepM_Diet_PC1_PC2 <- tune.diablo.Diet_PC1_PC2$choice.keepX
                                        #Final model
                                        diablo.metadiet_Diet_Dendogram_Final <- block.splsda(M, Group$Two_Diet.Type_Dendrogram, ncomp_Diet, 
                                                                              keepX = list.keepM_Diet_PC1_PC2, design = design1)
                                    
                                #Model final BY Metabolome
                                        tune.diablo.Metabolome_PC1_PC2_F <- tune.block.splsda(M, Group$Two_Metabolome.type, ncomp_Metabolome, 
                                                                          test.keepX = test.keepX, design = design1,
                                                                          validation = 'Mfold', folds = 10, nrepeat = 10, 
                                                                          BPPARAM = BiocParallel::SnowParam(workers = 2),
                                                                          dist = "centroids.dist")
                                        list.keepM_Metabolome_PC1_PC2 <- tune.diablo.Metabolome_PC1_PC2_F$choice.keepX
                                                
                                        #Final model
                                        diablo.metadiet_Metabolome_PC1_PC2_Final <- block.splsda(M, Group$Two_Metabolome.type, ncomp_Metabolome, 
                                                                                    keepX = list.keepM_Metabolome_PC1_PC2, design = design1)
          
                                                 #Permutation
                                              
                                                          perf.diablo.Metabolome_PC1_PC2 <- perf(
                                                            diablo.metadiet_Metabolome_PC1_PC2_Final, # Modelo final
                                                            validation = 'Mfold',                    # Validación cruzada M-fold
                                                            folds = 10,                              # Número de folds
                                                            nrepeat = 10,                            # Número de repeticiones
                                                            dist = "centroids.dist",                 # Distancia usada para clasificación
                                                            progressBar = TRUE,                      # Mostrar barra de progreso
                                                            perm = 100                               # Número de permutaciones (ajustar según recursos)
                                                          )
                                                          # Imprimir resultados
                                                          print(perf.diablo.Metabolome_PC1_PC2)
                                                          # Métricas de clasificación
                                                          perf.diablo.Metabolome_PC1_PC2$error.rate      # Tasas de error por componentes
                                                          perf.diablo.Metabolome_PC1_PC2$per.class.error # Error por clase
                                                          # Resultados de la permutación
                                                          perf.diablo.Metabolome_PC1_PC2$perm.dist 
                                                          # Comparar distancias permutadas con las observadas
                                                          plot(perf.diablo.Metabolome_PC1_PC2, col = "blue", main = "Performance Metrics with Permutation Test")
                                              
                                              
                                #FinalModel  BY Microbiota
                                       tune.diablo.Microbiota_PC1_PC2_F <- tune.block.splsda(M, Group$Microbiota.profile, ncomp_Microbiota, 
                                                                         test.keepX = test.keepX, design = design1,
                                                                         validation = 'Mfold', folds = 10, nrepeat = 10, 
                                                                         BPPARAM = BiocParallel::SnowParam(workers = 2),
                                                                         dist = "centroids.dist")
                                      list.keepM_Microbiota_PC1_PC2 <- tune.diablo.Microbiota_PC1_PC2_F$choice.keepX
                                              
                                      #Final model
                                      diablo.metadiet_Microbiota_PC1_PC2_Final <- block.splsda(M, Group$Microbiota.profile, ncomp_Microbiota, 
                                                                                  keepX = list.keepM_Microbiota_PC1_PC2, design = design1)
                            
                  
                  
          #3.5.Correlation in PLS1 and PLS2
          
                  setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Correlations_between_Matrices_By_Groups")
                                      
                  #Metadiet
                  plotDiablo(diablo.metadiet_Diet_Dendogram_Final, ncomp = 1) #"Figure.8a_DIABLO_Correlation_between_Matrices_By_Diet_Component_1"
                  plotDiablo(diablo.metadiet_Diet_Dendogram_Final, ncomp = 2) #Best model "Figure.8b_DIABLO_Correlation_between_Matrices_By_Diet_Component_2"
          
                              plotDiablo(
                                diablo.metadiet_Diet_Dendogram_Final,
                                ncomp = 1,
                                col = c("darkblue", "#fc4e07"),  # Cambia según tus grupos
                                cex = 5                                # Tamaño de los puntos
                            )
                  
                  
                  
                  plotDiablo(diablo.metadiet_Metabolome_PC1_PC2_Final, ncomp = 1) #Best "Figure.9a_DIABLO_Correlation_between_Matrices_By_Metabolome_Component_1"
                  plotDiablo(diablo.metadiet_Metabolome_PC1_PC2_Final, ncomp = 2) #"Figure.9b_DIABLO_Correlation_between_Matrices_By_Metabolome_Component_2"
          
                            plotDiablo(
                              diablo.metadiet_Metabolome_PC1_PC2_Final,
                              ncomp = 1,
                              col = c("darkblue", "#fc4e07"),  # Cambia según tus grupos
                              cex = 5                                # Tamaño de los puntos
                            )
                          
                  plotDiablo(diablo.metadiet_Microbiota_PC1_PC2_Final, ncomp = 1) #Best "Figure.10a_DIABLO_Correlation_between_Matrices_By_Microbiota_Component_1"
                  plotDiablo(diablo.metadiet_Microbiota_PC1_PC2_Final, ncomp = 2) #"Figure.10b_DIABLO_Correlation_between_Matrices_By_Microbiota_Component_2"
        
                  plotDiablo(
                    diablo.metadiet_Microbiota_PC1_PC2_Final,
                    ncomp = 1,
                    col = c("darkblue", "#fc4e07", "forestgreen"),  # Cambia según tus grupos
                    cex = 5                                # Tamaño de los puntos
                  )
                  
                  
          #3.6.Plot individuals  
                  setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Individuals_By_Groups")
                  #Metadiet
              
                  plotIndiv(diablo.metadiet_Diet_Dendogram_Final, ind.names = FALSE, legend = TRUE, 
                            title = 'Figure.11_By_Plant_DIABLO_individuals_comp_1-2')
                  
                  plotIndiv(diablo.metadiet_Diet_Dendogram_Final, ind.names = FALSE, legend = TRUE, 
                            ncomp = 1,
                            col = c("darkblue", "#fc4e07"), 
                            cex = 3, 
                            pch = 19,
                            title = 'Figure.12_By_Metabolome_DIABLO_individuals_comp_1-2')
                  
                  
                  
                  plotIndiv(diablo.metadiet_Metabolome_PC1_PC2_Final, ind.names = FALSE, legend = TRUE, 
                            ncomp = 1,
                            col = c("darkblue", "#fc4e07"), 
                            cex = 3, 
                            pch = 19,
                            title = 'Figure.12_By_Metabolome_DIABLO_individuals_comp_1-2')
                  
                  
                  plotIndiv(diablo.metadiet_Microbiota_PC1_PC2_Final, ind.names = FALSE, legend = TRUE, 
                            title = 'Figure.13_By_Microbiota_DIABLO_individuals_comp_1-2')

                  plotIndiv(diablo.metadiet_Microbiota_PC1_PC2_Final, ind.names = FALSE, legend = TRUE, 
                            ncomp = 1,
                            col = c("darkblue", "#fc4e07", "forestgreen"), 
                            cex = 3, 
                            pch = 19,
                            title = 'Figure.12_By_Metabolome_DIABLO_individuals_comp_1-2')
                  
    
          #3.7.Plot Arrow 
                  
                  setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Arrows_by_groups")
                  
                  #Metadiet
                  plotArrow(diablo.metadiet_Diet_Dendogram_Final, comp = c(1,2), ind.names = TRUE,ind.names.size = 3, legend = TRUE,arrow.size = 0.5,arrow.length = 0, pch.size = 3, 
                            title = 'Figure.14_By_Plant_DIABLO_Arrow_comp_1-2')
                              

                  plotArrow(
                    diablo.metadiet_Diet_Dendogram_Final,
                    comp = c(1, 2),
                    ind.names = FALSE,
                    ind.names.size = 4,
                    legend = TRUE,
                    arrow.size = 1.5,
                    arrow.length = 0.2,
                    col.per.group = c("darkblue", "#fc4e07"),
                    pch.size = 5,
                    title = 'Figure.15_By_Metabolome_DIABLO_Arrow_comp_1-2'
                  )
                  
                  
                  
                  
                  
                  
                  
                  plotArrow(
                    diablo.metadiet_Metabolome_PC1_PC2_Final,
                    comp = c(1, 2),
                    ind.names = FALSE,
                    ind.names.size = 4,
                    legend = TRUE,
                    arrow.size = 1.5,
                    arrow.length = 0.2,
                    col.per.group = c("darkblue", "#fc4e07", "forestgreen"),
                    pch.size = 5,
                    title = 'Figure.15_By_Metabolome_DIABLO_Arrow_comp_1-2'
                  )
                  
                  
                  
                  
                  
                  plotArrow(diablo.metadiet_Microbiota_PC1_PC2_Final, comp = c(1,2), ind.names = TRUE,ind.names.size = 3, legend = TRUE,arrow.size = 0.5,arrow.length = 0, pch.size = 3,  
                            title = 'Figure.16_By_Microbiota_DIABLO_Arrow_comp_1-2')
              
                  
                  plotArrow(
                    diablo.metadiet_Microbiota_PC1_PC2_Final,
                    comp = c(1, 2),
                    ind.names = FALSE,
                    ind.names.size = 4,
                    legend = TRUE,
                    arrow.size = 1.5,
                    arrow.length = 0.2,
                    col.per.group = c("darkblue", "#fc4e07", "forestgreen"),
                    pch.size = 5,
                    title = 'Figure.15_By_Metabolome_DIABLO_Arrow_comp_1-2'
                  )
                  
                  
                              
                  #permanova 
                  library(vegan)
                              
                              # Get the coordinates of the samples in the space projected by DIABLO
                              coordinates_Arrow_Diet <- diablo.metadiet_Diet_Dendogram_Final$variates$Y  # Bloque X (u otro bloque si es necesario)
                              coordinates_Arrow_metabolome <- diablo.metadiet_Metabolome_PC1_PC2_Final$variates$Y  # Bloque X (u otro bloque si es necesario)
                              coordinates_Arrow_microbiota <- diablo.metadiet_Microbiota_PC1_PC2_Final$variates$Y  # Bloque X (u otro bloque si es necesario)
                              
                              # Permanova
                              setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Arrows_by_groups/DIABLO_Adonis_Arrow")
                              
                              adonis_Arrow_Diet_diablo <- adonis2(coordinates_Arrow_Diet ~ Group$Two_Diet.Type_Dendrogram, method = "euclidean")
                              print(adonis_Arrow_Diet_diablo)
                              write.table(adonis_Arrow_Diet_diablo,file = "adonis_Figure.14_By_Plant_DIABLO_Arrow_comp_1-2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                              
                              
                              adonis_Arrow_Metabolome_diablo <- adonis2(coordinates_Arrow_metabolome ~ Group$Two_Metabolome.type, method = "euclidean")
                              print(adonis_Arrow_Metabolome_diablo )
                              write.table(adonis_Arrow_Metabolome_diablo,file = "adonis_Figure.15_By_Metabolome_DIABLO_Arrow_comp_1-2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 

                              adonis_Arrow_Microbiota_diablo <- adonis2(coordinates_Arrow_microbiota ~ Group$Microbiota.profile, method = "euclidean")
                              print(adonis_Arrow_Microbiota_diablo)
                              write.table(adonis_Arrow_Microbiota_diablo,file = "adonis_Figure.16_By_Microbiota_DIABLO_Arrow_comp_1-2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                              

          #3.8.Variable Plot
                setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Variables/DIABLO_PlotVar_METADIET")
                #Metadiet
                  plotVar(diablo.metadiet_Diet_Dendogram_Final, var.names = TRUE, style = 'graphics', legend = TRUE, 
                          pch = c(16, 17, 15), cex = c(0.4,0.4,0.4), 
                          col = c('darkgreen','darkorange', 'darkblue'),
                          title = 'Figure.17_By_Plant_DIABLO_plotVar_comp_1-2')
                                selected_Diet_vars_comp_1 <- selectVar(diablo.metadiet_Diet_Dendogram_Final, comp = 1) 
                                    write.table(selected_Diet_vars_comp_1,file = "Variables_Figure.17_By_Plant_DIABLO_plotVar_comp_1.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                selected_Diet_vars_comp_2 <- selectVar(diablo.metadiet_Diet_Dendogram_Final, comp = 2) 
                                    write.table(selected_Diet_vars_comp_2,file = "Variables_Figure.17_By_Plant_DIABLO_plotVar_comp_2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 

                                #the same but more bacano
                                plotLoadings(diablo.metadiet_Diet_Dendogram_Final, comp = 1, contrib = 'max', method = 'median')#Figure.18a_By_Plant_DIABLO_plotLoadings_comp_1
                                plotLoadings(diablo.metadiet_Diet_Dendogram_Final, comp = 2, contrib = 'max', method = 'median')#Figure.18b_By_Plant_DIABLO_plotLoadings_comp_2
                                
                                
                  plotVar(diablo.metadiet_Metabolome_PC1_PC2_Final, var.names = TRUE, style = 'graphics', legend = TRUE, 
                          pch = c(16, 17, 15), cex = c(0.4,0.4,0.4), 
                          col = c('darkgreen','darkorange', 'darkblue'),
                          title = 'Figure.19_By_Metabolome_DIABLO_plotVar_comp_1-2')
                                 selected_Metabolome_vars_comp_1  <- selectVar(diablo.metadiet_Metabolome_PC1_PC2_Final, comp = 1)
                                    write.table(selected_Metabolome_vars_comp_1,file = "Variables_Figure.19_By_Metabolome_DIABLO_plotVar_comp_1.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                 selected_Metabolome_vars_comp_2  <- selectVar(diablo.metadiet_Metabolome_PC1_PC2_Final, comp = 2)
                                    write.table(selected_Metabolome_vars_comp_2,file = "Variables_Figure.19_By_Metabolome_DIABLO_plotVar_comp_2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                 
                                 #the same but more bacano
                                 plotLoadings(diablo.metadiet_Metabolome_PC1_PC2_Final, comp = 1, contrib = 'max', method = 'median')#Figure.20a_By_Metabolome_DIABLO_plotLoadings_comp_1
                                 plotLoadings(diablo.metadiet_Metabolome_PC1_PC2_Final, comp = 2, contrib = 'max', method = 'median')#Figure.20b_By_Metabolome_DIABLO_plotLoadings_comp_2
                                 
                                 PLOT_LOADINGS <-     plotLoadings(diablo.metadiet_Metabolome_PC1_PC2_Final, comp = 1, contrib = 'max', method = 'median')#Figure.20b_By_Metabolome_DIABLO_plotLoadings_comp_2

                                 
                  plotVar(diablo.metadiet_Microbiota_PC1_PC2_Final, var.names = TRUE, style = 'graphics', legend = TRUE, 
                          pch = c(16, 17, 15), cex = c(0.7,0.7,0.7), 
                          col = c('darkorchid', 'brown1', 'lightgreen'),
                          title = 'Figure.21_By_Microbiota_DIABLO_plotVar_comp_1-2')
          
                                  selected_Microbiota_vars_comp_1  <- selectVar(diablo.metadiet_Microbiota_PC1_PC2_Final, comp = 1)
                                      write.table(selected_Microbiota_vars_comp_1,file = "Variables_Figure.21_By_Microbiota_DIABLO_plotVar_comp_1.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                  selected_Microbiota_vars_comp_2  <- selectVar(diablo.metadiet_Microbiota_PC1_PC2_Final, comp = 2)
                                     write.table(selected_Microbiota_vars_comp_2,file = "Variables_Figure.19_By_Microbiota_DIABLO_plotVar_comp_1.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                  
                                  #the same but more bacano
                                  plotLoadings(diablo.metadiet_Microbiota_PC1_PC2_Final, comp = 1, contrib = 'max', method = 'median') #Figure.22a_By_Microbiota_DIABLO_plotLoadings_comp_1
                                  plotLoadings(diablo.metadiet_Microbiota_PC1_PC2_Final, comp = 2, contrib = 'max', method = 'median') #Figure.22b_By_Microbiota_DIABLO_plotLoadings_comp_2
                                  
                        
                                  
                                  
                                  # Verificar si hay loadings
                                  str(diablo.metadiet_Diet_Dendogram_Final$loadings)
                                  
                                  # O ver directamente los nombres de las variables por bloque
                                  names(diablo.metadiet_Diet_Dendogram_Final$loadings)
                                  
                                  # También puedes revisar si hay datos en los componentes
                                  head(diablo.metadiet_Diet_Dendogram_Final$loadings$Diet)
                                  head(diablo.metadiet_Diet_Dendogram_Final$loadings$Metabolome)
                                  head(diablo.metadiet_Diet_Dendogram_Final$loadings$Microbiome)
                                  
                                  
                                  
                                  
                                  
                                  
            #3.9.Heatmaps_and hearchical variables_and_individuals
                                  setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Heatmap")
                                  
                         
                        cimDiablo(diablo.metadiet_Metabolome_PC1_PC2_Final, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), #Figure.24.DIABLO_Heatmap_Correlations Between Metabolomics_M1 vs. M2_, Diet, and Microbiota _Comp 1_
                                  comp = 1, margin=c(8,20), legend.position = "topright")
                        

                        
                        cimDiablo(
                          diablo.metadiet_Metabolome_PC1_PC2_Final,
                          comp = 1,
                          color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
                          margin = c(8, 20),
                          legend.position = "topright",
                          row.cex = 0.7,   # tamaño del texto de las variables en filas
                          col.cex = 0.7    # tamaño del texto de las muestras (columnas)
                        )
                        

                        
                        
                        cimDiablo(diablo.metadiet_Microbiota_PC1_PC2_Final, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), #Figure.25.DIABLO_Heatmap_Correlations Between Microbiota_B1 vs. D2_ vs B3, Metabolomics, and Diet _Comp 1_
                                  comp = 1, margin=c(8,20), legend.position = "right")
                  
                        
                        
                        
                        Fig.24._heatmap_metabo[["col.names"]] <- cimDiablo(diablo.metadiet_Metabolome_PC1_PC2_Final, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), #Figure.24.DIABLO_Heatmap_Correlations Between Metabolomics_M1 vs. M2_, Diet, and Microbiota _Comp 1_
                                  comp = 1, margin=c(8,20), legend.position = "topright")
                        
                        Fig.25._heatmap_Micro <- cimDiablo(diablo.metadiet_Microbiota_PC1_PC2_Final, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), #Figure.25.DIABLO_Heatmap_Correlations Between Microbiota_B1 vs. D2_ vs B3, Metabolomics, and Diet _Comp 1_
                                  comp = 1, margin=c(8,20), legend.position = "right")
                        
                        
                     write.table( Fig.24._heatmap_metabo[["col.names"]],file = "Fig.24._heatmap_metabo.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                     write.table( Fig.25._heatmap_Micro[["col.names"]],file = "Fig.25._heatmap_Micro.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                     

                                  
            #3.10.circosPlot               
                    
                     setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Correlations_Circle_blocks")
                    #Metadiet
                    CircosPlot_Diet <-circosPlot(diablo.metadiet_Diet_Dendogram_Final, cutoff = 0.15, line = TRUE, 
                                      color.blocks = c('darkgreen',  'brown1','darkorchid' ),
                                      color.cor = c("chocolate3","grey20"), size.labels = 3, size.variables = 1, size.legend = 1.5,start.degree = 90)
                                      #Figure.26_By_Dieta_CircosPlot                       
                    
                                      #Correlation 
                                      write.table( CircosPlot_Diet,file = "Correlations_Figure.26_By_Dieta_CircosPlot.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                      
                    
                                      #Top correlations
                                      cor_data_CircosPlot_Diet <- as.data.frame(as.table(as.matrix(CircosPlot_Diet)))
                                      colnames(cor_data_CircosPlot_Diet) <- c("Variable1", "Variable2", "Correlation")
                                      cor_data_CircosPlot_Diet$Plants <- ifelse(cor_data_CircosPlot_Diet$Variable1 %in% rownames(CircosPlot_Diet)[1:44], "Plants",
                                                                                      ifelse(cor_data_CircosPlot_Diet$Variable1 %in% rownames(CircosPlot_Diet)[45:54], "Metabolites", "Microbiota"))
                                      cor_data_CircosPlot_Diet$Metabolites <- ifelse(cor_data_CircosPlot_Diet$Variable2 %in% rownames(CircosPlot_Diet)[1:44], "Plants",
                                                                                           ifelse(cor_data_CircosPlot_Diet$Variable2 %in% rownames(CircosPlot_Diet)[45:54], "Metabolites", "Microbiota"))
                                      filtered_cor_data_CircosPlot_Diet <- cor_data_CircosPlot_Diet %>%
                                        filter(Plants != Metabolites & Variable1 != Variable2)
                                      
                                      filtered_cor_data_CircosPlot_Diet <- filtered_cor_data_CircosPlot_Diet %>%
                                        mutate(PairID = ifelse(as.character(Variable1) < as.character(Variable2),
                                                               paste(Variable1, Variable2, sep = "_"),
                                                               paste(Variable2, Variable1, sep = "_"))) %>%
                                        distinct(PairID, .keep_all = TRUE) %>%
                                        select(-PairID) 
                                      # Rename columns
                                      colnames(filtered_cor_data_CircosPlot_Diet) <- 
                                        gsub("Plants", "Variable1 Block", colnames(filtered_cor_data_CircosPlot_Diet))
                                      colnames(filtered_cor_data_CircosPlot_Diet) <- 
                                        gsub("Metabolites", "Variable2 Block", colnames(filtered_cor_data_CircosPlot_Diet))
                                      filtered_cor_data_CircosPlot_Diet <- filtered_cor_data_CircosPlot_Diet[order(-abs(filtered_cor_data_CircosPlot_Diet$Correlation)), ]
                                      write.csv(filtered_cor_data_CircosPlot_Diet, "Top_cor_data_CircosPlot_Diet.csv", row.names = FALSE)
                                      
                                      
                                      
                                      
                                      
                    CircosPlot_Metabolome <-circosPlot(diablo.metadiet_Metabolome_PC1_PC2_Final, cutoff = 0.15, line = TRUE, 
                                            color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
                                            color.cor = c("chocolate3","grey20"), size.labels = 1, size.variables = 0.7, size.legend = 1.5,start.degree = 90)
                                            #Figure.27_By_Metabolome_CircosPlot   
                    
                    CircosPlot_Metabolome <- circosPlot(
                      diablo.metadiet_Metabolome_PC1_PC2_Final,
                      cutoff = 0.15,
                      cor.value = "positive",  # <--- SOLO correlaciones positivas
                      line = TRUE,
                      color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
                      color.cor = c("chocolate3", "grey20"),
                      size.labels = 1,
                      size.variables = 0.5,
                      size.legend = 1.5,
                      start.degree = 90
                    )
                    
                                      #Correlation 
                                      write.table( CircosPlot_Metabolome,file = "Correlations_Figure.27_By_Metabolome_CircosPlot.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
                                          
                                      #Top correlations
                                            cor_data_CircosPlot_Metabolome <- as.data.frame(as.table(as.matrix(CircosPlot_Metabolome)))
                                            colnames(cor_data_CircosPlot_Metabolome) <- c("Variable1", "Variable2", "Correlation")
                                            cor_data_CircosPlot_Metabolome$Plants <- ifelse(cor_data_CircosPlot_Metabolome$Variable1 %in% rownames(CircosPlot_Metabolome)[1:12], "Plants",
                                                                           ifelse(cor_data_CircosPlot_Metabolome$Variable1 %in% rownames(CircosPlot_Metabolome)[13:197], "Metabolites", "Microbiota"))
                                            cor_data_CircosPlot_Metabolome$Metabolites <- ifelse(cor_data_CircosPlot_Metabolome$Variable2 %in% rownames(CircosPlot_Metabolome)[1:12], "Plants",
                                                                           ifelse(cor_data_CircosPlot_Metabolome$Variable2 %in% rownames(CircosPlot_Metabolome)[13:197], "Metabolites", "Microbiota"))
                                            filtered_cor_data_CircosPlot_Metabolome <- cor_data_CircosPlot_Metabolome %>%
                                              filter(Plants != Metabolites & Variable1 != Variable2)
                                            
                                            filtered_cor_data_CircosPlot_Metabolome <- filtered_cor_data_CircosPlot_Metabolome %>%
                                              mutate(PairID = ifelse(as.character(Variable1) < as.character(Variable2),
                                                                     paste(Variable1, Variable2, sep = "_"),
                                                                     paste(Variable2, Variable1, sep = "_"))) %>%
                                              distinct(PairID, .keep_all = TRUE) %>%
                                              select(-PairID) 
                                                  # Rename columns
                                                  colnames(filtered_cor_data_CircosPlot_Metabolome) <- 
                                                    gsub("Plants", "Variable1 Block", colnames(filtered_cor_data_CircosPlot_Metabolome))
                                                  colnames(filtered_cor_data_CircosPlot_Metabolome) <- 
                                                    gsub("Metabolites", "Variable2 Block", colnames(filtered_cor_data_CircosPlot_Metabolome))
                                                  filtered_cor_data_CircosPlot_Metabolome <- filtered_cor_data_CircosPlot_Metabolome[order(-abs(filtered_cor_data_CircosPlot_Metabolome$Correlation)), ]
                                                  write.csv(filtered_cor_data_CircosPlot_Metabolome, "Top_cor_data_CircosPlot_Metabolome.csv", row.names = FALSE)

                                                  
                                                  




                                                  
                                            
            
                    CircosPlot_Microbiota <-circosPlot(diablo.metadiet_Microbiota_PC1_PC2_Final, cutoff = 0.5, line = TRUE, 
                                            color.blocks = c('darkgreen',  'brown1','darkorchid' ),
                                            color.cor = c("chocolate3","grey20"), size.labels = 3, size.variables = 0.7, size.legend = 1.5,start.degree = 90)
                                      
                                      #Figure.28_By_Microbiota_CircosPlot                       
                                      #Correlation 
                                      write.table( CircosPlot_Microbiota,file = "Correlations_Figure.28_By_Microbiota_CircosPlot.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
        
                                      
                                      
                                      #Top correlations
                                      cor_data_CircosPlot_Microbiota <- as.data.frame(as.table(as.matrix(CircosPlot_Microbiota)))
                                      colnames(cor_data_CircosPlot_Microbiota) <- c("Variable1", "Variable2", "Correlation")
                                      cor_data_CircosPlot_Microbiota$Plants <- ifelse(cor_data_CircosPlot_Microbiota$Variable1 %in% rownames(CircosPlot_Microbiota)[1:45], "Plants",
                                                                                ifelse(cor_data_CircosPlot_Microbiota$Variable1 %in% rownames(CircosPlot_Microbiota)[46:60], "Metabolites", "Microbiota"))
                                      cor_data_CircosPlot_Microbiota$Metabolites <- ifelse(cor_data_CircosPlot_Microbiota$Variable2 %in% rownames(CircosPlot_Microbiota)[1:45], "Plants",
                                                                                     ifelse(cor_data_CircosPlot_Microbiota$Variable2 %in% rownames(CircosPlot_Microbiota)[46:60], "Metabolites", "Microbiota"))
                                      filtered_cor_data_CircosPlot_Microbiota <- cor_data_CircosPlot_Microbiota %>%
                                        filter(Plants != Metabolites & Variable1 != Variable2)
                                      
                                      filtered_cor_data_CircosPlot_Microbiota <- filtered_cor_data_CircosPlot_Microbiota %>%
                                        mutate(PairID = ifelse(as.character(Variable1) < as.character(Variable2),
                                                               paste(Variable1, Variable2, sep = "_"),
                                                               paste(Variable2, Variable1, sep = "_"))) %>%
                                        distinct(PairID, .keep_all = TRUE) %>%
                                        select(-PairID) 
                                      # Rename columns
                                      colnames(filtered_cor_data_CircosPlot_Microbiota) <- 
                                        gsub("Plants", "Variable1 Block", colnames(filtered_cor_data_CircosPlot_Microbiota))
                                      colnames(filtered_cor_data_CircosPlot_Microbiota) <- 
                                        gsub("Metabolites", "Variable2 Block", colnames(filtered_cor_data_CircosPlot_Microbiota))
                                      filtered_cor_data_CircosPlot_Microbiota <- filtered_cor_data_CircosPlot_Microbiota[order(-abs(filtered_cor_data_CircosPlot_Microbiota$Correlation)), ]
                                      write.csv(filtered_cor_data_CircosPlot_Microbiota, "Top_cor_data_CircosPlot_Microbiota.csv", row.names = FALSE)
                                      
                                      
                                      
                                      
                        
                                      
                                      #Top variables and venn
                                      ## replace the interesting file from circosplot correlation
                                      
                                      #top 3 blocks
                                      
                                      library(dplyr)
                                      library(tidyr)
                                      library(readr)
                                      
                                      setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Correlations_Circle_blocks")
                                      # Cargar los datos desde el archivo Excel
                                      file_path <- "Correlations_Figure.26_By_Dieta_CircosPlot.xlsx"
                                      cor_data1 <- read.xlsx(file_path)
                                      
                                      
                                      # Definir los rangos de variables para cada bloque ### MIRAR CADA ARCHIVO XLSX########
                                      plants_range <- 1:44
                                      metabolites_range <- 45:54
                                      otus_range <- 55:102
                                      
                                      # Convertir a formato largo para filtrar más fácilmente
                                      cor_long <- cor_data1 %>%
                                        pivot_longer(cols = -1, names_to = "Variable2", values_to = "Correlation") %>%
                                        rename(Variable1 = 1)
                                      
                                      # Eliminar correlaciones intrabloque y mantener solo interbloques
                                      cor_long_filtered <- cor_long %>%
                                        filter(Correlation > 0.25) %>%  # Filtrar correlaciones positivas con R2 > 0.25
                                        mutate(Block1 = case_when(
                                          Variable1 %in% names(cor_data1)[plants_range] ~ "Plants",
                                          Variable1 %in% names(cor_data1)[metabolites_range] ~ "Metabolites",
                                          Variable1 %in% names(cor_data1)[otus_range] ~ "OTUs",
                                          TRUE ~ NA_character_
                                        )) %>%
                                        mutate(Block2 = case_when(
                                          Variable2 %in% names(cor_data1)[plants_range] ~ "Plants",
                                          Variable2 %in% names(cor_data1)[metabolites_range] ~ "Metabolites",
                                          Variable2 %in% names(cor_data1)[otus_range] ~ "OTUs",
                                          TRUE ~ NA_character_
                                        )) %>%
                                        filter(!is.na(Block1) & !is.na(Block2) & Block1 != Block2)  # Mantener solo interbloques
                                      
                                      # ⚠ MODIFICACIÓN: Evitar duplicados bidireccionales y ordenar de mayor a menor correlación
                                      cor_long_unique <- cor_long_filtered %>%
                                        mutate(Var1 = pmin(Variable1, Variable2),  # Ordenar nombres para evitar duplicados
                                               Var2 = pmax(Variable1, Variable2)) %>%
                                        distinct(Var1, Var2, .keep_all = TRUE) %>%  # Eliminar duplicados bidireccionales
                                        select(Var1, Var2, Correlation, Block1, Block2) %>%
                                        arrange(desc(Correlation))  # Ordenar de mayor a menor correlación
                                      
                                      # Guardar el archivo sin redundancias y ordenado
                                      write_csv(cor_long_unique, "Top_Correlations_3Blocks_Figure.26.csv")                                                  # Vista previa del resultado
                                      
                                      
                                      
                                      #top
                                      
                                      # Cargar paquetes necesarios
                                      library(dplyr)
                                      library(readr)
                                      library(tidyr)
                                      
                                      #Read the data
                                      top_correlations <- read.csv("Top_Correlations_3Blocks_Figure.27.csv")
                                      
                                      # Count the number of significant correlations for each variable
                                      correlation_counts <- top_correlations %>%
                                        select(Var1, Var2) %>%
                                        pivot_longer(cols = c(Var1, Var2), names_to = "Var_type", values_to = "Var") %>%
                                        group_by(Var) %>%
                                        summarise(Significant_Correlations = n()) %>%
                                        arrange(desc(Significant_Correlations))
                                      
                                      # Find related variables for each variable (Corrected)
                                      related_variables <- top_correlations %>%
                                        group_by(Var1) %>%
                                        summarise(Related_Variables = paste(Var2, collapse = ", "))
                                      
                                      # Rename Var1 to Var for merging
                                      related_variables <- related_variables %>%
                                        rename(Var = Var1)
                                      
                                      # Merge the two tables
                                      result_table <- merge(correlation_counts, related_variables, by = "Var") %>%
                                        select(Var, Significant_Correlations, Related_Variables)
                                      
                                      # Write the result to a new CSV file
                                      write.csv(result_table, "Variable_Correlation_Summary_figure.26.csv", row.names = FALSE)
                                      
                                      
                                      ## VENN
                                      #Var	Significant_Correlations	Related_Variables
                                      
                                      file_name6 <- "Variable_Correlation_Summary_Figure.27.xlsx"
                                      datavenn6 <- read.xlsx(file_name6, sheet = 1)  # Leer la primera hoja del archivo Excel
                                      
                                      # Seleccionar las METABOLITOS en las filas 2, 3 y 16
                                      Top_3_venn_diablo <- datavenn6[c(27, 28, 117), ]  # SELECCIONAR LAS VARIABLES DE INTERES
                                      
                                      # Separar las variables relacionadas para cada METABOLITO
                                      related_variables_list4 <- Top_3_venn_diablo %>%
                                        mutate(Variables = str_split(Related_Variables, ";")) %>%
                                        pull(Variables)  # Extraer la lista de variables relacionadas
                                      
                                      # Nombres de las Metabolites
                                      variable_names <- Top_3_venn_diablo$Var
                                      
                                      # Calcular las variables comunes entre las tres Metabolites
                                      common_variables_all <- Reduce(intersect, related_variables_list4)
                                      
                                      # Calcular las variables compartidas entre cada par de Metabolites
                                      pairwise_common_variables3 <- combn(seq_along(variable_names), 2, function(idx) {
                                        Variable.1 <- variable_names[idx[1]]
                                        Variable.2 <- variable_names[idx[2]]
                                        common_vars <- intersect(related_variables_list4[[idx[1]]], related_variables_list4[[idx[2]]])
                                        data.frame(METABOLYTE1 = Variable.1 , METABOLYTE2 = Variable.2, Variables_Comunes = paste(common_vars, collapse = "; "))
                                      }, simplify = FALSE)
                                      
                                      # Calcular las variables exclusivas para cada Metabolites
                                      exclusive_variables <- lapply(related_variables_list4, function(x) setdiff(x, common_variables_all))
                                      
                                      # Crear un archivo Excel
                                      output_file <- "Variables_Comunes_y_Exclusivas_DIABLO_Figure.26.xlsx"
                                      wb <- createWorkbook()
                                      
                                      # Agregar hoja para variables comunes entre las tres Metabolites
                                      addWorksheet(wb, "Comunes")
                                      writeData(wb, "Comunes", data.frame(Variables_Comunes = common_variables_all))
                                      
                                      # Agregar hoja para variables compartidas entre pares de Metabolites
                                      addWorksheet(wb, "Comunes1")
                                      pairwise_df <- do.call(rbind, pairwise_common_variables3)
                                      writeData(wb, "Comunes1", pairwise_df)
                                      
                                      # Agregar hojas para variables exclusivas por Metabolites
                                      for (i in seq_along(variable_names)) {
                                        sheet_name <- paste0("Ex_", variable_names[i])
                                        addWorksheet(wb, sheet_name)
                                        writeData(wb, sheet_name, data.frame(Variables_Exclusivas = exclusive_variables[[i]]))
                                      }
                                      
                                      # Guardar el archivo Excel
                                      saveWorkbook(wb, output_file, overwrite = TRUE)
                                      cat("Archivo exportado exitosamente como:", output_file, "\n")
                                      
                                      # Preparar los datos para el diagrama de Venn
                                      venn_data5 <- setNames(related_variables_list4, variable_names)
                                      
                                      # Crear el diagrama de Venn
                                      venn_plot5 <- ggvenn(venn_data5, fill_color = c("#FF9999", "#99CCFF", "#99FF99"))
                                      
                                      # Visualizar el diagrama en la ventana de gráficos
                                      print(venn_plot5)
                                      
                                      # Guardar el diagrama como archivo PNG
                                      output_venn4 <- "Venn_Diagrama_Metabolitos_Microbiota.png"
                                      ggsave(output_venn4, venn_plot4, width = 8, height = 6)
                                      cat("Diagrama de Venn guardado como:", output_venn4, "\n")
                                      
                                      
                                      
                                      
                                      
                                                    
        
                          
            #3.11.Networks              
            
                    #Metadiet
            
                    network(diablo.metadiet_Diet_Dendogram_Final, blocks = c(1,2,3), 
                            cutoff = 0.4,
                            color.node = c('#8a2be2','#c6e2ff', 'lightgreen'),
                            # To save the plot, uncomment below line
                            #save = 'png', name.save = 'diablo-network'
                            size.node = 0.03,
                            keysize.label = 0.3, cex.node.name =0.8, cex.edge.label = 0.3,color.edge = color.GreenRed(50))
                            #Figure.29_By_Diet_Network
                    
                    
                    network(diablo.metadiet_Metabolome_PC1_PC2_Final, blocks = c(1,2,3), 
                            cutoff = 0.4,
                            color.node = c('#8a2be2','#c6e2ff', 'lightgreen'),
                            # To save the plot, uncomment below line
                            #save = 'png', name.save = 'diablo-network'
                            size.node = 0.03,
                            keysize.label = 0.3, cex.node.name =0.5, cex.edge.label = 0.2,color.edge = color.GreenRed(50))
                            #Figure.30_By_Metabolome_Network
                    
                    
                    network(diablo.metadiet_Microbiota_PC1_PC2_Final, blocks = c(1,2,3), 
                            cutoff = 0.4,
                            color.node = c('#8a2be2','#c6e2ff', 'lightgreen'),
                            # To save the plot, uncomment below line
                            #save = 'png', name.save = 'diablo-network'
                            size.node = 0.05,
                            keysize.label = 0.3, cex.node.name =0.5, cex.edge.label = 0.2)
                            #Figure.31_By_Microbiota_Network
                    


#4. Univariate analysis - Correlations


                    Microbiota <- list (E_OTU_Microbiota_CLR[c(-1,-2)])
                    
                    

          #4.1 Data 
                PlantsCor <- data.frame (B_OTU_Plants_CLR[c(-1,-2)])
                MetabolitesCor <- data.frame (log10 (G_Metabolome_IsoMS[c(-1,-2)])) 
                MicrobiotaCor <- data.frame (E_OTU_Microbiota_CLR[c(-1,-2)])
                NutritionalCor <- data.frame (B_OTU_Plants_CLR_nutritional_class)
                
                library(caret)
                # To Identify plants with zero variance
                zero_var_plants <- nearZeroVar(PlantsCor, saveMetrics = TRUE)
                  # View plants with zero variance
                  print(zero_var_plants)
                  # Filter plants with zero variance and define it in PlantsCor
                  PlantsCor <- PlantsCor[, !zero_var_plants$nzv]
                  colnames(MetabolitesCor) 
          
          #4.1 normality test 
            
                  # Function to check normality
                  check_normality <- function(df) {
                  # Select only numeric columns
                  df <- df %>% select(where(is.numeric))
                  
                  # Apply Shapiro-Wilk Test to each variable
                  results <- df %>%
                  summarise(across(everything(), ~ shapiro.test(.)$p.value)) %>%
                  pivot_longer(cols = everything(), names_to = "Variable", values_to = "p_value") %>%
                  mutate(Normal_Distribution = if_else(p_value > 0.05, "Yes", "No"))
                  
                  return(results)
                  }
                  
                  # Apply to dataframes
                  normality_plants <- check_normality(PlantsCor)
        
                  # View results
                  normality_plants
                  normality_Metabolites
                  normality_microbiota
                  
                  # Export results to Excel
                  library(openxlsx)
                  write.xlsx(normality_plants, "Normality_Check_Plants.xlsx")
                  
                  ##fin normality test
                  
                  #Q-Qplot
                  library(ggplot2)
          
                  # Función para crear Q-Q plots
                  generate_qq_plot <- function(data, variable_name) {
                  ggplot(data, aes(sample = .data[[variable_name]])) +
                  stat_qq() +
                  stat_qq_line(color = "red") +
                  ggtitle(paste("Q-Q Plot:", variable_name)) +
                  theme_minimal()
                  }
                  
                  # Generar Q-Q plots para todas las variables numéricas
                  numeric_Metabolites_vars <- names(select(MetabolitesCor, where(is.numeric)))
                  qq_plot_Metabolites <- lapply(numeric_Metabolites_vars, function(var) generate_qq_plot(MetabolitesCor, var))
                  
                  # Mostrar los gráficos
                  qq_plot_Metabolites[[261]]  # Mostrar el primero como ejemplo
                  
          
          #4.2 Correlations
              #Due transformation CLRs and log in databases, we assumed normal distribution and applied pearson correlation
              #Top 1000 correlations between plants and metabolome
              #top 10000 correlations between plants and microbiota
              #top 10000 correlations between metabolome and microbiota
          
                library(tidyverse)
                library(openxlsx)
                  setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/Correlations_Pearson")
                
                # Función para calcular correlaciones, incluir p-values, seleccionar las mil más fuertes, calcular FDR y guardar resultados
                calculate_and_save_top_correlations_with_fdr <- function(df1, df2, name1, name2) {
                # Seleccionar solo columnas numéricas
                df1 <- df1 %>% select(where(is.numeric))
                df2 <- df2 %>% select(where(is.numeric))
                
                # Calcular las correlaciones y los valores p
                results <- expand.grid(
                Var1 = colnames(df1),
                Var2 = colnames(df2)
                ) %>%
                mutate(
                  Correlation = map2_dbl(
                    Var1, Var2,
                    ~ cor(df1[[.x]], df2[[.y]], method = "pearson", use = "pairwise.complete.obs")
                  ),
                  P_value = map2_dbl(
                    Var1, Var2,
                    ~ cor.test(df1[[.x]], df2[[.y]], method = "pearson")$p.value
                  )
                ) %>%
                arrange(desc(abs(Correlation))) %>% # Ordenar por correlaciones más fuertes (absolutas)
                slice(1:10000) %>% # Seleccionar las primeras mil correlaciones más fuertes
                mutate(
                  FDR = p.adjust(P_value, method = "BH") # Ajustar los valores p seleccionados usando Benjamini-Hochberg
                )
                
                # Exportar resultados a un archivo Excel
                file_name <- paste0("Top_1000_Correlaciones_", name1, "_vs_", name2, "_with_FDR.xlsx")
                write.xlsx(results, file_name)
                
                return(results)
                }
                
                # Realizar las correlaciones con FDR para los 1000 valores más fuertes y guardar los resultados
                results_plants_Metabolites_fdr <- calculate_and_save_top_correlations_with_fdr(
                PlantsCor, MetabolitesCor, "Plants", "Metabolites"
                )
                
                results_plants_microbiota_pearson_fdr <- calculate_and_save_top_correlations_with_fdr(
                PlantsCor, MicrobiotaCor, "Plants", "Microbiota"
                )
                
                results_Metabolites_microbiota_pearson_fdr <- calculate_and_save_top_correlations_with_fdr(
                MetabolitesCor, MicrobiotaCor, "Metabolites", "Microbiota"
                )
      
                
                plantas_interes <- c("s__Theobroma_cacao", "g__Ficus", "g__Passiflora", "g__Sorghum","sf__Nepetoideae.1", "s__Solanum_lycopersicum", "s__Ananas_comosus", "s__Ocimum_basilicum","s__Vicia_faba")
                

                # Usamos solo las correlaciones significativas y positivas de tus plantas
                filtered_data <- results_plants_Metabolites_fdr %>%
                  filter(Var1 %in% plantas_interes, FDR < 0.05, Correlation > 0)
                
                # Graficamos con ggplot como un heatmap extendido
                ggplot(filtered_data, aes(x = Var2, y = Var1, fill = Correlation)) +
                  geom_tile() +
                  scale_fill_gradient2(low = "white", high = "red", midpoint = 0,
                                       limit = c(0, 1), space = "Lab",
                                       name = "Correlación") +
                  theme_minimal() +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                  labs(title = "Correlaciones positivas significativas (FDR < 0.05)",
                       x = "Metabolitos", y = "Plantas")
                
                
                plantas_interes <- c(
                  "s__Theobroma_cacao", "g__Ficus", "g__Passiflora", "s__Vicia_faba", "g__Sorghum",
                  "s__Solanum_lycopersicum", "s__Coriandrum_sativum", "t__Triticeae", "s__Spinacia_tetrandra",
                  "s__Oryza_sativa.CP", "sf__Brassiceae", "s__Cocos_nucifera", "g__Musa", "f__Cucurbitaceae.2",
                  "sf__Solanoideae", "g__Phaseolus", "s__Allium_rubellum", "s__Pimpinella_anisum", "g__Oryza",
                  "t__Anserineae", "f__Lauraceae", "g__Solanum", "g__Allium", "f__Anacardiaceae", "s__Carica_papaya",
                  "g__Mentha.2", "s__Linum_usitatissimum", "s__Annona.muricata", "s__Brassica_oleracea.AB", "g__Cissus"
                )
                
                # Usamos solo las correlaciones significativas y positivas de tus plantas
                filtered_data <- results_plants_Metabolites_fdr %>%
                  filter(Var1 %in% plantas_interes, FDR < 0.05, Correlation > 0.39)
                

                # Ordenamos Var2 por la correlación promedio con las plantas de interés
                ordered_vars <- filtered_data %>%
                  group_by(Var2) %>%
                  summarise(mean_corr = mean(Correlation, na.rm = TRUE)) %>%
                  arrange(desc(mean_corr)) %>%
                  pull(Var2)
                
                # Convertimos Var2 en factor con el nuevo orden
                filtered_data$Var2 <- factor(filtered_data$Var2, levels = ordered_vars)
                
                # Graficamos con leyendas en inglés
                ggplot(filtered_data, aes(x = Var2, y = Var1, fill = Correlation)) +
                  geom_tile() +
                  scale_fill_gradient(low = "lightblue", high = "darkred", name = "Correlation") +
                  theme_minimal() +
                  theme(
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),  # Cambia el tamaño aquí
                    axis.text.y = element_text(size = 10) 
                    ) +
                  labs(
                    title = "Significant Positive Correlations (FDR < 0.05)",
                    x = "Metabolites",
                    y = "Plants"
                  )
                
                #en mi orden
                
                plantas_interes <- c(
                  # Vegetables
                  "f__Cucurbitaceae.2",
                  "sf__Brassiceae",
                  "t__Anserineae",
                  "g__Allium",
                  "g__Solanum",
                  "s__Allium_rubellum",
                  "s__Brassica_oleracea.AB",
                  "s__Solanum_lycopersicum",
                  "s__Spinacia_tetrandra",
                  
                  # Fruits
                  "f__Anacardiaceae",
                  "sf__Solanoideae",
                  "g__Cissus",
                  "g__Ficus",
                  "g__Musa",
                  "g__Passiflora",
                  "s__Annona.muricata",
                  "s__Carica_papaya",
                  
                  # plant-based fats
                  "s__Cocos_nucifera",
                  "s__Linum_usitatissimum",
                  "s__Theobroma_cacao",
                  
                  # Aromatic, spices and coffee
                  "f__Lauraceae",
                  "g__Mentha.2",
                  "s__Coriandrum_sativum",
                  "s__Pimpinella_anisum",
                  
                  # Cereals
                  "t__Triticeae",
                  "g__Oryza",
                  "s__Oryza_sativa.CP",
                  "g__Sorghum",
                  
                  # Legumes
                  "g__Phaseolus",
                  "s__Vicia_faba"
                )
                
                
                filtered_data <- results_plants_Metabolites_fdr %>%
                  filter(Var1 %in% plantas_interes, FDR < 0.05, Correlation > 0.39)
                
                # Orden de metabolitos por correlación
                ordered_vars <- filtered_data %>%
                  group_by(Var2) %>%
                  summarise(mean_corr = mean(Correlation, na.rm = TRUE)) %>%
                  arrange(desc(mean_corr)) %>%
                  pull(Var2)
                
                # Ordenar factores
                filtered_data$Var2 <- factor(filtered_data$Var2, levels = ordered_vars)
                filtered_data$Var1 <- factor(filtered_data$Var1, levels = plantas_interes)
                
                # Plot
                ggplot(filtered_data, aes(x = Var2, y = Var1, fill = Correlation)) +
                  geom_tile() +
                  scale_fill_gradient(low = "lightblue", high = "darkred", name = "Correlation") +
                  theme_minimal() +
                  theme(
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
                    axis.text.y = element_text(size = 10)
                  ) +
                  labs(
                    title = "Significant Positive Correlations (FDR < 0.05)",
                    x = "Metabolites",
                    y = "Plants"
                  )
                
                
                ### opciones visualizacion
                
                plant_nutritional_class <- data.frame(
                  Var1 = plantas_interes,
                  NutritionalGroup = rep(c(
                    "Vegetables", "Fruits", "Plant-based fats", 
                    "Aromatic/spices", "Cereals", "Legumes"
                  ), times = c(9, 8, 3, 4, 4, 2))
                )
                
                filtered_data <- filtered_data %>%
                  left_join(plant_nutritional_class, by = "Var1")
                
                ggplot(filtered_data, aes(x = Var2, y = Var1, fill = Correlation)) +
                  geom_tile() +
                  scale_fill_gradient(low = "lightblue", high = "darkred", name = "Correlation") +
                  facet_grid(NutritionalGroup ~ ., scales = "free_y", space = "free_y") +
                  theme_minimal() +
                  theme(
                    strip.text.y = element_text(angle = 0, size = 12),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    axis.text.y = element_text(size = 10)
                  ) +
                  labs(
                    title = "Significant Positive Correlations (FDR < 0.05)",
                    x = "Metabolites",
                    y = "Plants"
                  )
                
                
                #######
                #metabolitos microbiota
                
                library(tidyverse)
                library(readxl)
                library(readr)
                library(ggplot2)
                
                # 1. Leer resultados de correlación
                cor_results <- read.xlsx("Top_1000_Correlaciones_Metabolites_vs_Microbiota_with_FDR.xlsx")
                
                # 2. Filtrar solo correlaciones positivas y con FDR significativo
                cor_pos_sig <- cor_results %>%
                  filter(Correlation > 0, FDR < 0.05)
                
                # 3. Leer archivo de taxonomía
                taxonomy <- F_Taxonomy_Microbiota
                
                # 4. Unir taxonomía con correlaciones
                cor_pos_sig_tax <- cor_pos_sig %>%
                  left_join(taxonomy, by = c("Var2" = "OTU"))
                
                # 5. Contar correlaciones positivas por metabolito
                top_metabolites <- cor_pos_sig_tax %>%
                  group_by(Var1) %>%
                  summarise(n_corr = n()) %>%
                  arrange(desc(n_corr)) %>%
                  slice_head(n = 50) %>%
                  pull(Var1)
                
                # 6. Filtrar solo esos 50 metabolitos
                cor_top50 <- cor_pos_sig_tax %>%
                  filter(Var1 %in% top_metabolites)
                
                # 7. Gráfico
                # 1. Contar cuántas correlaciones hay por metabolito y género
                cor_top50_summary <- cor_top50 %>%
                  group_by(Var1, Genus) %>%
                  summarise(n = n(), .groups = "drop")
                
                # 2. Reordenar metabolitos por total de correlaciones
                cor_top50_summary <- cor_top50_summary %>%
                  mutate(Var1 = fct_reorder(Var1, n, .fun = sum))
                
                # 3. Graficar
                ggplot(cor_top50_summary, aes(x = Var1, y = n, fill = Genus)) +
                  geom_col() +
                  coord_flip() +
                  theme_minimal() +
                  labs(
                    title = "Genera from the Lachnospiraceae family",
                    x = "Metabolyte",
                    y = "# of positive correlations",
                    fill = "Genus"
                  ) +
                  theme(
                    axis.text.y = element_text(size = 4),
                    legend.position = "right"
                  )
                
          
                ######      
                ##Gráfico excluyendo Lachnospiraceae
                
                library(viridis)
                cor_top50_no_lachno <- cor_top50 %>%
                  filter(Family != "Lachnospiraceae") %>%
                  group_by(Var1, Genus) %>%
                  filter(n() > 3) %>%  # Excluye pares que aparecen solo una vez
                  ungroup() %>%
                  group_by(Var1, Genus) %>%
                  summarise(n = n(), .groups = "drop") %>%
                  mutate(Var1 = fct_reorder(Var1, n, .fun = sum))
                
                ggplot(cor_top50_no_lachno, aes(x = Var1, y = n, fill = Genus)) +
                  geom_col() +
                  coord_flip() +
                  theme_minimal() +
                  labs(
                    title = "Non-Lachnospiraceae genera",
                    x = "Metabolyte",
                    y = "# of positive correlations",
                    fill = "Genus"
                  ) +
                  theme(
                    axis.text.y = element_text(size = 8),
                    legend.position = "right"
                  )

                write.xlsx(cor_top50_no_lachno, "cor_top50_no_lachno.xlsx")
                
                
                                ######
                
                #2. Gráfico solo con Lachnospiraceae
                cor_top50_lachno <- cor_top50 %>%
                  filter(Family == "Lachnospiraceae") %>%
                  group_by(Var1, Genus) %>%
                  filter(n() > 2) %>%
                  ungroup() %>%
                  group_by(Var1, Genus) %>%
                  summarise(n = n(), .groups = "drop") %>%
                  mutate(Var1 = fct_reorder(Var1, n, .fun = sum))
                
                ggplot(cor_top50_lachno, aes(x = Var1, y = n, fill = Genus)) +
                  geom_col() +
                  coord_flip() +
                  theme_minimal() +
                  labs(
                    title = "Genera from the Lachnospiraceae family",
                    x = "Metabolyte",
                    y = "# of positive correlations",
                    fill = "Genus"
                  ) +
                  theme(
                    axis.text.y = element_text(size = 8),
                    legend.position = "right"
                  )
                
                
                write.xlsx(cor_top50_lachno, "cor_top50_lachno.xlsx")
                
                
                ###
                ######
                library(svglite)
                
                # OTUs de interés
                otus_interes <- c(
                  "Otu1358", "Otu1961", "Otu0071", "Otu0423", 
                  "Otu0024", "Otu0131", "Otu0482", "Otu0038", "Otu0361",
                  "Otu0017", "Otu0052", "Otu1261", "Otu0358"
                )                
                # Filtramos correlaciones significativas y positivas para esos OTUs
                filtered_data_otus <- results_Metabolites_microbiota_pearson_fdr %>%
                  filter(Var2 %in% otus_interes, FDR < 0.05, Correlation > 0)
                
                # Ordenamos metabolitos por la correlación promedio con los OTUs
                ordered_metabolites <- filtered_data_otus %>%
                  group_by(Var1) %>%
                  summarise(mean_corr = mean(Correlation, na.rm = TRUE)) %>%
                  arrange(desc(mean_corr)) %>%
                  pull(Var1)
                
                # Convertimos Var1 (metabolitos) en factor con el nuevo orden
                filtered_data_otus$Var1 <- factor(filtered_data_otus$Var1, levels = ordered_metabolites)
                
                # Graficamos con metabolitos en el eje X
                ggplot(filtered_data_otus, aes(x = Var1, y = Var2, fill = Correlation)) +
                  geom_tile() +
                  scale_fill_gradient(low = "lightblue", high = "darkred", name = "Correlation") +
                  theme_minimal() +
                  theme(
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    axis.text.y = element_text(size = 10)
                  ) +
                  labs(
                    title = "Significant Positive Correlations (FDR < 0.05)",
                    x = "Metabolites",
                    y = "Microbiota (OTUs)"
                  )
                
                
                #burbuja
                

                library(ggplot2)
                library(dplyr)
                
                # Contamos OTUs por metabolito y filtramos los que tienen más de 5
                otu_counts <- results_Metabolites_microbiota_pearson_fdr %>%
                  group_by(Var1) %>%
                  summarise(n_otus = n()) %>%
                  filter(n_otus > 3) %>%
                  arrange(n_otus)
                
                # Factor ordenado para mantener el orden en el eje Y
                otu_counts$Var1 <- factor(otu_counts$Var1, levels = otu_counts$Var1)
                
                # Gráfico: OTUs en X, Metabolitos en Y, círculos azules con número dentro
                ggplot(otu_counts, aes(x = n_otus, y = Var1)) +
                  geom_point(aes(size = n_otus), shape = 21, fill = "#4682B4", color = "black", alpha = 0.9) +
                  geom_text(aes(label = n_otus), color = "white", size = 3.5) +
                  scale_size(range = c(4, 15), guide = FALSE) +
                  theme_minimal() +
                  theme(
                    axis.text.y = element_text(size = 10),
                    axis.title.x = element_text(size = 12),
                    plot.title = element_text(hjust = 0.5)
                  ) +
                  labs(
                    title = "Number of OTUs Correlated per Metabolite (r > 0.5, n > 5)",
                    x = "Number of OTUs",
                    y = "Metabolites"
                  )
                
                
                
                ##Nutritional
                library(tidyverse)
                library(openxlsx)
                
                setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/Correlations_Pearson")
                
                # Función para calcular correlaciones, incluir p-values, seleccionar las mil más fuertes, calcular FDR y guardar resultados
                calculate_and_save_top_correlations_with_fdr <- function(df1, df2, name1, name2) {
                  # Seleccionar solo columnas numéricas
                  df1 <- df1 %>% select(where(is.numeric))
                  df2 <- df2 %>% select(where(is.numeric))
                  
                  # Calcular las correlaciones y los valores p
                  results <- expand.grid(
                    Var1 = colnames(df1),
                    Var2 = colnames(df2)
                  ) %>%
                    mutate(
                      Correlation = map2_dbl(
                        Var1, Var2,
                        ~ cor(df1[[.x]], df2[[.y]], method = "pearson", use = "pairwise.complete.obs")
                      ),
                      P_value = map2_dbl(
                        Var1, Var2,
                        ~ cor.test(df1[[.x]], df2[[.y]], method = "pearson")$p.value
                      )
                    ) %>%
                    arrange(desc(abs(Correlation))) %>% # Ordenar por correlaciones más fuertes (absolutas)
                    slice(1:100) %>% # Seleccionar las primeras mil correlaciones más fuertes
                    mutate(
                      FDR = p.adjust(P_value, method = "BH") # Ajustar los valores p seleccionados usando Benjamini-Hochberg
                    )
                  
                  # Exportar resultados a un archivo Excel
                  file_name <- paste0("Top_1000_Correlaciones_", name1, "_vs_", name2, "_with_FDR.xlsx")
                  write.xlsx(results, file_name)
                  
                  return(results)
                }
                
                # Realizar las correlaciones con FDR para los 1000 valores más fuertes y guardar los resultados (versión Nutritional)
                results_nutritional_Metabolites_fdr <- calculate_and_save_top_correlations_with_fdr(
                  NutritionalCor, MetabolitesCor, "Nutritional", "Metabolites"
                )
                
                results_nutritional_microbiota_pearson_fdr <- calculate_and_save_top_correlations_with_fdr(
                  NutritionalCor, MicrobiotaCor, "Nutritional", "Microbiota"
                )
                
                
                
                    ###plot
                # Definimos los grupos de interés nutricional
                nutritional_interest <- c(
                  "Aromatic..spices.and.coffee",
                  "Cereals",
                  "Fruits",
                  "Legumes",
                  "Vegetables",
                  "plant.based.fats"
                )                
                # Filtramos correlaciones significativas y positivas para los grupos de interés
                filtered_data <- results_nutritional_Metabolites_fdr %>%
                  filter(Var1 %in% nutritional_interest, FDR < 0.05, Correlation > 0)
                
                # Ordenamos los metabolitos (Var2) por la correlación promedio con los grupos de interés
                ordered_vars <- filtered_data %>%
                  group_by(Var2) %>%
                  summarise(mean_corr = mean(Correlation, na.rm = TRUE)) %>%
                  arrange(desc(mean_corr)) %>%
                  pull(Var2)
                
                # Convertimos Var2 en factor con el nuevo orden
                filtered_data$Var2 <- factor(filtered_data$Var2, levels = ordered_vars)
                
                # Graficamos
                ggplot(filtered_data, aes(x = Var2, y = Var1, fill = Correlation)) +
                  geom_tile() +
                  scale_fill_gradient(low = "lightblue", high = "darkred", name = "Correlation") +
                  theme_minimal() +
                  theme(
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    axis.text.y = element_text(size = 10)
                  ) +
                  labs(
                    title = "Significant Positive Correlations (FDR < 0.05)",
                    x = "Metabolites",
                    y = "Nutritional Groups"
                  )
                
                
                
            #4.3 Top variables with highest correlations with other variables.
                
                calculate_top_variables_custom <- function(df1, df2, name1, name2, p_threshold = 0.05, top_n_df1, top_n_df2) {
                  # Seleccionar solo columnas numéricas
                  df1 <- df1 %>% select(where(is.numeric))
                  df2 <- df2 %>% select(where(is.numeric))
                  
                  # Calcular las correlaciones y los valores p
                  results <- expand.grid(
                    Var1 = colnames(df1),
                    Var2 = colnames(df2)
                  ) %>%
                    mutate(
                      Correlation = map2_dbl(
                        Var1, Var2,
                        ~ cor(df1[[.x]], df2[[.y]], method = "pearson", use = "pairwise.complete.obs")
                      ),
                      P_value = map2_dbl(
                        Var1, Var2,
                        ~ cor.test(df1[[.x]], df2[[.y]], method = "pearson")$p.value
                      )
                    )
                  
                  # Filtrar correlaciones significativas por valor p
                  significant_results <- results %>%
                    filter(P_value < p_threshold, , Correlation > 0)
                  
                  # Asegurarse de que las columnas estén en formato correcto
                  significant_results <- significant_results %>%
                    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
                  
                  # Contar correlaciones significativas por variable en df1 y obtener las variables relacionadas
                  top_variables_df1 <- significant_results %>%
                    group_by(Var1) %>%
                    summarise(
                      Significant_Correlations = n(),
                      Related_Variables = paste(Var2, collapse = "; ")
                    ) %>%
                    arrange(desc(Significant_Correlations)) %>%
                    slice(1:top_n_df1)
                  
                  # Contar correlaciones significativas por variable en df2 y obtener las variables relacionadas
                  top_variables_df2 <- significant_results %>%
                    group_by(Var2) %>%
                    summarise(
                      Significant_Correlations = n(),
                      Related_Variables = paste(Var1, collapse = "; ")
                    ) %>%
                    arrange(desc(Significant_Correlations)) %>%
                    slice(1:top_n_df2)
                  
                  # Guardar resultados en Excel
                  file_name_df1 <- paste0("Top_", top_n_df1, "_Variables_", name1, "_Significant_Correlations_with_", name2, ".xlsx")
                  file_name_df2 <- paste0("Top_", top_n_df2, "_Variables_", name2, "_Significant_Correlations_with_", name1, ".xlsx")
                  
                  write.xlsx(top_variables_df1, file_name_df1)
                  write.xlsx(top_variables_df2, file_name_df2)
                  
                  return(list(Top_Variables_df1 = top_variables_df1, Top_Variables_df2 = top_variables_df2))
                }
                
                # Aplicar la función con los valores específicos para tus necesidades
                # 1. Top 15 para plantas con metaboloma
                top_plants_metabolites <- calculate_top_variables_custom(
                  PlantsCor, MetabolitesCor, "Plants", "Metabolites", 
                  top_n_df1 = 48, top_n_df2 = 50
                )
                
                # 2. Top 15 para plantas con microbiota
                top_plants_microbiota <- calculate_top_variables_custom(
                  PlantsCor, MicrobiotaCor, "Plants", "Microbiota", 
                  top_n_df1 = 48, top_n_df2 = 500
                )
                
                # 3. Top 100 para microbiota con metaboloma
                top_microbiota_metabolites <- calculate_top_variables_custom(
                  MicrobiotaCor, MetabolitesCor, "Microbiota", "Metabolites", 
                  top_n_df1 = 500, top_n_df2 = 50
                  
                  
                )
                
                
                
                
 
                
            #4.4 VennDiagram
                
                setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/Correlations_Pearson")

                library(openxlsx)
                library(dplyr)
                library(stringr)
                library(ggvenn)
                
                # Diet and metabolites
                # Leer el archivo de Excel
                file_name <- "Top_48_Variables_Plants_Significant_Correlations_with_Metabolites.xlsx"
                datavenn <- read.xlsx(file_name, sheet = 1)  # Leer la primera hoja del archivo Excel
                
                # Seleccionar las plantas en las filas 1, 7 y 9 (modificación realizada aquí)
                selected_plants <- datavenn[c(1,4, 7), ]  # Seleccionar las plantas de las filas 1, 7 y 9
                
                # Separar las variables relacionadas para cada planta
                related_variables_list <- selected_plants %>%
                  mutate(Variables = str_split(Related_Variables, ";")) %>%
                  pull(Variables)  # Extraer la lista de variables relacionadas
                
                # Nombres de las plantas seleccionadas
                plant_names <- selected_plants$Var1
                
                # Calcular las variables comunes entre las tres plantas
                common_variables_all <- Reduce(intersect, related_variables_list)
                
                # Calcular las variables compartidas entre cada par de plantas
                pairwise_common_variables <- combn(seq_along(plant_names), 2, function(idx) {
                  plant1 <- plant_names[idx[1]]
                  plant2 <- plant_names[idx[2]]
                  common_vars <- intersect(related_variables_list[[idx[1]]], related_variables_list[[idx[2]]])
                  data.frame(Planta1 = plant1, Planta2 = plant2, Variables_Comunes = paste(common_vars, collapse = "; "))
                }, simplify = FALSE)
                
                # Calcular las variables exclusivas para cada planta
                exclusive_variables <- lapply(related_variables_list, function(x) setdiff(x, common_variables_all))
                
                # Crear un archivo Excel
                output_file <- "Variables_Comunes_y_Exclusivas_Ampliado.xlsx"
                wb <- createWorkbook()
                
                # Agregar hoja para variables comunes entre las tres plantas
                addWorksheet(wb, "Comunes_3_Plantas")
                writeData(wb, "Comunes_3_Plantas", data.frame(Variables_Comunes = common_variables_all))
                
                # Agregar hoja para variables compartidas entre pares de plantas
                addWorksheet(wb, "Comunes_2_Plantas")
                pairwise_df <- do.call(rbind, pairwise_common_variables)
                writeData(wb, "Comunes_2_Plantas", pairwise_df)
                
                # Agregar hojas para variables exclusivas por planta
                for (i in seq_along(plant_names)) {
                  sheet_name <- paste0("Ex_", plant_names[i])
                  addWorksheet(wb, sheet_name)
                  writeData(wb, sheet_name, data.frame(Variables_Exclusivas = exclusive_variables[[i]]))
                }
                
                # Guardar el archivo Excel
                saveWorkbook(wb, output_file, overwrite = TRUE)
                
                cat("Archivo exportado exitosamente como:", output_file, "\n")
                
                # Preparar los datos para el diagrama de Venn
                venn_data1 <- setNames(related_variables_list, plant_names)
                
                # Crear el diagrama de Venn
                venn_plot1 <- ggvenn(venn_data1, fill_color = c("#FF9999", "#99CCFF", "#99FF99"))
                
                # Visualizar el diagrama en la ventana de gráficos
                print(venn_plot1)
                
                # Guardar el diagrama como archivo PNG
                output_venn1 <- "Venn_Diagrama_Plantas_Metabolites.png"
                ggsave(output_venn1, venn_plot1, width = 8, height = 6)
                cat("Diagrama de Venn guardado como:", output_venn1, "\n")
                
                
                
            #Diet and microbiota
                    
                    
                    # Leer el archivo de Excel
                    file_name3 <- "Top_48_Variables_Plants_Significant_Correlations_with_Microbiota.xlsx"
                    datavenn3 <- read.xlsx(file_name3, sheet = 1)  # Leer la primera hoja del archivo Excel
                    
                    # Seleccionar las 3 plantas con más correlaciones significativas
                    top_3_plants1 <- datavenn3 %>%
                      arrange(desc(Significant_Correlations)) %>%
                      slice(1:3)  # Seleccionar las 3 primeras
                    
                    # Separar las variables relacionadas para cada planta
                    related_variables_list1 <- top_3_plants1 %>%
                      mutate(Variables = str_split(Related_Variables, ";")) %>%
                      pull(Variables)  # Extraer la lista de variables relacionadas
                    
                    # Nombres de las plantas
                    plant_names1 <- top_3_plants1$Var1
                    
                    # Calcular las variables comunes entre las tres plantas
                    common_variables_all <- Reduce(intersect, related_variables_list1)
                    
                    # Calcular las variables compartidas entre cada par de plantas
                    pairwise_common_variables <- combn(seq_along(plant_names1), 2, function(idx) {
                      plant1 <- plant_names1[idx[1]]
                      plant2 <- plant_names1[idx[2]]
                      common_vars <- intersect(related_variables_list1[[idx[1]]], related_variables_list1[[idx[2]]])
                      data.frame(Planta1 = plant1, Planta2 = plant2, Variables_Comunes = paste(common_vars, collapse = "; "))
                    }, simplify = FALSE)
                    
                    # Calcular las variables exclusivas para cada planta
                    exclusive_variables <- lapply(related_variables_list1, function(x) setdiff(x, common_variables_all))
                    
                    # Crear un archivo Excel
                    output_file <- "Variables_Comunes_y_Exclusivas_Microbiota.xlsx"
                    wb <- createWorkbook()
                    
                    # Agregar hoja para variables comunes entre las tres plantas
                    addWorksheet(wb, "Comunes_4_Plantas")
                    writeData(wb, "Comunes_4_Plantas", data.frame(Variables_Comunes = common_variables_all))
                    
                    # Agregar hoja para variables compartidas entre pares de plantas
                    addWorksheet(wb, "Comunes_3_Plantas")
                    pairwise_df <- do.call(rbind, pairwise_common_variables)
                    writeData(wb, "Comunes_3_Plantas", pairwise_df)
                    
                    # Agregar hojas para variables exclusivas por planta
                    for (i in seq_along(plant_names)) {
                      sheet_name <- paste0("Exclusivas_", plant_names[i])
                      addWorksheet(wb, sheet_name)
                      writeData(wb, sheet_name, data.frame(Variables_Exclusivas = exclusive_variables[[i]]))
                    }
                    
                    # Guardar el archivo Excel
                    saveWorkbook(wb, output_file, overwrite = TRUE)
                    cat("Archivo exportado exitosamente como:", output_file, "\n")
                    
                    # Preparar los datos para el diagrama de Venn
                    venn_data2 <- setNames(related_variables_list1, plant_names1)
                    
                    # Crear el diagrama de Venn
                    venn_plot2 <- ggvenn(venn_data2, fill_color = c("#FF9999", "#99CCFF", "#99FF99"))
                    
                    # Visualizar el diagrama en la ventana de gráficos
                    print(venn_plot2)
                    
                    # Guardar el diagrama como archivo PNG
                    output_venn2 <- "Venn_Diagrama_Plantas_Microbiota.png"
                    ggsave(output_venn2, venn_plot2, width = 8, height = 6)
                    cat("Diagrama de Venn guardado como:", output_venn2, "\n")
                    
                    
                    
                #Microbiota y metablytes
                    
                    # Leer el archivo de Excel
                    file_name4 <- "Top_500_Variables_Microbiota_Significant_Correlations_with_Metabolites.xlsx"
                    datavenn4 <- read.xlsx(file_name4, sheet = 1)  # Leer la primera hoja del archivo Excel
                    
                    # Seleccionar las 3 OTUs con más correlaciones significativas
                    top_3_otus <- datavenn4 %>%
                      arrange(desc(Significant_Correlations)) %>%
                      slice(1:3)  # Seleccionar las 3 primeras
                    
                    # Separar las variables relacionadas para cada OTU
                    related_variables_list2 <- top_3_otus %>%
                      mutate(Variables = str_split(Related_Variables, ";")) %>%
                      pull(Variables)  # Extraer la lista de variables relacionadas
                    
                    # Nombres de las OTUs
                    otu_names <- top_3_otus$Var1
                    
                    # Calcular las variables comunes entre las tres OTUs
                    common_variables_all <- Reduce(intersect, related_variables_list2)
                    
                    # Calcular las variables compartidas entre cada par de OTUs
                    pairwise_common_variables <- combn(seq_along(otu_names), 2, function(idx) {
                      otu1 <- otu_names[idx[1]]
                      otu2 <- otu_names[idx[2]]
                      common_vars <- intersect(related_variables_list2[[idx[1]]], related_variables_list2[[idx[2]]])
                      data.frame(OTU1 = otu1, OTU2 = otu2, Variables_Comunes = paste(common_vars, collapse = "; "))
                    }, simplify = FALSE)
                    
                    # Calcular las variables exclusivas para cada OTU
                    exclusive_variables <- lapply(related_variables_list2, function(x) setdiff(x, common_variables_all))
                    
                    # Crear un archivo Excel
                    output_file <- "Variables_Comunes_y_Exclusivas_Microbiota_Metabolitos.xlsx"
                    wb <- createWorkbook()
                    
                    # Agregar hoja para variables comunes entre las tres OTUs
                    addWorksheet(wb, "Comunes_3_OTUs")
                    writeData(wb, "Comunes_3_OTUs", data.frame(Variables_Comunes = common_variables_all))
                    
                    # Agregar hoja para variables compartidas entre pares de OTUs
                    addWorksheet(wb, "Comunes_2_OTUs")
                    pairwise_df <- do.call(rbind, pairwise_common_variables)
                    writeData(wb, "Comunes_2_OTUs", pairwise_df)
                    
                    # Agregar hojas para variables exclusivas por OTU
                    for (i in seq_along(otu_names)) {
                      sheet_name <- paste0("Exclusivas_", otu_names[i])
                      addWorksheet(wb, sheet_name)
                      writeData(wb, sheet_name, data.frame(Variables_Exclusivas = exclusive_variables[[i]]))
                    }
                    
                    # Guardar el archivo Excel
                    saveWorkbook(wb, output_file, overwrite = TRUE)
                    cat("Archivo exportado exitosamente como:", output_file, "\n")
                    
                    # Preparar los datos para el diagrama de Venn
                    venn_data3 <- setNames(related_variables_list2, otu_names)
                    
                    # Crear el diagrama de Venn
                    venn_plot3 <- ggvenn(venn_data3, fill_color = c("#FF9999", "#99CCFF", "#99FF99"))
                    
                    # Visualizar el diagrama en la ventana de gráficos
                    print(venn_plot3)
                    
                    # Guardar el diagrama como archivo PNG
                    output_venn3 <- "Venn_Diagrama_Microbiota_Metabolitos.png"
                    ggsave(output_venn3, venn_plot3, width = 8, height = 6)
                    cat("Diagrama de Venn guardado como:", output_venn3, "\n")
                    

            # Metabolites Y Microbiota 
                    
                    # Leer el archivo de Excel
                    # Leer el archivo de Excel
                    file_name5 <- "Top_50_Variables_Metabolites_Significant_Correlations_with_Microbiota.xlsx"
                    datavenn5 <- read.xlsx(file_name5, sheet = 1)  # Leer la primera hoja del archivo Excel
                    
                    # Seleccionar las METABOLITOS en las filas 2, 3 y 16
                    top_3_Metabolites <- datavenn5[c(24, 15, 19), ]  # Seleccionar las filas 2, 3 y 16
                    
                    # Separar las variables relacionadas para cada METABOLITO
                    related_variables_list3 <- top_3_Metabolites %>%
                      mutate(Variables = str_split(Related_Variables, ";")) %>%
                      pull(Variables)  # Extraer la lista de variables relacionadas
                    
                    # Nombres de las Metabolites
                    metabolyte_names <- top_3_Metabolites$Var1
                    
                    # Calcular las variables comunes entre las tres Metabolites
                    common_variables_all <- Reduce(intersect, related_variables_list3)
                    
                    # Calcular las variables compartidas entre cada par de Metabolites
                    pairwise_common_variables2 <- combn(seq_along(metabolyte_names), 2, function(idx) {
                      Metabolyte1 <- metabolyte_names[idx[1]]
                      Metabolyte2 <- metabolyte_names[idx[2]]
                      common_vars <- intersect(related_variables_list3[[idx[1]]], related_variables_list3[[idx[2]]])
                      data.frame(METABOLYTE1 = Metabolyte1, METABOLYTE2 = Metabolyte2, Variables_Comunes = paste(common_vars, collapse = "; "))
                    }, simplify = FALSE)
                    
                    # Calcular las variables exclusivas para cada Metabolites
                    exclusive_variables <- lapply(related_variables_list3, function(x) setdiff(x, common_variables_all))
                    
                    # Crear un archivo Excel
                    output_file <- "Variables_Comunes_y_Exclusivas_Metabolitos_Microbiota.xlsx"
                    wb <- createWorkbook()
                    
                    # Agregar hoja para variables comunes entre las tres Metabolites
                    addWorksheet(wb, "Comunes_3_Metabolites")
                    writeData(wb, "Comunes_3_Metabolites", data.frame(Variables_Comunes = common_variables_all))
                    
                    # Agregar hoja para variables compartidas entre pares de Metabolites
                    addWorksheet(wb, "Comunes_2_OTUs")
                    pairwise_df <- do.call(rbind, pairwise_common_variables2)
                    writeData(wb, "Comunes_2_OTUs", pairwise_df)
                    
                    # Agregar hojas para variables exclusivas por Metabolites
                    for (i in seq_along(metabolyte_names)) {
                      sheet_name <- paste0("Ex_", metabolyte_names[i])
                      addWorksheet(wb, sheet_name)
                      writeData(wb, sheet_name, data.frame(Variables_Exclusivas = exclusive_variables[[i]]))
                    }
                    
                    # Guardar el archivo Excel
                    saveWorkbook(wb, output_file, overwrite = TRUE)
                    cat("Archivo exportado exitosamente como:", output_file, "\n")
                    
                    # Preparar los datos para el diagrama de Venn
                    venn_data4 <- setNames(related_variables_list3, metabolyte_names)
                    
                    # Crear el diagrama de Venn
                    venn_plot4 <- ggvenn(venn_data4, fill_color = c("#FF9999", "#99CCFF", "#99FF99"))
                    
                    # Visualizar el diagrama en la ventana de gráficos
                    print(venn_plot4)
                    
                    # Guardar el diagrama como archivo PNG
                    output_venn4 <- "Venn_Diagrama_Metabolitos_Microbiota.png"
                    ggsave(output_venn4, venn_plot4, width = 8, height = 6)
                    cat("Diagrama de Venn guardado como:", output_venn4, "\n")
                    
                    
        #PLS-DA (eXPLORAR PARA PERFILAR LOS METABOLOMAS M1 Y M2)
                    
                    library(mixOmics)
                    library(ggplot2)
                    library(dplyr)
                    
                    # Aplicar PLS-DA
                    Metabolites_scaled <- scale(Metabolites[[1]])
                    plsda_model <- plsda(Metabolites_scaled, Group$Two_Metabolome.type, ncomp = 2)  # Seleccionamos 2 componentes
                    
                    dim(Metabolites[[1]])  # Número de filas y columnas de X
                    length(Group$Two_Metabolome.type)
                    
                    # Gráfico de proyección PLS-DA
                    plotIndiv(plsda_model, 
                              comp = c(1,2), 
                              group = Group$Two_Metabolome.type, 
                              legend = TRUE, 
                              title = "PLS-DA de Metabolitos por Tipo de Metaboloma")
                    
                    # Cargar importancia de variables (VIP)
                    vip_scores <- vip(plsda_model)  # Scores de importancia
                    vip_df <- data.frame(Variable = rownames(vip_scores), VIP = vip_scores[,1])
                    
                    # Seleccionar variables con VIP > 1 (más relevantes para la discriminación)
                    vip_top <- vip_df %>% arrange(desc(VIP)) %>% filter(VIP > 1)
                    print(vip_top)
                    
                    
                    # Calcular VIP scores
                    vip_scores <- vip(plsda_model)
                    
                    # Seleccionar variables con VIP > 1
                    top_variables <- names(vip_scores[vip_scores > 1])
                    
                    # Filtrar los loadings1 de estas variables
                    loadings1_top <- plsda_model$loadings$X[vip_top, ]
                    
                    # Graficar solo las variables más discriminantes
                    ggplot() +
                      geom_point(data = scores_df, aes(x = comp1, y = comp2, color = Group), size = 3, alpha = 0.7) +
                      geom_segment(data = as.data.frame(loadings1_top), aes(x = 0, y = 0, xend = comp1, yend = comp2),
                                   arrow = arrow(length = unit(0.3, "cm")), color = "red") +
                      geom_text(data = as.data.frame(loadings1_top), aes(x = comp1, y = comp2, label = rownames(loadings1_top)),
                                vjust = -0.5, hjust = -0.5, color = "red", size = 3) +
                      labs(title = "Biplot de PLS-DA con Variables Discriminantes (VIP > 1)",
                           x = "Componente 1", y = "Componente 2") +
                      theme_minimal()
                    
                    
                    # Ajustar el modelo PLS-DA con 2 componentes
                    plsda_model <- plsda(Metabolites[[1]], Group$Two_Metabolome.type, ncomp = 2)
                    
                    # Verificar estructura del modelo
                    print(summary(plsda_model))
                    
                    # Calcular VIP scores
                    vip_scores <- as.data.frame(vip(plsda_model))
                    
                    # Verificar si los VIP scores son correctos
                    print(summary(vip_scores))
                    
                    # Convertir en data frame con nombres correctos
                    vip_df <- data.frame(Variable = rownames(vip_scores), VIP = vip_scores[, 1])
                    
                    # Ordenar por importancia
                    vip_df <- vip_df %>% arrange(desc(VIP))
                    
                    # Seleccionar las 15 variables más importantes
                    top_variables <- vip_df[1:15, ]
                    
                    # Graficar los VIP Scores de las 15 variables más importantes
                    ggplot(top_variables, aes(x = reorder(Variable, VIP), y = VIP)) +
                      geom_bar(stat = "identity", fill = "steelblue") +
                      coord_flip() +
                      theme_minimal() +
                      labs(title = "Top 15 Variables Discriminantes (VIP Scores)",
                           x = "Variable",
                           y = "VIP Score")
                    
                    
        ######
        ######
                    
                    
                    library(openxlsx)
                    library(dplyr)
                    library(tidyr)
                    library(writexl)
                    
                    # Leer los datos desde el archivo de entrada
                    input_file <- "Summary_correlations_3_Models_Plants_Metabolites.xlsx"
                    data <- read.xlsx(input_file, sheet = "Summary")
                    


                    
                    
                    
          #Funciona          
          #Summary All models correlations plants vs and Metabolites 
                    
                    setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Correlations_Circle_blocks")
                  
                    
                    library(openxlsx)
                    library(dplyr)
                    library(tidyr)
                    library(writexl)
                    
                    # Leer los datos desde el archivo de entrada
                    input_file <- "Summary_correlations_3_Models_Plants_Metabolites.xlsx"
                    data <- read.xlsx(input_file, sheet = "Summary")
                    
                    
                    
                    # Filtrar solo las correlaciones > 0.2
                    data_filtered <- data %>%
                      filter(Correlation > 0.3)
                    
                    # Crear una lista única de plantas
                    plants <- unique(data_filtered$Variable2)
                    
                    # Generar el resumen
                    summary_list <- lapply(plants, function(plant) {
                      # Filtrar los metabolitos asociados a la planta por modelo
                      metab_metabolome <- data_filtered %>%
                        filter(Variable2 == plant, Model == "Metabolome") %>%
                        pull(Variable1) %>%
                        unique() %>%
                        paste(collapse = ";")
                      
                      metab_diet <- data_filtered %>%
                        filter(Variable2 == plant, Model == "Diet") %>%
                        pull(Variable1) %>%
                        unique() %>%
                        paste(collapse = ";")
                      
                      metab_microbiota <- data_filtered %>%
                        filter(Variable2 == plant, Model == "Microbiota") %>%
                        pull(Variable1) %>%
                        unique() %>%
                        paste(collapse = ";")
                      
                      # Identificar los metabolitos comunes a los tres modelos
                      common_metabolites <- intersect(
                        intersect(
                          strsplit(metab_metabolome, ";")[[1]],
                          strsplit(metab_diet, ";")[[1]]
                        ),
                        strsplit(metab_microbiota, ";")[[1]]
                      )
                      
                      model_all <- paste(common_metabolites[common_metabolites != ""], collapse = ";")
                      
                      # Devolver la fila para la planta actual
                      tibble(
                        Plants = plant,
                        Metabolome = metab_metabolome,
                        Diet = metab_diet,
                        Microbiota = metab_microbiota,
                        Model.all = ifelse(model_all == "", "*", model_all)
                      )
                    })
                    
                    # Combinar las filas en un único dataframe
                    summary_df <- bind_rows(summary_list)
                    
                    # Guardar el resultado en un archivo de Excel
                    output_file <- "Summary_all_blocks_Plants_Metabolites.xlsx"
                    write_xlsx(summary_df, output_file)
                    
                    cat("Archivo generado exitosamente: ", output_file, "\n")
                    ")
    }
  ]
}

#no corre, copiar y pegar en otra hoja
library(openxlsx)
library(dplyr)
library(tidyr)
library(writexl)

                    setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/Output/MixOMICS_DIABLO_MODEL/DIABLO_Correlations_Circle_blocks")

                    library(openxlsx)
                    library(dplyr)
                    library(tidyr)
                    library(writexl)
                    
                    # Leer los datos desde el archivo de entrada
                    input_file <- "Summary_correlations_3_Models_OTUs_Metabolites.xlsx"
                    data <- read.xlsx(input_file, sheet = "Summary")
                    
                    # Filtrar solo las correlaciones > 0.2
                    data_filtered <- data %>%
                      filter(Correlation > 0.2)
                    
                    # Crear una lista única de OTUs
                    otus <- unique(data_filtered$Variable2)
                    
                    # Identificar OTUs del modelo Univariate
                    univariate_otus <- data %>%
                      filter(Model == "Univariate") %>%
                      pull(Variable1) %>%
                      unique()
                    
                    # Generar el resumen
                    summary_list <- lapply(otus, function(otu) {
                      # Filtrar los OTUs asociados al metabolito por modelo
                      otus_metabolome <- data_filtered %>%
                        filter(Variable2 == otu, Model == "Metabolome") %>%
                        pull(Variable1) %>%
                        unique() %>%
                        paste(collapse = ";")
                      
                      otus_diet <- data_filtered %>%
                        filter(Variable2 == otu, Model == "Diet") %>%
                        pull(Variable1) %>%
                        unique() %>%
                        paste(collapse = ";")
                      
                      otus_microbiota <- data_filtered %>%
                        filter(Variable2 == otu, Model == "Microbiota") %>%
                        pull(Variable1) %>%
                        unique() %>%
                        paste(collapse = ";")
                      
                      # Identificar los OTUs comunes a los tres modelos
                      common_otus <- intersect(
                        intersect(
                          strsplit(otus_metabolome, ";")[[1]],
                          strsplit(otus_diet, ";")[[1]]
                        ),
                        strsplit(otus_microbiota, ";")[[1]]
                      )
                      
                      model_all <- paste(common_otus[common_otus != ""], collapse = ";")
                      
                      # Identificar OTUs en Univariate
                      all_otus <- unique(c(
                        strsplit(otus_metabolome, ";")[[1]],
                        strsplit(otus_diet, ";")[[1]],
                        strsplit(otus_microbiota, ";")[[1]]
                      ))
                      univariate_hits <- intersect(all_otus, univariate_otus)
                      univariate_column <- paste(univariate_hits, collapse = ";")
                      
                      # Devolver la fila para el OTU actual
                      tibble(
                        OTUs = otu,
                        Metabolome = otus_metabolome,
                        Diet = otus_diet,
                        Microbiota = otus_microbiota,
                        Model.all = ifelse(model_all == "", "*", model_all),
                        Univariate = ifelse(univariate_column == "", "-", univariate_column)
                      )
                    })
                    
                    # Combinar las filas en un único dataframe
                    summary_df <- bind_rows(summary_list)
                    
                    # Guardar el archivo de salida
                    output_file <- "Summary_all_blocks_OTUs_Metabolites.xlsx"
                    write_xlsx(summary_df, output_file)
                    
                    cat("Archivo generado exitosamente con las columnas 'Univariate' y los OTUs comunes en 'Model.all': ", output_file, "\n")")
    }
    ]
    }
    

###riqueza  y diversidad
    
    # Función corregida para calcular riqueza 
    calcular_riqueza <- function(matriz) {
      apply(matriz, 1, function(x) sum(x != 0 & !is.na(x)))  # Cuenta valores distintos de 0 y no NA
    }

    # Aplicar la transformación en tertiles
    
    #Originales
    DB_metabolomics_TENTATIVE_ID_uchuva[c(-1)]    

    data_quartiles <- apply(DB_metabolomics_TENTATIVE_ID_uchuva[c(-1)], 2, function(x) {
      cut(x, 
          breaks = quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 
          labels = c(0, 1, 2, 3), 
          include.lowest = TRUE)
    })
    
    # Convertir a data frame
    data_quartiles <- as.data.frame(data_quartiles)
    data_quartiles1 <- apply(data_quartiles, 2, as.numeric)
    
    
    data_tertiles <- apply(DB_metabolomics_TENTATIVE_ID_uchuva[c(-1)], 2, function(x) {
      cut(x, 
          breaks = quantile(x, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE), 
          labels = c(0, 1, 2), 
          include.lowest = TRUE)
    })
    
    data_tertiles <- as.data.frame(data_tertiles)
    data_tertiles1 <- apply(data_tertiles, 2, as.numeric)
    
    
    data_Mod_25 <- apply(DB_metabolomics_TENTATIVE_ID_uchuva[c(-1)], 2, function(x) {
      p25 <- quantile(x, 0.25, na.rm = TRUE)  # Calcular percentil 25
      x[x < p25] <- 0  # Reemplazar valores menores al percentil 25 con 0
      return(x)
    })
    
    
    # Instalar y cargar paquetes si no los tienes
    install.packages("vegan")
    library(vegan)

    #Plants#veganPlants
    Plants[[1]]
    range(Plants[[1]])
    riqueza_plantas <- calcular_riqueza(Plants[[1]])
    shannon_plantas <- diversity(Plants[[1]], index = "shannon")
    pielou_plantas <- shannon_plantas / log(specnumber(Plants[[1]]))
    Simpson_plantas <- diversity(Plants[[1]], index = "simpson")
    
    
    #Microbiota
    range(Microbiota[[1]])
    riqueza_microbiota <- calcular_riqueza(Microbiota[[1]])
    shannon_Microbiota <- diversity(Microbiota[[1]], index = "shannon")
    pielou_microbiota <- shannon_Microbiota / log(specnumber(Microbiota[[1]]))
    Simpson_Microbiota <- diversity(Microbiota[[1]], index = "simpson")
    
    
    #Metabolitos
    range(data_quartiles1)
    riqueza_metabolitos <- calcular_riqueza(data_Mod_25)  # Quita primera columna si es ID
    shannon_metabolitos <- diversity(DB_metabolomics_TENTATIVE_ID_uchuva[c(-1)], index = "shannon")
    pielou_metabolitos <- shannon_metabolitos / log(specnumber(DB_metabolomics_TENTATIVE_ID_uchuva[c(-1)]))
    Simpson_metabolitos <- diversity(DB_metabolomics_TENTATIVE_ID_uchuva[c(-1)], index = "simpson")
    


    # DIversidad y riqueza integrados
    datos_Diversity <- data.frame(Simpson_plantas, Simpson_Microbiota,Simpson_metabolitos, riqueza_plantas, riqueza_microbiota,riqueza_metabolitos,  shannon_Microbiota,pielou_microbiota, shannon_metabolitos, pielou_metabolitos, shannon_plantas,pielou_plantas)
    
    
    
    cor.test(datos_Diversity$shannon_Microbiota, datos_Diversity$shannon_metabolitos, method = "spearman")
    cor.test(datos_Diversity$shannon_plantas, datos_Diversity$shannon_metabolitos, method = "spearman") #mejor p<0.13
    cor.test(datos_Diversity$shannon_plantas, datos_Diversity$shannon_Microbiota, method = "spearman") #mejor p<0.18  #mejó a p<0.18  con CLR2
    
    cor.test(datos_Diversity$pielou_plantas, datos_Diversity$pielou_microbiota, method = "spearman")
    cor.test(datos_Diversity$pielou_plantas, datos_Diversity$pielou_metabolitos, method = "spearman") #mejor con CLR2p<0.12
    cor.test(datos_Diversity$pielou_microbiota, datos_Diversity$pielou_metabolitos, method = "spearman")
    
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$pielou_microbiota, method = "spearman") #mejor con CLR2p<0.02
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$pielou_metabolitos, method = "spearman") #mejor p<0.12 mejoró con metaboloma tertilado 0.08347
    cor.test(datos_Diversity$riqueza_microbiota, datos_Diversity$pielou_metabolitos, method = "spearman")
    
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$shannon_metabolitos, method = "spearman") #mejor p<0.12
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$pielou_microbiota, method = "spearman") #mejor con CLR2p<0.011
    cor.test(datos_Diversity$riqueza_microbiota, datos_Diversity$shannon_metabolitos, method = "spearman")
    
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$riqueza_microbiota, method = "spearman") #mejor con CLR2 p<0.001
    cor.test(datos_Diversity$riqueza_microbiota, datos_Diversity$shannon_plantas, method = "spearman") #mejor con CLR2p<0.005

    cor.test(datos_Diversity$Simpson_Microbiota, datos_Diversity$Simpson_metabolitos, method = "spearman")
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$Simpson_metabolitos, method = "spearman") #mejor p<0.07
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$Simpson_Microbiota, method = "spearman")#mejoró con CLR2p<0.03
    
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$shannon_metabolitos, method = "spearman") 
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$shannon_Microbiota, method = "spearman")  #mejoró con CLR2p<0.02
    
    cor.test(datos_Diversity$shannon_plantas, datos_Diversity$Simpson_metabolitos, method = "spearman") #mejor p<0.05
    cor.test(datos_Diversity$shannon_plantas, datos_Diversity$Simpson_Microbiota, method = "spearman")  #mejoró con CLR2p<0.009
    cor.test(datos_Diversity$shannon_Microbiota, datos_Diversity$Simpson_metabolitos, method = "spearman") 
    
    
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$riqueza_metabolitos, method = "spearman") 
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$riqueza_microbiota, method = "spearman") #mejoró con CLR2p<0.02
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$riqueza_metabolitos, method = "spearman") 
    cor.test(datos_Diversity$riqueza_microbiota, datos_Diversity$riqueza_metabolitos, method = "spearman") 
    
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$Simpson_metabolitos, method = "spearman") #mejor p<0.04
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$Simpson_Microbiota, method = "spearman")  #mejor p<0.16  #mejoró con CLR2p<0.002
    cor.test(datos_Diversity$riqueza_microbiota, datos_Diversity$Simpson_metabolitos, method = "spearman") 
    
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$pielou_metabolitos, method = "spearman")  #Mejoro con metabolitos tertiladas
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$pielou_microbiota, method = "spearman") #mejor CLR2 p<0.05
    
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$Simpson_Microbiota, method = "spearman")  #Mejoro con metabolitos tertiladas
    cor.test(datos_Diversity$Simpson_plantas, datos_Diversity$Simpson_Microbiota, method = "spearman")  #Mejoro con metabolitos tertiladas
    
    
    cor.test(datos_Diversity$pielou_plantas, datos_Diversity$Simpson_metabolitos, method = "spearman")  #mejor CLR2 p<0.12
    cor.test(datos_Diversity$pielou_plantas, datos_Diversity$Simpson_Microbiota, method = "spearman") 
    cor.test(datos_Diversity$pielou_microbiota, datos_Diversity$Simpson_metabolitos, method = "spearman") 
    
    
    
    modelo_glm1 <- glm(shannon_plantas ~ shannon_Microbiota, data = datos_Diversity, family = gaussian())
    modelo_glm2 <- glm(shannon_plantas ~ shannon_metabolitos, data = datos_Diversity, family = gaussian())
    modelo_glm3 <- glm(shannon_metabolitos ~ shannon_Microbiota, data = datos_Diversity, family = gaussian())
    
    modelo_glm4 <- glm(pielou_plantas ~ pielou_microbiota, data = datos_Diversity, family = gaussian())
    modelo_glm5 <- glm(pielou_plantas ~ pielou_metabolitos, data = datos_Diversity, family = gaussian())
    modelo_glm6 <- glm(pielou_metabolitos ~ pielou_microbiota, data = datos_Diversity, family = gaussian())
    
    modelo_glm7 <- glm(riqueza_plantas ~ shannon_Microbiota, data = datos_Diversity, family = gaussian())
    modelo_glm8 <- glm(riqueza_plantas ~ shannon_metabolitos, data = datos_Diversity, family = gaussian())
    modelo_glm9 <- glm(riqueza_microbiota ~ shannon_metabolitos, data = datos_Diversity, family = gaussian())
    
    modelo_glm10 <- glm(riqueza_plantas ~ pielou_microbiota, data = datos_Diversity, family = gaussian())
    modelo_glm11 <- glm(riqueza_plantas ~ pielou_metabolitos, data = datos_Diversity, family = gaussian())
    modelo_glm12 <- glm(riqueza_microbiota ~ pielou_metabolitos, data = datos_Diversity, family = gaussian())
    
    modelo_glm13 <-glm(Simpson_plantas ~ Simpson_metabolitos, data = datos_Diversity, family = gaussian())
      
    modelo_glm14 <- glm(shannon_plantas ~ Simpson_metabolitos, data = datos_Diversity, family = gaussian())
    
    modelo_glm15 <- glm(riqueza_plantas ~ Simpson_metabolitos, data = datos_Diversity, family = gaussian())

    
    
    
    summary(modelo_glm1) #mejor p<0.14 platas vs microbiota
    summary(modelo_glm2)
    summary(modelo_glm3)
    summary(modelo_glm4)
    summary(modelo_glm5)
    summary(modelo_glm6)
    summary(modelo_glm7) #mejor p<0.15 platas vs microbiota
    summary(modelo_glm8)
    summary(modelo_glm9)
    summary(modelo_glm10)
    summary(modelo_glm11)
    summary(modelo_glm12)
    summary(modelo_glm13)
    summary(modelo_glm14)
    summary(modelo_glm15)
    
    
    summary(datos_RICH$shannon_plantas)
    summary(datos_RICH$shannon_Microbiota)
    
    library(ggplot2)
    datos_Diversity$Muestra_ID <- rownames(datos_Diversity)
    
    ggplot(datos_Diversity, aes(x = riqueza_plantas , y = riqueza_microbiota,label = Muestra_ID)) +
      geom_point() +
      geom_text(vjust = -0.5, size = 3) +  # Agrega los nombres de las observaciones
      geom_smooth(method = "lm", se = FALSE) +
      theme_minimal()
    
  cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$riqueza_microbiota, method = "spearman") #mejor con CLR2 p<0.001

    
    datos_Diversity$Muestra_ID <- rownames(datos_Diversity)
    
    ggplot(datos_Diversity, aes(x = riqueza_plantas[-c(2,10)], y = shannon_metabolitos[-c(2,10)], label = Muestra_ID[-c(2,10)])) +
      geom_point() +
      geom_text(vjust = -0.5, size = 3) +  # Agrega los nombres de las observaciones
      geom_smooth(method = "lm", se = FALSE) +
      theme_minimal()
    
    
    cor.test(datos_Diversity$riqueza_plantas, datos_Diversity$Simpson_metabolitos, method = "spearman") #mejor p<0.04
    
    
    
    #figura metaboloma clasificacion
    
    # Instalar y cargar paquetes necesarios
    library(ggplot2)
    library(openxlsx)
    setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/RAW")
    
    # Leer el archivo Excel (reemplaza con la ruta correcta)
    
    MetCla <- read.xlsx("Tabla_Metabolomics_results.xlsx", sheet = 3)  # Puedes cambiar el número de hoja    

    # Contar la cantidad de compuestos por categoría
    MetCla_summary <- as.data.frame(table(MetCla$class_hmdb))
    colnames(MetCla_summary) <- c("class_hmdb", "count")
    
    library(RColorBrewer)
    num_colors <- length(unique(MetCla_summary$class_hmdb))  # Número de categorías
    palette_colors <- colorRampPalette(brewer.pal(8, "Set2"))(num_colors)  # Generar colores extra
    
    ggplot(MetCla_summary, aes(x = reorder(class_hmdb, count), y = count, fill = class_hmdb)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      coord_flip() +
      scale_fill_manual(values = palette_colors) +  # Aplicar la paleta extendida
      theme_minimal(base_size = 14) +
      labs(title = "Distribución de Class_HMDB en los Compounds",
           x = "Class_HMDB", y = "Cantidad de Compounds") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank())


    library(ggplot2)
    library(dplyr)
    library(viridis)  # Para colores atractivos
    
    # Paleta de colores
    num_colors <- length(unique(MetCla_summary$class_hmdb))
    palette_colors <- colorRampPalette(viridis::viridis(num_colors))(num_colors)
    
    # Gráfico con eje X segmentado
    ggplot(MetCla_summary, aes(y = reorder(class_hmdb, count), x = count, fill = class_hmdb)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      scale_fill_manual(values = palette_colors) +
      scale_x_continuous(
        breaks = c(0, 30, 60, 90, 120, 160),  # Puntos de referencia en el eje
        labels = c("0", "30", "60", "90", "120", "160"),  # Etiquetas visibles
        limits = c(0, 160)  # Rango completo
      ) +
      theme_minimal(base_size = 14) +
      labs(
        title = "Distribución de Class_HMDB en los Compounds",
        x = "Cantidad de Compounds",
        y = "Class_HMDB"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        panel.grid.major.x = element_line(color = "gray80"),  # Mantener líneas guía
        panel.grid.minor = element_blank()
      )
    
    library(ggbreak)
    library(yulab.utils)
    library(aplot)
    library(ggfun)
    library(ggplotify)
    library(gridGraphics)
    
    # Contar la cantidad de compuestos por categoría
    MetCla_summary <- as.data.frame(table(MetCla$super_class_hmdb))
    colnames(MetCla_summary) <- c("super_class_hmdb", "count")
    
    
    
    # Paleta de colores personalizada
    num_colors <- length(unique(MetCla_summary$class_hmdb))
    palette_colors <- colorRampPalette(viridis::viridis(num_colors))(num_colors)
    
    # Gráfico con eje X segmentado e invertido
    ggplot(MetCla_summary, aes(y = reorder(class_hmdb, -count), x = count, fill = class_hmdb)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      scale_fill_manual(values = palette_colors) +
      scale_x_break(c(20, 100), scales = c(0.8, 0.2)) +  # 📌 Segmentación 80% - 20%
      coord_flip() +  # 🔄 Invierte el gráfico
      theme_minimal(base_size = 14) +
      labs(
        title = "",
        x = "Number of Metabolites",
        y = "Class by HMDB"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(size = 12, hjust = 1 ),  # 🔥 Texto vertical
        axis.text.x = element_text(size = 12,angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    
    
    library(ggplot2)
    library(ggbreak)
    library(viridis)
    
    # Crear tabla de frecuencias
    MetCla <- read.xlsx("Tabla_Metabolomics_results.xlsx", sheet = 3)  # Puedes cambiar el número de hoja    
    
    MetCla_summary <- as.data.frame(table(MetCla$super_class_hmdb))
    colnames(MetCla_summary) <- c("super_class_hmdb", "count")
    
    # Asegurar que la paleta tiene suficientes colores
    num_colors <- length(unique(MetCla_summary$super_class_hmdb))
    palette_colors <- viridis(num_colors)  # Genera exactamente la cantidad de colores necesaria
    
    # Gráfico con segmentación y orden invertido
    ggplot(MetCla_summary, aes(y = reorder(super_class_hmdb, -count), x = count, fill = super_class_hmdb)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      scale_fill_manual(values = palette_colors) +  # Usa la paleta corregida
      scale_x_break(c(40, 100), scales = c(0.8, 0.2)) +  # 📌 Segmentación 80% - 20%
      coord_flip() +  # 🔄 Invierte el gráfico
      theme_minimal(base_size = 14) +
      labs(
        title = "",
        x = "Number of Metabolites",
        y = "Super class by HMDB"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(size = 12, hjust = 1),  # 🔥 Mejor alineación
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    
    library(ggplot2)
    library(dplyr)
    
    library(ggplot2)
    library(dplyr)
    
    # Crear tabla de frecuencias
    MetCla_summary <- as.data.frame(table(MetCla$super_class_hmdb))
    colnames(MetCla_summary) <- c("super_class_hmdb", "count")
    
    # Ordenar los datos por cantidad
    MetCla_summary <- MetCla_summary %>%
      arrange(desc(count)) %>%
      mutate(label_pos = count + max(count) * 0.05)  # Posición de etiquetas
    
    # Crear gráfico
    ggplot(MetCla_summary, aes(x = reorder(super_class_hmdb, -count), y = count)) +
      geom_segment(aes(xend = super_class_hmdb, yend = 0), color = "black", size = 1) +  # Líneas en negro
      geom_point(aes(y = count), size = 9, color = "black", fill = "#ffa500", shape = 21) +  # Círculos gris claro
      geom_text(aes(y = count, label = count), color = "black", size = 4, fontface = "bold") +  # Números en negro
      coord_flip() +  # Girar el gráfico
      theme_minimal(base_size = 14) +
      labs(
        title = "",
        x = "Super class by HMDB",
        y = "Number of Metabolites"
      ) +
      theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    
    
    #### Metabolomics classificatioin
    
    #figura metaboloma clasificacion
    
    # Instalar y cargar paquetes necesarios
    
    library(ggplot2)
    library(openxlsx)
    library(RColorBrewer)
    
    setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/RAW")
    
    # Leer el archivo Excel (reemplaza con la ruta correcta)
    
    MetCla <- read.xlsx("H_Metabolomics_Taxonomy_&_Classification.xlsx", sheet = 1)  # Puedes cambiar el número de hoja    
    
    # Contar la cantidad de compuestos por categoría
    MetCla_summary <- as.data.frame(table(MetCla$classhmdb))
    colnames(MetCla_summary) <- c("classhmdb", "count")
    
    num_colors <- length(unique(MetCla_summary$classhmdb))  # Número de categorías
    palette_colors <- colorRampPalette(brewer.pal(8, "Set2"))(num_colors)  # Generar colores extra
    
    ggplot(MetCla_summary, aes(x = reorder(classhmdb, count), y = count, fill = classhmdb)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      coord_flip() +
      scale_fill_manual(values = palette_colors) +  # Aplicar la paleta extendida
      theme_minimal(base_size = 14) +
      labs(title = "Distribución de classhmdb en los Compounds",
           x = "classhmdb", y = "Cantidad de Compounds") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank())
    
    
    library(ggplot2)
    library(dplyr)
    library(viridis)  # Para colores atractivos
    
    # Paleta de colores
    num_colors <- length(unique(MetCla_summary$classhmdb))
    palette_colors <- colorRampPalette(viridis::viridis(num_colors))(num_colors)
    
    # Gráfico con eje X segmentado
    ggplot(MetCla_summary, aes(y = reorder(classhmdb, count), x = count, fill = classhmdb)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      scale_fill_manual(values = palette_colors) +
      scale_x_continuous(
        breaks = c(0, 30, 60, 90, 120, 150),  # Puntos de referencia en el eje
        labels = c("0", "30", "60", "90", "120", "150"),  # Etiquetas visibles
        limits = c(0, 150)  # Rango completo
      ) +
      theme_minimal(base_size = 14) +
      labs(
        title = "Distribución de classhmdb en los Compounds",
        x = "Cantidad de Compounds",
        y = "classhmdb"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        panel.grid.major.x = element_line(color = "gray80"),  # Mantener líneas guía
        panel.grid.minor = element_blank()
      )
    
    
    
    
    library(ggbreak)
    library(yulab.utils)
    library(aplot)
    library(ggfun)
    library(ggplotify)
    library(gridGraphics)
    library(ggplot2)
    library(openxlsx)
    library(RColorBrewer)
    library(viridis) 
    library(ggplot2)
    library(dplyr)
    
    
    
    setwd("C:/INFORMACION D/R/MixOmics/Metadiet/Paper_Metadiet/RAW")
    
    # Leer el archivo Excel (reemplaza con la ruta correcta)
    
    MetCla <- read.xlsx("H_Metabolomics_Taxonomy_&_Classification.xlsx", sheet = 1)  # Puedes cambiar el número de hoja    
    
    # Contar la cantidad de compuestos por categoría
    MetCla_summary <- as.data.frame(table(MetCla$classhmdb))
    colnames(MetCla_summary) <- c("classhmdb", "count")
    
    
    
    # Paleta de colores personalizada
    num_colors <- length(unique(MetCla_summary$classhmdb))
    palette_colors <- colorRampPalette(viridis::viridis(num_colors))(num_colors)
    
    # Gráfico con eje X segmentado e invertido
    ggplot(MetCla_summary, aes(y = reorder(classhmdb, -count), x = count, fill = classhmdb)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      scale_fill_manual(values = palette_colors) +
      scale_x_break(c(30, 70), scales = c(0.8, 0.2)) +  # 📌 Segmentación 80% - 20%
      coord_flip() +  # 🔄 Invierte el gráfico
      theme_minimal(base_size = 14) +
      labs(
        title = "",
        x = "Number of Metabolites",
        y = "Class by HMDB"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(size = 12, hjust = 1 ),  # 🔥 Texto vertical
        axis.text.x = element_text(size = 12,angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    
    
    library(ggplot2)
    library(ggbreak)
    library(viridis)
    
    # Crear tabla de frecuencias
    
    MetCla1 <- read.xlsx("H_Metabolomics_Taxonomy_&_Classification.xlsx", sheet = 1)  # Puedes cambiar el número de hoja    
    MetCla2 <- read.xlsx("H_Metabolomics_Taxonomy_&_Classification.xlsx", sheet = 2)  # Puedes cambiar el número de hoja    
    
    MetCla_summary_S1 <- as.data.frame(table(MetCla1$superclasshmdb))
    MetCla_summary_S2 <- as.data.frame(table(MetCla2$Tentative.superclass.classification.no.HMDB))
    
    library(dplyr)
    
    # Filtrar los "Not available" solo de S2
    na_S2 <- MetCla_summary_S2 %>% filter(Var1 == "Not available")
    
    # Filtrar los datos útiles (excluyendo los "Not available" de S1)
    MetCla_S1_filtered <- MetCla_summary_S1 %>% filter(Var1 != "Not available")
    
    # Filtrar los datos útiles de S2 (sin tocar los "Not available")
    MetCla_S2_filtered <- MetCla_summary_S2 %>% filter(Var1 != "Not available")
    
    # Combinar todo: S1 (sin NA), S2 (sin NA), y los "Not available" una sola vez (de S2)
    MetCla_combined <- bind_rows(MetCla_S1_filtered, MetCla_S2_filtered, na_S2)
    
    # Agrupar y sumar frecuencias
    MetCla_summary_total <- MetCla_combined %>%
      group_by(Var1) %>%
      summarise(Freq = sum(Freq)) %>%
      arrange(desc(Freq))
    
    # Ver resultado
    print(MetCla_summary_total)
    sum(MetCla_summary_total$Freq)
    
    
    colnames(MetCla_summary_total) <- c("superclasshmdb", "count")
    
    # Asegurar que la paleta tiene suficientes colores
    num_colors <- length(unique(MetCla_summary_total$superclasshmdb))
    palette_colors <- viridis(num_colors)  # Genera exactamente la cantidad de colores necesaria
    
    # Gráfico con segmentación y orden invertido
    ggplot(MetCla_summary_total, aes(y = reorder(superclasshmdb, -count), x = count, fill = superclasshmdb)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      scale_fill_manual(values = palette_colors) +  # Usa la paleta corregida
      scale_x_break(c(55, 100), scales = c(0.8, 0.2)) +  # 📌 Segmentación 80% - 20%
      coord_flip() +  # 🔄 Invierte el gráfico
      theme_minimal(base_size = 14) +
      labs(
        title = "",
        x = "Number of Metabolites",
        y = "Super class by HMDB"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(size = 12, hjust = 1),  # 🔥 Mejor alineación
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    
    
    library(ggplot2)
    library(dplyr)
    
    # Crear tabla de frecuencias
    colnames(MetCla_summary_total) <- c("superclasshmdb", "count")
    
    # Ordenar los datos por cantidad
    MetCla_summary_total <- MetCla_summary_total %>%
      arrange(desc(count)) %>%
      mutate(label_pos = count + max(count) * 0.05)  # Posición de etiquetas
    
    # Crear gráfico
    ggplot(MetCla_summary_total, aes(x = reorder(superclasshmdb, -count), y = count)) +
      geom_segment(aes(xend = superclasshmdb, yend = 0), color = "black", size = 1) +  # Líneas en negro
      geom_point(aes(y = count), size = 9, color = "black", fill = "#ffa500", shape = 21) +  # Círculos gris claro
      geom_text(aes(y = count, label = count), color = "black", size = 4, fontface = "bold") +  # Números en negro
      coord_flip() +  # Girar el gráfico
      theme_minimal(base_size = 14) +
      labs(
        title = "",
        x = "Super class by HMDB",
        y = "Number of Metabolites"
      ) +
      theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    
    #Canonical Correlation Analysis (CCA)
    
    PlantsCCA <- data.frame (B_OTU_Plants_CLR[c(-1,-2)])
    MetabolitesCCA <- data.frame (scale (G_Metabolome_IsoMS[c(-1,-2)])) 
    MicrobiotaCCA <- data.frame (E_OTU_Microbiota_CLR[c(-1,-2)])
    
    
    library(CCA)
    library(fda)
    library(fds)
    library(rainbow)
    library(ks)
    library(hdrcde)
    library(RCurl)
    library(deSolve)
    library(fields)
    library(Rcpp)
    
    library(spam)
    
    
    
    
    
    library(pcaPP)
    
    
    library(ggplot2)
    
    # Asegúrate de que las matrices tengan el mismo número de filas
    dim(PlantsCCA)
    dim(MetabolitesCCA)
    
    # CCA entre ingesta de plantas y metabolitos
    cca_plants_metabolites <- cc(PlantsCCA, MetabolitesCCA)
    
    # Ver las correlaciones canónicas
    cca_plants_metabolites$cor
    