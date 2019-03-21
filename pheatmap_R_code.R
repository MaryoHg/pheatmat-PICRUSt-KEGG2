#######################################
## HEATMAP WITH KEGG2 FUNCTIONALITY
#######################################

## 1. Packages and functions
options(scipen = 100000000) ## avoid scientific annotation (e+) values.
library(pheatmap)
source('/home/mario/Dropbox/R_functions/summarySE.R')
source('/home/mario/Dropbox/R_functions/tranpose.any.table.R')

## 2. dataframe: complete data frame with KEGG2 values.
K2 <- read.delim('/home/mario/Documents/Doctorado/PROYECTOS/Maryo_Hg/ammonia/downstream/picrust_3000/KEGG2_LONG_FINAL_DATA.txt')
colnames(K2)


## 3. Sorting data to better visualize the tratments.
        
        ##Sort data with dynamic, days, and dynTreat:
        K2.sorted <- K2[order(K2$dynamic, K2$days, K2$dynTREATday),]
        View(K2.sorted)
        

## 4. AVERAGE, SD and SE for all treatments:
K2.sorted.average <- summarySE(data = K2.sorted,  measurevar = "log10x",
                               groupvars = c("dynamic", "Level_2", "Treatment", "days"),
                               na.rm = T)
View(K2.sorted.average)
        
        ## a. create a unique label:
        K2.sorted.average$label <- factor(paste (K2.sorted.average$dynamic, K2.sorted.average$Treatment, K2.sorted.average$days, sep = '.'))
        View(K2.sorted.average); 
        length(levels(K2.sorted.average$label)) # 96 valores únicos.
        
        ## b. adding LABEL (unique and ordered label to heatmap) with CONDITIONALS:
        #indata <- c("am.CT0.0", "am.CT0.1", "am.CT0.3", "am.CT0.5", "am.CT0.7", "am.CT0.14", "am.CT0.28", "am.CT0.56", "am.CT3.0", "am.CT3.1", "am.CT3.3", "am.CT3.5", "am.CT3.7", "am.CT3.14", "am.CT3.28", "am.CT3.56", "am.PBB0.0", "am.PBB0.1", "am.PBB0.3", "am.PBB0.5", "am.PBB0.7", "am.PBB0.14", "am.PBB0.28", "am.PBB0.56", "am.PBB3.0", "am.PBB3.1", "am.PBB3.3", "am.PBB3.5", "am.PBB3.7", "am.PBB3.14", "am.PBB3.28", "am.PBB3.56", "am.PBC0.0", "am.PBC0.1", "am.PBC0.3", "am.PBC0.5", "am.PBC0.7", "am.PBC0.14", "am.PBC0.28", "am.PBC0.56", "am.PBC3.0", "am.PBC3.1", "am.PBC3.3", "am.PBC3.5", "am.PBC3.7", "am.PBC3.14", "am.PBC3.28", "am.PBC3.56", "un.CT0.0", "un.CT0.1", "un.CT0.3", "un.CT0.5", "un.CT0.7", "un.CT0.14", "un.CT0.28", "un.CT0.56", "un.CT3.0", "un.CT3.1", "un.CT3.3", "un.CT3.5", "un.CT3.7", "un.CT3.14", "un.CT3.28", "un.CT3.56", "un.PBB0.0", "un.PBB0.1", "un.PBB0.3", "un.PBB0.5", "un.PBB0.7", "un.PBB0.14", "un.PBB0.28", "un.PBB0.56", "un.PBB3.0", "un.PBB3.1", "un.PBB3.3", "un.PBB3.5", "un.PBB3.7", "un.PBB3.14", "un.PBB3.28", "un.PBB3.56", "un.PBC0.0", "un.PBC0.1", "un.PBC0.3", "un.PBC0.5", "un.PBC0.7", "un.PBC0.14", "un.PBC0.28", "un.PBC0.56", "un.PBC3.0", "un.PBC3.1", "un.PBC3.3", "un.PBC3.5", "un.PBC3.7", "un.PBC3.14", "un.PBC3.28", "un.PBC3.56")
        #newlabel <- c("A1", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "C1", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "E1", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "G1", "G11", "G12", "G13", "G14", "G15", "G16", "G17", "I1", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "K1", "K11", "K12", "K13", "K14", "K15", "K16", "K17", "B1", "B11", "B12", "B13", "B14", "B15", "B16", "B17", "D1", "D11", "D12", "D13", "D14", "D15", "D16", "D17", "F1", "F11", "F12", "F13", "F14", "F15", "F16", "F17", "H1", "H11", "H12", "H13", "H14", "H15", "H16", "H17", "J1", "J11", "J12", "J13", "J14", "J15", "J16", "J17", "L1", "L11", "L12", "L13", "L14", "L15", "L16", "L17")
        indata <- c("am.CT0.0", "am.CT0.1", "am.CT0.3", "am.CT0.5", "am.CT0.7", "am.CT0.14", "am.CT0.28", "am.CT0.56", "un.CT0.0", "un.CT0.1", "un.CT0.3", "un.CT0.5", "un.CT0.7", "un.CT0.14", "un.CT0.28", "un.CT0.56", "am.CT3.0", "am.CT3.1", "am.CT3.3", "am.CT3.5", "am.CT3.7", "am.CT3.14", "am.CT3.28", "am.CT3.56", "un.CT3.0", "un.CT3.1", "un.CT3.3", "un.CT3.5", "un.CT3.7", "un.CT3.14", "un.CT3.28", "un.CT3.56", "am.PBB0.0", "am.PBB0.1", "am.PBB0.3", "am.PBB0.5", "am.PBB0.7", "am.PBB0.14", "am.PBB0.28", "am.PBB0.56", "un.PBB0.0", "un.PBB0.1", "un.PBB0.3", "un.PBB0.5", "un.PBB0.7", "un.PBB0.14", "un.PBB0.28", "un.PBB0.56", "am.PBB3.0", "am.PBB3.1", "am.PBB3.3", "am.PBB3.5", "am.PBB3.7", "am.PBB3.14", "am.PBB3.28", "am.PBB3.56", "un.PBB3.0", "un.PBB3.1", "un.PBB3.3", "un.PBB3.5", "un.PBB3.7", "un.PBB3.14", "un.PBB3.28", "un.PBB3.56", "am.PBC0.0", "am.PBC0.1", "am.PBC0.3", "am.PBC0.5", "am.PBC0.7", "am.PBC0.14", "am.PBC0.28", "am.PBC0.56", "un.PBC0.0", "un.PBC0.1", "un.PBC0.3", "un.PBC0.5", "un.PBC0.7", "un.PBC0.14", "un.PBC0.28", "un.PBC0.56", "am.PBC3.0", "am.PBC3.1", "am.PBC3.3", "am.PBC3.5", "am.PBC3.7", "am.PBC3.14", "am.PBC3.28", "am.PBC3.56", "un.PBC3.0", "un.PBC3.1", "un.PBC3.3", "un.PBC3.5", "un.PBC3.7", "un.PBC3.14", "un.PBC3.28", "un.PBC3.56")
        
        newlabel <- c("A1", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "B1", "B11", "B12", "B13", "B14", "B15", "B16", "B17", "C1", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "D1", "D11", "D12", "D13", "D14", "D15", "D16", "D17", "E1", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "F1", "F11", "F12", "F13", "F14", "F15", "F16", "F17", "G1", "G11", "G12", "G13", "G14", "G15", "G16", "G17", "H1", "H11", "H12", "H13", "H14", "H15", "H16", "H17", "I1", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "J1", "J11", "J12", "J13", "J14", "J15", "J16", "J17", "K1", "K11", "K12", "K13", "K14", "K15", "K16", "K17", "L1", "L11", "L12", "L13", "L14", "L15", "L16", "L17")
        
        
        K2.sorted.average$LABEL <- newlabel [match(K2.sorted.average$label, indata)]
        View(K2.sorted.average)
        
        
## Long to wide data to be able to use pheatmap:
        library (tidyr)
        colnames(K2.sorted.average)
        data_mat <- spread (data = K2.sorted.average[,c(2,6,11)],
                           key = "LABEL",
                           value = "log10x")
        
        data_mat <- data_mat[-4,] ## QUITAMOS 'Cell Communication'.
        rownames(data_mat) <- data_mat[,1]
        View(data_mat)
        
        ## Is it necessary ordered by rowSum? ¡NO!
        #data_mat_order <- data_mat[order(-rowSums(data_mat[,2:ncol(data_mat)])),]
        #View(data_mat_order)
        #rownames(data_mat_order) <- data_mat[,1]
        
        ## Heatmap with annotation: loading annotation mapping file
        annotation <- read.csv('/home/mario/Documents/Doctorado/PROYECTOS/Maryo_Hg/ammonia/downstream/picrust_3000/annotation.csv',
                               header = T, sep = ',')
        colnames(annotation)
        rownames(annotation) <- annotation$LABEL
        annotation <- annotation[c(6, 4, 5)]

        ## Specify colors editing annotation ¡VARIABLES!:
        dyn <- c('black', 'white')
        names (dyn) <- c('200 mg N kg-1','0 mg N kg-1')
        dyn
        
        #days = colorRampPalette(c("white", "yellow"))(8)
        #days <- c('#000000', '#696969', '#808080', '#A9A9A9', '#C0C0C0', '#D3D3D3', '#DCDCDC', '#F5F5F5')
        #days <- c('#ffffff', '#ffcccc', '#ff9999', '#ff6666', '#ff3333', '#ff0000', '#cc0000', '#800000')
        days <- c('#E6E6FA', '#DDA0DD', '#EE82EE', '#BA55D3', '#8A2BE2', '#9932CC', '#800080', '#4B0082')
        names(days) = c(0,1,3,5,7,14,28,56)
        days
        
        #Treats = c('red', 'blue', 'yellow', 'black', "lightgreen", "navy")
        Treats <- c('#C0C0C0', '#FF0000', '#FFFF00', '#00FF00', '#00FFFF', '#FF00FF')
        names (Treats) <- c("CT0", 'CT3', 'PBB0', 'PBB3', 'PBC0', 'PBC3')
        Treats
        
        ann_colors = list(treatment = Treats, ammonia=dyn, days = days)
        ann_colors
        
                ## This heatmap uses the default color palette for pheatmap()
        pheatmap(data_mat[,2:ncol(data_mat)],
                 annotation = annotation, annotation_colors = ann_colors,
                 cluster_cols = F, cluster_rows = T, fontsize = 8,
                 cellwidth = 6.35, cellheight = 16, facefont='italic',
                 show_colnames = F, cutree_rows = 3,
                 filename = 'KEGG2_HEATMAP_FINAL.png') #to save it as pdf, png, etc.)
                # options to control the color range for heatmap value:
                #legend_breaks=0:6,legend_labels=c(0,1,2,3,4,5,6))
                #
        dev.off()

                ## Editing the color palette for the same heatmap.
        pheatmap(data_mat[,2:ncol(data_mat)], 
                 annotation = annotation, annotation_colors = ann_colors,
                 cluster_cols = F, cluster_rows = T, fontsize = 8,
                 cellwidth = 6.35, cellheight = 16, facefont='italic',
                 show_colnames = F,
                 cutree_rows = 3, color = colorRampPalette(c("navy", 'white', "firebrick3"))(50),
                 filename = 'KEGG2_HEATMAP2.png')
        dev.off()
        
        
        #my_palette <- colorRampPalette(c("white","yellow","red"))(n=599)
       
        ## This is, too, a beautiful way of colors.
        #png('/home/mario/Documents/Doctorado/PROYECTOS/Maryo_Hg/ammonia/downstream/picrust_3000/heatmap.png',
        #    width = 1390, height = 720, units = 'px')
        #pheatmap(data_mat[,2:ncol(data_mat)],
        #         cluster_rows = T, cluster_cols = F,
        #         fontsize = 8, fontface = 'italic', cellwidth = 6.35, cellheight = 16
        #         color=colorRampPalette(c("navy", "white", "firebrick3"))(50))
        #dev.off()
        

#########################
#########################
        ## this does not work anymore.
## 5. Creating Factor column to heatmap:
        ## divide data into deciles: 10 groups
        library (dplyr)
        data <- K2.sorted.average
        data$DECILES <- ntile (data$log10x, 10)
        
        ## Create breaks for colors:
        data$FACTOR <- cut (data$log10x,
                            breaks = c(-1,0,1,2,3,4,5,max(data$log10x)),
                            labels = c('0', '0 - 1.0', '1.01 - 2.0', '2.01 - 3.0', '3.01 - 4.0',
                                                  '4.01 - 5.0', '5.01 - 6.0'))
        ## Change level order: abundanceFactor:
        data$FACTOR <- factor (as.character(data$FACTOR),
                               levels = rev (levels(data$FACTOR)))
        colnames (data)
        
        ## GRAPH IT:
        data.am <- subset(data, dynamic=='am')
        indata <- c("am.CT0.0", "am.CT0.1", "am.CT0.3", "am.CT0.5", "am.CT0.7", "am.CT0.14", "am.CT0.28", "am.CT0.56", "am.CT3.0", "am.CT3.1", "am.CT3.3", "am.CT3.5", "am.CT3.7", "am.CT3.14", "am.CT3.28", "am.CT3.56", "am.PBB0.0", "am.PBB0.1", "am.PBB0.3", "am.PBB0.5", "am.PBB0.7", "am.PBB0.14", "am.PBB0.28", "am.PBB0.56", "am.PBB3.0", "am.PBB3.1", "am.PBB3.3", "am.PBB3.5", "am.PBB3.7", "am.PBB3.14", "am.PBB3.28", "am.PBB3.56", "am.PBC0.0", "am.PBC0.1", "am.PBC0.3", "am.PBC0.5", "am.PBC0.7", "am.PBC0.14", "am.PBC0.28", "am.PBC0.56", "am.PBC3.0", "am.PBC3.1", "am.PBC3.3", "am.PBC3.5", "am.PBC3.7", "am.PBC3.14", "am.PBC3.28", "am.PBC3.56")
        new.value <- c("c1", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c2", "c21", "c22", "c23", "c24", "c25", "c26", "c27", "c3", "c32", "c33", "c34", "c35", "c36", "c37", "c38", "c4", "c41", "c42", "c43", "c44", "c45", "c46", "c47", "c5", "c51", "c52", "c53", "c54", "c55", "c56", "c57", "c6", "c61", "c62", "c63", "c64", "c65", "c66", "c67")
        data.am$NEW <- new.value [match(data.am$LABEL, indata)]
        #data.am$unik <- paste ('c', 1:9, 21:29, 31:39, 41:49, sep = '.')
        View(data.am)
        
        library (ggplot2)
        library(RColorBrewer)
        mypalette<-brewer.pal(n = 7, name = "Oranges")
        mypalette <- c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0", "#EDF8E9")
        mypalette <- c("#8C2D04", "#D94801", "#F16913", "#FD8D3C", "#FDAE6B", "#FDD0A2", "#FEEDDE")
        image(1:7,1,as.matrix(1:7),col=mypalette,xlab="Greens (sequential)",
              ylab="",xaxt="n",yaxt="n",bty="n")
        
        ggplot (data = data.am, aes(x=NEW, y=Level_2, fill=FACTOR)) + 
                geom_tile() +
                scale_fill_manual(values = mypalette)
        
                scale_color_gradient(low = "white", high = "red")
        
        
        textcol <- "black"
        ggplot (data.am, aes(x=NEW, y=Level_2, fill=FACTOR)) +
                geom_tile() +
                geom_tile (colour="white", size=0.25, show.legend = FALSE) +
                labs (x="", y="", title="") +
                scale_y_discrete (expand=c(0,0)) +
                scale_x_discrete (expand=c(0,0)) +
                scale_fill_manual (name = "log10(pi)",
                                   values=c("#d53e4f","#f46d43","#fdae61", "#fee08b","#e6f598","#abdda4", "#ddf1da"),
                                   #values = mypalette,
                                   #values = colorRampPalette(c('black','white'))(7),
                                   na.value="grey90") +
                coord_fixed() +
                theme_grey (base_size=12) +
                theme(  legend.title=element_text(size = 8, colour = "black"),
                        legend.margin = grid::unit(0,"cm"),
                        legend.text=element_text(colour=textcol,size=7,face="bold"),
                        legend.key.height=grid::unit(0.8,"cm"),
                        legend.key.width=grid::unit(0.3,"cm"),
                        legend.justification = "top",
                        axis.text.x=element_text(size=10,colour=textcol,angle = 90),
                        axis.text.y=element_text(vjust = 0.2,colour=textcol),
                        axis.ticks=element_line(size=0.4),
                        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
                        plot.background=element_blank(),
                        panel.border=element_blank())


#################################################
#################################################        
## LONG to wide format for KEGG L2:

## a. Separated by dynamic
am.average <- subset(K2.sort.average, dynamic == "am")
am.average$label <- paste (am.average$Treatment, am.average$days, sep = '.')

un.average <- subset(K2.sort.average, dynamic == "un")
un.average$label <- paste (am.average$Treatment, am.average$days, sep = '.')

## b. long to wide for both dynamics:
am.K2_mat <- spread (data = am.average[,c(2,6,10)],
                     key = "label",
                     value = "log10x")

un.K2_mat <- spread (data = un.average[,c(2,6,10)],
                     key = 'label',
                     value = 'log10x')
View(un.K2_mat)

## Calculate the total assigned reads for all samples.
sum(rowSums(KEGG_L2_READS[,2:ncol(KEGG_L2_READS)]))

## save the MEAN dataframe
write.table(KEGG_L2_READS, '/home/mario/Documents/Doctorado/PROYECTOS/Maryo_Hg/ammonia/downstream/picrust_3000/LUC_PICRUSt/KEGGL2_dynam_day_treatm_ASSIGNED_READS.txt',
            sep = '\t', row.names = F)



#################################################
#################################################
## Editing to final data to create heatmap.
##6. Working with amended dynamic:
am.average$observations <- paste(am.average$Treatment, am.average$days, sep = '.')
View(am.average)
am.average <- am.average[,c(2,6,10)]
library(tidyr)
am.data <- spread (data = am.average,
                   key = "Level_2",
                   value = "logproptrans")
View(am.data)
getwd()
#        write.table(am.data, '/home/mario/Documents/Doctorado/PROYECTOS/Maryo_Hg/ammonia/downstream/picrust_3000/me_working/am_kegg2_log+1_to_heatmap.csv',
#                    sep = ',', row.names = F) #was manually edited; be careful.

## 7. Working with unamended dynamic:
un.average$Observations <- paste (un.average$Treatment, un.average$days, sep = '.')
View(un.average)
un.average <- un.average[,c(2,6,10)]
library(tidyr)                        
un.data <- spread (data = un.average,
                   key = 'Level_2',
                   value = 'logproptrans')        
View (un.data)
#        write.table(am.data, '/home/mario/Documents/Doctorado/PROYECTOS/Maryo_Hg/ammonia/downstream/picrust_3000/me_working/un_kegg2_log+1_to_heatmap.csv',
#                    sep = ',', row.names = F) #was manually edited; be careful.
