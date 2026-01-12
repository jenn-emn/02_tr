#/bin/bash/Rscript

############################################################################
#
# Plot switch blocks into chromossome 6
#
############################################################################

library(tidyverse)
library(ggplot2)
library(tidyr)



# PATHS

outpath="/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/switch"

#pathkin=paste0(outpath,'/kin')
#dir.create(pathkin, recursive = TRUE)


# names lists
names=c("HG00733", "HG01109", "HG01243", "HG02055", "HG02080", "HG02145", "HG02723", "HG02818", "HG03098", "HG03486", "HG03492", "NA18906", "NA19240", "NA20129")
names=c("HG00733")

for (i in seq(1,length(names))) {

    path_plot = paste0(outpath,'/',names[i],'/',names[i],'.plot.png')
    path_switch_sizes = paste0(outpath,'/',names[i],'/',names[i],'.switch.sizes.tsv')
    # head /home/jennifer/02_datas/04_data_processing_trios/01_intermediate/switch/HG00733/HG00733.switch.sizes.tsv 
    # start	end	block_id	block_size
    # chr6:29948260:C:T	chr6:29948260:C:T	1	1
    # chr6:29949558:G:GTAAA	chr6:29949558:G:GTAAA	2	1
    
    # data
    data <- read.table(file= path_switch_sizes, header=T, sep="\t" )

    parse_coord <- function(coord_string) { as.numeric(sapply(strsplit(as.character(coord_string), ":"), `[`, 2)) }

    data_clean <- data %>% mutate( pos_start = parse_coord(start), pos_end = parse_coord(end) )

    # plot
    y_pos <- 1 
    p <- ggplot(data_clean) +
        geom_segment(aes(x = min(pos_start) - 100000, 
                        xend = max(pos_end) + 100000, 
                        y = y_pos, yend = y_pos), 
                    color = "gray80", size = 2) +
        geom_rect(aes(xmin = pos_start - 500, xmax = pos_end + 500, 
                        ymin = y_pos - 0.1, ymax = y_pos + 0.1), 
                    fill = "red", color = "darkred") +
        geom_text(aes(x = (pos_start + pos_end)/2, 
                        y = y_pos + 0.3, 
                        label = paste0(block_id, ":", block_size)), 
                    size = 3, vjust = 0, angle = 45) +
        scale_x_continuous(labels = scales::comma_format(), name = "Posição no Chr6 (bp)") +
        scale_y_continuous(limits = c(0.5, 2), name = "Block: switches") +
        labs(title = paste0("Localization of switched blocks in the Chromosome 6 (",names[i],")"),
            #subtitle = names[i],
            #caption = "Note"
            ) +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank())

    ggsave(path_plot, width = 12, height = 4, dpi = 300)
    
    # scp -P 2205  jennifer@143.107.244.187:/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/switch/HG*/*.plot.png
}


# end