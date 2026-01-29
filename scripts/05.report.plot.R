#!/usr/bin/env Rscript

############################################################################
#
# Phasing Accuracy Report - Chromosome 6
#
############################################################################

library(tidyverse)
library(ggplot2)
library(scales)

# 1. Path Settings
parent_dir <- "/home/jennifer/02_datas/04_data_processing_trios/01_intermediate/switch"
out_dir <- file.path(parent_dir, "report")
html_path <- file.path(out_dir, "phasing_report.html")
combined_plot_name <- "combined_phasing_plot.png"
combined_plot_name2 <- "combined_phasing_plot2.png"

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# 2. Table Data (with Standardized Labels)
metrics_df <- data.frame(
  Sample = c("HG00733", "HG01109", "HG01243",
            "HG02055", "HG02080", "HG02145",
            "HG02723", "HG02818", "HG03098",
            "HG03486", "HG03492", "NA18906",
            "NA19240", "NA20129"),
  WHswitchrate = c("0/8323 (0%)", "2/9502 (0.02%)", "0/7931 (0%)",
            "3/8174 (0.04%)", "2/8737 (0.02%)", "2/8896 (0.02%)",
            "8/8827 (0.09%)", "2/9048 (0.02%)", "14/8990 (0.16%)",
            "2/7489 (0.03%)", "4/8621 (0.05%)", "2/8766 (0.02%)",
            "3/6311 (0.05%)", "3/8353 (0.04%)"),
  WhatsHap = c("0/8323 (0%)", "0/9502 (0%)", "0/7931 (0%)",
            "1233/8174 (15.08%)", "2118/8737 (24.24%)", "1/8896 (0.01%)",
            "4/8827 (0.05%)", "1/9048 (0.01%)", "12/8990 (0.13%)",
            "1/7489 (0.01%)", "2/8621 (0.02%)", "1/8766 (0.01%)",
            "1222/6311 (13.36%)", "1191/8353 (14.26%)"),
  InHouse = c("0/8340 (0%)", "0/9533 (0%)", "0/7953 (0%)",
            "1240/8193 (15%)", "2128/8773 (24.25%)", "1/8918 (0.01%)",
            "4/8866 (0.045%)", "1/9069 (0.011%)", "12/8926 (0.134%)",
            "1/7507 (0.013%)", "2/8645 (0.02%)", "1/8786 (0.01%)",
            "1231/6334 (19.43%)", "1193/8363 (14.26%)")
  ) %>%
  mutate(
    # Extract the %
    sort_key = as.numeric(str_extract(WhatsHap, "[0-9.]+(?=%)")),
    
    # Extract the number of the %
    inhouse_val = as.numeric(str_extract(InHouse, "[0-9.]+(?=%)")),
    Pct_Label = paste0(Sample, " (", sprintf("%.2f", inhouse_val), "%)")
  ) %>%
  arrange(sort_key) %>%
  select(-sort_key, -inhouse_val)

# 3. Start HTML Building
html_lines <- c(
  "<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'>",
  "<title>Phasing Accuracy Report - Chr6</title>",
  "<style>",
  "  body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; background: #f8f9fa; color: #333; line-height: 1.6; }",
  "  .container { max-width: 1200px; margin: auto; }",
  "  h1 { color: #2c3e50; text-align: center; border-bottom: 2px solid #3498db; padding-bottom: 10px; }",
  "  table { border-collapse: collapse; width: 100%; margin: 30px 0; background: white; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }",
  "  th, td { border: 1px solid #dee2e6; padding: 12px; text-align: center; }",
  "  th { background-color: #34495e; color: white; font-weight: 600; }",
  "  tr:nth-child(even) { background-color: #f2f2f2; }",
  "  .description-box { background-color: #e8f4fd; border-left: 5px solid #3498db; padding: 15px; margin: 20px 0; border-radius: 0 4px 4px 0; }",
  "  .plot-container { background: white; padding: 25px; border: 1px solid #e1e4e8; margin-bottom: 30px; border-radius: 10px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }",
  "  img { width: 100%; height: auto; border-radius: 4px; border: 1px solid #eee; }",
  "  .footer { text-align: center; font-size: 0.9em; color: #7f8c8d; margin-top: 50px; }",
  "</style></head><body><div class='container'>",
  "<h1>Phasing Accuracy Analysis: MHC region, chromosome 6</h1>",
  "<h3>Summary Metrics</h3>",
  "<div class='description-box'>",
  "<strong>Hamming Distance:</strong> Is the number of heterozygous positions where the phase orientation differs between the predicted (HLA-mapper) and the truth (HRPC) dataset.",
  "</div>"
)

# 4. Add the Table
html_lines <- c(html_lines,
                "<table><thead><tr><th>Sample</th>",
                "<th>WhatsHap (Switch error/common heterozygous variants, %)</th>",
                "<th>WhatsHap (Hamming distance/common heterozygous variants, %)</th>",
                "<th>InHouse (Hamming distance/common heterozygous variants, %)</th>",
                "</tr></thead><tbody>")

for(i in 1:nrow(metrics_df)) {
  row_html <- paste0("<tr><td><strong>", metrics_df$Sample[i], "</strong></td>",
                     "<td>", metrics_df$WHswitchrate[i], "</td>",
                     "<td>", metrics_df$WhatsHap[i], "</td>",
                     "<td>", metrics_df$InHouse[i], "</td></tr>")
  html_lines <- c(html_lines, row_html)
}
html_lines <- c(html_lines, "</tbody></table><hr>")

# 5. Load and Combine Data
data_list <- list()
for(sample_id in metrics_df$Sample) {
  sample_file <- file.path(parent_dir, sample_id, paste0(sample_id, ".hrpc.hlamapper.switch.errors.tsv"))
  
  if (file.exists(sample_file)) {
    current_label <- metrics_df$Pct_Label[metrics_df$Sample == sample_id]
    
    temp_df <- read.table(sample_file, header=T, sep="\t") %>%
      mutate(pos = as.numeric(str_extract(idcomp, "(?<=:)[0-9]+")),
             PlotLabel = current_label) %>%
      select(PlotLabel, pos, status_error)
    data_list[[sample_id]] <- temp_df
  }
}

master_df <- bind_rows(data_list)
master_df$PlotLabel <- factor(master_df$PlotLabel, levels = metrics_df$Pct_Label)

# 6. Generate Plot
global_min_pos <- min(master_df$pos, na.rm = TRUE)
global_max_pos <- max(master_df$pos, na.rm = TRUE)

p_combined <- ggplot(master_df, aes(x = pos, y = 1)) +
  annotate("segment", x = global_min_pos, xend = global_max_pos, y = 1, yend = 1, 
           color = "grey90", linewidth = 0.5) +
  geom_point(aes(color = status_error), size = 2, alpha = 0.8, shape = 124) + 
  facet_grid(PlotLabel ~ ., switch = "y") + 
  scale_color_manual(values = c("MATCH" = "#27ae60", "DIFF" = "#e74c3c")) +
  scale_x_continuous(labels = label_comma(), expand = c(0.01, 0)) +
  labs(title = "Phasing error map", x = "Position on Chromosome 6 (bp)", y = "Sample (Hamming distance %)", color = "Status") +
  theme_minimal() +
  theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      # O uso de family='mono' garante que os números ocupem o mesmo espaço
      strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, 
                                       face = "bold", size = 9, family = "mono"),
      strip.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      panel.spacing = unit(0.1, "lines")
  )

plot_height <- max(4, length(unique(master_df$Sample)) * 1)
ggsave(filename = combined_plot_name, path = out_dir, plot = p_combined, 
       width = 12, height = plot_height, dpi = 300, bg = "white")

# 7. Generate Plot with genotyping errors
global_min_pos <- min(master_df$pos, na.rm = TRUE)
global_max_pos <- max(master_df$pos, na.rm = TRUE)

p_combined <- ggplot(master_df, aes(x = pos, y = 1)) +
  annotate("segment", x = global_min_pos, xend = global_max_pos, y = 1, yend = 1, 
           color = "grey90", linewidth = 0.5) +
  geom_point(aes(color = status_error), size = 2, alpha = 0.8, shape = 19) + 
  facet_grid(PlotLabel ~ ., switch = "y") + 
  scale_color_manual(values = c("MATCH" = "#27ae60", "DIFF" = "#e74c3c", "ERROR" = "#000000")) +
  scale_x_continuous(labels = label_comma(), expand = c(0.01, 0)) +
  labs(title = "Phasing error map", x = "Position on Chromosome 6 (bp)", y = "Sample (Hamming distance %)", color = "Status") +
  theme_minimal() +
  theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, 
                                       face = "bold", size = 9, family = "mono"),
      strip.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      panel.spacing = unit(0.1, "lines")
  )

plot_height <- max(4, length(unique(master_df$Sample)) * 1)
ggsave(filename = combined_plot_name2, path = out_dir, plot = p_combined, 
       width = 12, height = plot_height, dpi = 300, bg = "white")

# 8. Finalize HTML
html_lines <- c(html_lines, 
  "<h3>Comparative plot</h3>",
  "<div class='plot-container'>",
  paste0("<img src='", combined_plot_name2, "' alt='Combined Switch Error Plot'>"),
  "</div>"
)

html_lines <- c(html_lines, 
                "<div class='footer'>Report generated on ", as.character(Sys.Date()), "</div>",
                "</div></body></html>")

writeLines(html_lines, con = html_path)

cat("\n--- Report successfully generated ---\n")
cat("Location:", html_path, "\n")