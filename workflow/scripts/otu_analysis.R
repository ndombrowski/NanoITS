#!/usr/bin/env Rscript

#rm(list=ls())


#___________________________________
#  Defaults (change if needed)
#___________________________________
#samples with less than x reads are discarded
max_reads_per_sample <- 1

#otus that occur just once across all samples are discarded
otu_cutoff <- 1

#___________________________________
#  General setup
#___________________________________

# Use commandArgs() to get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to print help information
print_help <- function() {
  cat("Usage: Rscript otu_analysis.R <file_paths> <marker_id> <output_folder>\n")
  cat("Arguments:\n")
  cat("  <file_paths>: File paths for input files (wildcards allowed)\n")
  cat("  <marker_id>: Name of the marker investigated \n")
  cat("  <output_folder>: Output folder for results\n")
  cat("  -h, --help: Print this help message\n")
  quit(status = 0)
}

# Check if help option is provided
if (any(grepl("--help|-h", args))) {
  print_help()
}

# Check if the correct number of arguments is provided
if (length(args) < 3) {
  stop("Usage: script.R <otu_table1> <otu_table2> ... <marker_id> <output_folder> ")
}


# Extract file paths from command line arguments
file_paths <- args[1:(length(args) - 2)]
marker_id <- args[length(args) - 1]
output_folder <- args[length(args)]


#___________________________________
#  Loading Libraries
#___________________________________
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(microbiome))
suppressPackageStartupMessages(library(wesanderson))



#___________________________________
#  Define custom functions
#___________________________________
# Define a custom labeller function
custom_labeller <- function(variable, value) {
  if (variable == "barcode") {
    # Shorten barcode labels
    return(paste0("bc", substr(value, nchar(value)-1, nchar(value))))
  } else {
    return(value)
  }
}

#define a custom theme for plotting
custom_theme <- function() {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "#FFFFFF", color = "black"),
    axis.text.x = element_text(size = 7)
  )
}

#define some color palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette2 <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
              "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
              "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
              "#8A7C64", "#599861")
col_vector <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 
                '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
                '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')



#___________________________________
#  Read in otu tables
#___________________________________
# Read otu tables
methods <- character(0)

for (i in 1:length(file_paths)){
  table <- read.table(file_paths[i], header = T, sep = '\t', comment = "", check.names = FALSE)
  print(colnames(table))
  method <- gsub(".*/([^/]+)/.*", "\\1", file_paths[i])
  methods <- c(methods, method)
  colnames(table) <- paste0(colnames(table), "_", method)
  names(table)[1] <- "taxid"
  ifelse (i == 1, merged_otu_table <- table, 
          merged_otu_table <- merge(merged_otu_table, table, by = "taxid",  all=TRUE))
} 

# Replace NA with 0
merged_otu_table[is.na(merged_otu_table)] <- 0

# Restore row names
rownames(merged_otu_table) <- merged_otu_table$taxid
merged_otu_table$taxid <- NULL
dim(merged_otu_table)

write.table(merged_otu_table, file = paste0(output_folder, "/tables/", marker_id, "_otu_table.txt"), sep = "\t", col.names=NA)


#___________________________________
#  Generate sample mapping file
#___________________________________
#extract sample ids
metadata_combined <- as.data.frame(colnames(merged_otu_table))
colnames(metadata_combined) <- "name"

#create column to separate methods
metadata_combined <- metadata_combined |> 
  separate(name, c("barcode", "method"), sep = "_", extra = "merge", remove = FALSE )

#generate rownames
row.names(metadata_combined) <- metadata_combined$name


#___________________________________
#  Generate taxonomy mapping file
#___________________________________
#Generate tax file
#extract taxonomy string
temp <- as.data.frame(rownames(merged_otu_table))
colnames(temp) <- "OTU"

#separate the taxonomic headers                      
taxonomy_file <- temp |> 
  distinct(OTU) |> 
  separate(OTU,
           c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = ";", remove = FALSE) |>
  column_to_rownames(var = "OTU")


#___________________________________
#  Prepare phyloseq data object
#___________________________________
#convert to phyloseq object
OTU = otu_table(merged_otu_table, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy_file))
physeq = phyloseq(OTU, TAX, sample_data(metadata_combined))
physeq

#___________________________________
#  Data filtering
#___________________________________
#remove taxa without tax assignment
physeq <- subset_taxa(physeq, Phylum != "NA")

#remove samples with less than 20 reads
physeq <- prune_samples(sample_sums(physeq) >= max_reads_per_sample, physeq)

#remove singletons
physeq <- filter_taxa(physeq, function(x){sum(x > 0) > otu_cutoff}, prune = TRUE)

#print 
physeq_printing <- psmelt(physeq)
physeq_printing <- spread(physeq_printing[,c("Sample", "OTU", "Abundance")], Sample, Abundance, fill = 0)

write.table(physeq_printing, file = paste0(output_folder, "/tables/", marker_id, "_otu_table_filtered.txt"), sep = "\t", row.names = FALSE)



#___________________________________
#  Data normalization
#___________________________________
#Perform three different normalization methods
#right now, only the relative abundance table is used
physeq_rel_abundance <- microbiome::transform(physeq, "compositional")
physeq_clr <- microbiome::transform(physeq, "clr")
physeq_rarified <- rarefy_even_depth(physeq)


#___________________________________
#  Generate barplots for count data
#___________________________________
#plot relative abundance data for each taxonomic level
for (level in colnames(taxonomy_file)){
  #condense at specific tax rank
  grouped_taxa <- tax_glom(physeq, level)

  #find topx taxa and group remaining into "other"
  top_taxa <- names(sort(taxa_sums(grouped_taxa), TRUE)[1:19])
  other_taxa <- names(taxa_sums(grouped_taxa))[which(!names(taxa_sums(grouped_taxa)) %in% top_taxa)]

  merged_physeq = merge_taxa(grouped_taxa, other_taxa, 2)

  #transform to dataframe 
  df <- psmelt(merged_physeq) 
  names(df)[names(df) == level] <- "tax_level"

  #Add NAs to an Other taxonomy level
  df$tax_level[which(is.na(df$tax_level))] <- "Other"
  df$tax_level <- ifelse(df$tax_level %in% c("", "NA", "N/A", "missing", "other_missing_representation"),
                       "Other", df$tax_level)

  #sort taxa labels by abundance
  df$tax_level <- as.factor(df$tax_level)

  labels <- df |>
    group_by(tax_level) |>
    summarise(sum = sum(Abundance)) |>
    arrange(desc(sum))

  #Sort levels and ensure that Other category is the last category
  desired_levels <- setdiff(labels$tax_level, "Other")
  df$tax_level2 <- factor(df$tax_level, levels = c(desired_levels, "Other"))

  #get max value for axis height
  max <- df |>
    group_by(Sample) |>
    summarise(max_value = sum(Abundance))

  max <- max(max$max_value)

  #generate color scheme
  cols <- palette2[1:length(unique(df$tax_level2))]
  cols[levels(df$tax_level2) == "Other"] <- "#CCCCCC"

  #generate barplots
  counts_barplot <-
    ggplot(df, aes(x = Sample, y = Abundance, fill = tax_level2) ) +
      geom_bar(position = "stack", stat = "identity") +
      labs(y = "Total counts", x = "", title = paste0("Total counts at ", level, " rank for ", marker_id)) +
      scale_fill_manual(name = paste0(level,"_rank"), labels = levels(df$tax_level2), values = cols) +
      facet_wrap(~barcode, scales = "free_x", ncol = 12) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max * 1.01)) +
      custom_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

  ggsave(paste0(output_folder,"/plotting/", marker_id, "/", marker_id, "_Barplot_counts_", level, "_rank.pdf"), plot=counts_barplot, device="pdf")
}


#___________________________________
#  Generate barplots for relative abundance data
#___________________________________
#plot relative abundance data for each taxonomic level
for (level in colnames(taxonomy_file)){
  #condense at specific tax rank
  grouped_taxa <- tax_glom(physeq_rel_abundance, level)

  #find topx taxa and group remaining into "other"
  top_taxa <- names(sort(taxa_sums(grouped_taxa), TRUE)[1:19])
  other_taxa <- names(taxa_sums(grouped_taxa))[which(!names(taxa_sums(grouped_taxa)) %in% top_taxa)]

  merged_physeq = merge_taxa(grouped_taxa, other_taxa, 2)

  #transform to dataframe 
  df <- psmelt(merged_physeq) 
  names(df)[names(df) == level] <- "tax_level"

  #Add NAs to an Other taxonomy level
  df$tax_level[which(is.na(df$tax_level))] <- "Other"
  df$tax_level <- ifelse(df$tax_level %in% c("", "NA", "N/A", "missing", "other_missing_representation"),
                       "Other", df$tax_level)

  #sort taxa labels by abundance
  df$tax_level <- as.factor(df$tax_level)

  labels <- df |>
    group_by(tax_level) |>
    summarise(sum = sum(Abundance)) |>
    arrange(desc(sum))

  #Sort levels and ensure that Other category is the last category
  desired_levels <- setdiff(labels$tax_level, "Other")
  df$tax_level2 <- factor(df$tax_level, levels = c(desired_levels, "Other"))

  #generate color scheme
  cols <- palette2[1:length(unique(df$tax_level2))]
  cols[levels(df$tax_level2) == "Other"] <- "#CCCCCC"

  #generate barplots
  ra_barplot <-
    ggplot(df, aes(x = Sample, y = Abundance, fill = tax_level2) ) +
      geom_bar(position = "stack", stat = "identity") +
      labs(y = "Relative abundance", x = "", title = paste0("Relative abundance at ", level, " rank for ", marker_id)) +
      scale_fill_manual(name = paste0(level,"_rank"), labels = levels(df$tax_level2), values = cols) +
      facet_wrap(~barcode, scales = "free_x", ncol = 12) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
      custom_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
      
  ggsave(paste0(output_folder,"/plotting/", marker_id, "/", marker_id, "_Barplot_ra_", level, "_rank.pdf"), plot=ra_barplot, device="pdf")
}

