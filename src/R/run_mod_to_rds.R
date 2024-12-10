# Library
library("dplyr")
library("stringr")
library("data.table")

# Set input directory path
BED_DIR<- "/nfs/research/sds/sds-1kg/GRCh38/haplotypesv2"

# List bed files in directory
list_bed_files <- list.files(path = BED_DIR, pattern = "mod.bed$", full.names = TRUE) 

# Create metadata of bed files
samples.df<- data.frame(bed_file = list_bed_files)
samples.df <- samples.df %>% 
    mutate(
        sample_name = list_bed_files %>% str_remove(".+haplotypesv2/"),
        sample_name = sample_name     %>% str_remove(".genes.+"),
        hap         = sample_name %>% str_remove(".+hap"),
        sample_name = sample_name %>% str_remove(".hap.+")
    )

# Read bed files into list of dfs
bed.list <- samples.df %>% apply(1,function(sample_bed){
    tmp.df <- read.table(file = sample_bed["bed_file"], header = FALSE, sep = "\t", quote = "") %>% mutate(hap = sample_bed["hap"],sample_name = sample_bed["sample_name"] )
}) 

# Bind list of dfs
bed.df <- bed.list %>% dplyr::bind_rows()

# Select only relevant columns of df
bed.df <- bed.df %>% dplyr::select(V1,V2,V3,V4,V5,V6,V10,V11,V12,V13,V14,hap,sample_name) 

# Rename columns
bed.df <- bed.df %>% 
    dplyr::rename( 
        "chr"       = "V1", 
        "start"     = "V2", 
        "end"       = "V3", 
        "mod"       = "V4", 
        "cov"       = "V5", 
        "strand"    = "V6", 
        "mod_stats" = "V10", 
        "gene_chr"  = "V11",        
        "gene_start"= "V12",        
        "gene_end"  = "V13",        
        "gene"      = "V14"
    )

# Split string from modification stats
mod_stats.df <- bed.df %>% 
    pull("mod_stats") %>%
    str_split(" ", simplify = TRUE) %>% 
    as.data.frame(   ) %>% 
    dplyr::select(-V1) %>%
    dplyr::rename( 
        "mod_ratio"            ="V2",
        "mod_count"            ="V3",
        "canonical_count"      ="V4",
        "othermod_count"       ="V5",
        "read_deletion_count"  ="V6",
        "failed_count"         ="V7",
        "ntsubstitution_count" ="V8",
        "nocall_count"         ="V9" 
    )

# Append modification stats
bed.df <- cbind( bed.df, mod_stats.df )

# Drop mod_stats column
bed.df <- bed.df %>% dplyr::select(-mod_stats) %>% mutate(mod_ratio = as.numeric(mod_ratio))

# Save RDS object 
saveRDS(bed.df, file = "samples_cpg_pileup.rds")
