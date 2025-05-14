# create_database.R


# Libraries Load 

library(DBI)        
library(RSQLite)   
library(readr)     
library(dplyr)      
library(stringr)  
library(readxl)

# Configuration (Settings)

# 1. CSV Files List:

xlsx_files <- list(
  
  S1 = "data/Table_S1.xlsx",
  S2 = "data/Table_S2.xlsx",
  S3 = "data/Table_S3.xlsx",
  S4 = "data/Table_S4.xlsx",
  S5 = "data/Table_S5.xlsx",
  S6 = "data/Table_S6.xlsx",
  S7 = "data/Table_S7.xlsx",
  S8 = "data/Table_S8.xlsx",
  S9 = "data/Table_S9.xlsx",
  S10 = "data/Table_S10.xlsx",
  S11 = "data/Table_S11.xlsx",
  S12 = "data/Table_S12.xlsx"
)

database_path <- "data/chd_supplementary_data.sqlite"

message("Creating/Updating database: ", database_path)
con <- dbConnect(RSQLite::SQLite(), database_path)

success_count <- 0
for (table_name in names(xlsx_files)) {
  file_path <- xlsx_files[[table_name]]
  message("Processing: ", file_path, " into table: ", table_name)
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path, ". Skipping.")
    next
  }
  
  tryCatch({
    df <- read_excel(
      file_path,
      skip = 1,
      col_names = TRUE
    )
    
    if (nrow(df) == 0 && ncol(df) == 0) { 
      message(" -> Table ", table_name, " is empty or unreadable after skipping the first row. Skipping database write.")
      next
    }
    
    # Standardize chromosome columns
    chr_cols_pattern <- "^(CHR|Chr|chr_enh|CHD_Enh_chr|HEDD_chr|Chr_rCHD_SNP)$" 
    chr_cols <- names(df)[grepl(chr_cols_pattern, names(df), ignore.case = TRUE)]
    
    for(col in chr_cols) {
      if(any(!is.na(df[[col]]))) {
        current_col_values <- as.character(df[[col]])
        needs_prefix <- !startsWith(current_col_values, "chr") & !is.na(current_col_values)
        current_col_values[needs_prefix] <- paste0("chr", current_col_values[needs_prefix])
        df[[col]] <- current_col_values
      }
    }
    
    # Clean rsIDs
    if(table_name == "S2" && "SNP_ID_CURRENT" %in% names(df)) {
      df <- df %>%
        mutate(SNP_ID_CURRENT = gsub("-\\?$", "", as.character(SNP_ID_CURRENT))) %>%
        rename(rs_ID = SNP_ID_CURRENT)
    } else if ("rs_ID" %in% names(df)) {
      df <- df %>% mutate(rs_ID = as.character(rs_ID)) 
    }
    
    
    # Add Combined Coordinate Columns
    if (table_name == "S4" && all(c("Chr", "Start_enh", "End_enh") %in% names(df))) {
      df <- df %>% mutate(Enhancer_Coord = paste0(Chr, ":", Start_enh, "-", End_enh)) %>% relocate(Enhancer_Coord, .after = End_enh)
    } else if (table_name == "S6" && all(c("Chr", "Start_enh", "End_enh") %in% names(df))) {
      df <- df %>% mutate(Enhancer_Coord = paste0(Chr, ":", Start_enh, "-", End_enh)) %>% relocate(Enhancer_Coord, .after = End_enh)
    } else if (table_name == "S7" && all(c("CHD_Enh_chr", "CHD_Enh_start", "CHD_Enh_end") %in% names(df))) {
      df <- df %>% mutate(CHD_Enh_Coord = paste0(CHD_Enh_chr, ":", CHD_Enh_start, "-", CHD_Enh_end)) %>% relocate(CHD_Enh_Coord, .after = CHD_Enh_end)
      if (all(c("HEDD_chr", "HEDD_start", "HEDD_end") %in% names(df))) {
        df <- df %>% mutate(HEDD_Coord = paste0(HEDD_chr, ":", HEDD_start, "-", HEDD_end)) %>% relocate(HEDD_Coord, .after = HEDD_end)
      }
    } else if (table_name == "S8" && all(c("chr_enh", "start_enh", "end_enh") %in% names(df))) {
      df <- df %>% mutate(Enhancer_Coord = paste0(chr_enh, ":", start_enh, "-", end_enh)) %>% relocate(Enhancer_Coord, .after = end_enh)
    } else if (table_name == "S12") {
      if (all(c("Chr_enh", "Start_enh", "End_enh") %in% names(df))) {
        df <- df %>% mutate(Enhancer_Coord_S12 = paste0(Chr_enh, ":", Start_enh, "-", End_enh)) %>% relocate(Enhancer_Coord_S12, .after = End_enh)
      }
      if (all(c("Chr_rCHD_SNP", "rCHD_SNP_Position") %in% names(df))) {
        df <- df %>% mutate(rCHD_SNP_Coord = paste0(Chr_rCHD_SNP, ":", `rCHD_SNP_Position`)) %>% relocate(rCHD_SNP_Coord, .after = `rCHD_SNP_Position`)
      }
    } else if (table_name == "S2" && all(c("CHR_ID", "CHR_POS") %in% names(df))) {
      df <- df %>% mutate(SNP_Coord = paste0(CHR_ID, ":", CHR_POS)) %>% relocate(SNP_Coord, .after = CHR_POS)
    } else if (table_name == "S3" && all(c("CHR", "Start") %in% names(df))) {
      df <- df %>% mutate(Variant_Coord = paste0(CHR, ":", Start)) %>% relocate(Variant_Coord, .after = End)
    } else if (table_name == "S11" && all(c("Chrom", "SNP_start") %in% names(df))) {
      df <- df %>% mutate(Conserved_SNP_Coord = paste0(Chrom, ":", SNP_start)) %>% relocate(Conserved_SNP_Coord, .after = SNP_end)
    }
    
    if(table_name == "S10") {
      original_name_col_s10  <- names(df)[grepl("Name.*conserved.*CHD.*enhancer", names(df), ignore.case = TRUE)][1]
      
      
      if (!is.na(original_name_col_s10) && original_name_col_s10 %in% names(df)) {
        message("-> Extracting coords for S10 from column: '", original_name_col_s10, "'")
        df <- df %>%
          mutate(
            Extracted_Part_S10 = stringr::str_extract(!!sym(original_name_col_s10), "chr[^-]+-\\d+-\\d+"),
            Enhancer_Coord_S10 = stringr::str_replace(Extracted_Part_S10, "-", ":")
          ) %>%
          relocate(Extracted_Part_S10, Enhancer_Coord_S10, .after = all_of(original_name_col_s10))
        
        if(all(is.na(df$Enhancer_Coord_S10[!is.na(df[[original_name_col_s10]])]))) {
          warning("Coordinate extraction for S10 might have failed for non-NA original names. Check column '", original_name_col_s10, "' and regex pattern.")
        } else {
          message(" -> Successfully processed Enhancer_Coord_S10 column for S10.")
        }
      } else {
        warning("Original name column for S10 (expected like 'Name_of_conserved_CHD_enhancer') not found. Cannot extract coordinates.")
      }
    }
    
    dbWriteTable(con, table_name, df, overwrite = TRUE, row.names = FALSE)
    success_count <- success_count + 1
    message(" -> Successfully wrote table: ", table_name)
    
  }, error = function(e) {
    warning("Error processing ", file_path, ": ", e$message)
  })
}

message("\nCreating/Updating database indexes...")
index_commands <- list(
  "CREATE INDEX IF NOT EXISTS idx_S1_Term ON S1 (Term);",
  "CREATE INDEX IF NOT EXISTS idx_S2_rsID ON S2 (rs_ID);",
  "CREATE INDEX IF NOT EXISTS idx_S2_SNPCoord ON S2 (SNP_Coord);",
  "CREATE INDEX IF NOT EXISTS idx_S2_MAPPED_GENE ON S2 (MAPPED_GENE);",
  "CREATE INDEX IF NOT EXISTS idx_S3_rsID ON S3 (rs_ID);",
  "CREATE INDEX IF NOT EXISTS idx_S3_VariantCoord ON S3 (Variant_Coord);",
  "CREATE INDEX IF NOT EXISTS idx_S3_GeneSymbol ON S3 (GeneSymbol);",
  "CREATE INDEX IF NOT EXISTS idx_S4_EnhCoord ON S4 (Enhancer_Coord);",
  "CREATE INDEX IF NOT EXISTS idx_S4_rsID ON S4 (rs_ID);",
  "CREATE INDEX IF NOT EXISTS idx_S6_EnhCoord ON S6 (Enhancer_Coord);",
  "CREATE INDEX IF NOT EXISTS idx_S6_rsID ON S6 (rs_ID);",
  "CREATE INDEX IF NOT EXISTS idx_S7_CHDEnhCoord ON S7 (CHD_Enh_Coord);",
  "CREATE INDEX IF NOT EXISTS idx_S8_EnhCoord ON S8 (Enhancer_Coord);",
  "CREATE INDEX IF NOT EXISTS idx_S9_seqname ON S9 (`sequence_name (CHD-enhancer coordinates)`);",
  "CREATE INDEX IF NOT EXISTS idx_S10_name ON S10 (`Name of conserved CHD-enhancer`);",
  "CREATE INDEX IF NOT EXISTS idx_S10_EnhCoord ON S10 (Enhancer_Coord_S10);",
  "CREATE INDEX IF NOT EXISTS idx_S11_rsID ON S11 (rs_ID);",
  "CREATE INDEX IF NOT EXISTS idx_S11_ConsSNPCoord ON S11 (Conserved_SNP_Coord);",
  "CREATE INDEX IF NOT EXISTS idx_S12_EnhCoord ON S12 (Enhancer_Coord_S12);",
  "CREATE INDEX IF NOT EXISTS idx_S12_Genes ON S12 (Genes);"
)

index_success <- 0
for (cmd in index_commands) {
  message("Executing: ", cmd)
  tryCatch({
    dbExecute(con, cmd)
    index_success <- index_success + 1
  }, error = function(e) {
    warning("Index creation/update failed for: ", cmd, " Error: ", e$message)
  })
}

dbDisconnect(con)
message("\nDatabase creation/update complete.")
message("Successfully wrote ", success_count, " tables.")
message("Successfully created/updated ", index_success, " indexes.")
message("Database file saved to: ", database_path)