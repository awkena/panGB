#' Build minimal sample list for variant validation
#' This function takes a variant summary and an inventory of available samples to construct a minimal sample list for validating
#' specific variants. It prioritizes samples based on their classification (recurrent parents, donors, de novo references) and ensures a minimum number of samples per variant and allele class.
#' @param variant_summary A list containing variant summary data, including accessions and their associated variants.
#' @param inventory_path A URL pointing to a CSV file containing the inventory of available samples, including their PI numbers and classifications.
#' @param metadata_columns A character vector specifying the columns in the variant summary that contain metadata (default: c("LIB", "Sample", "PInumber")).
#' @param min_req An integer specifying the minimum number of samples required per variant and allele class (default: 3).
#' @param tech_rep An integer specifying the number of technical replicates to include for each sample (default: 2).
#' @param het_tech_rep An integer specifying the number of technical replicates to include for artificial heterozygous samples (default: 3).
#' @param blank_reps An integer specifying the number of empty wells to include for each sample (default: 2). library(shiny)
#' @param sorghum_PCIL_donors A character vector of PI numbers for PCIL donor accessions.
#' @param sorghum_PCIL_recurrent_parents A character vector of PI numbers for PCIL recurrent parent accessions.
#' @param sorghum_de_novo_references A character vector of PI numbers for de novo reference accessions.
#' @param verbose A logical value indicating whether to print detailed messages during processing (default: TRUE).
#' @return A list containing the distribution of samples by variant, the selected samples per variant, and the final minimal sample list for validation.
#' 
#' @importFrom readxl read_excel
#' @importFrom utils read.csv
#' 
build_min_sample_list <- function(
    variant_summary,
    inventory_path,
    metadata_columns = c("LIB", "Sample", "PInumber"),
    min_req = 3,
    tech_rep=2,
    het_tech_rep=3,
    blank_reps=2,
    sorghum_PCIL_donors = c(
      "PI656005","PI656041","PI656012","PI656019","PI655992","PI655993","PI552861","PI656029","PI561072","PI655995",
      "PI24969","PI601816","PI656031","PI576350","PI653617","PI613536","PI656044","PI656046","PI656058","PI655976",
      "PI641836","PI655977","PI561071","PI576376","PI656071","PI656072","PI595739","PI576387","PI656073","PI595714",
      "PI576401","PI576435","PI533965","PI595720","PI533759","PI533961","PI595745","PI597967","PI597980","PI597982",
      "PI656080","PI656081","PI533921","PI656082","PI656085","PI656086","PI656089","PI533924","PI656068","PI533799",
      "PI534127","PI656093","PI656094","PI534037","PI533754","PI533833","PI533957","PI533986","PI534075","PI656096",
      "PI533866","PI534079","PI533769","PI533758","PI533762","PI597950","PI656098","PI533996","PI533830","PI656101",
      "PI533788","PI533755","PI533987","PI533876","PI534092","PI533980","PI656102","PI576337","PI533936","PI597946",
      "PI576332","PI656104","PI533956","PI576366","PI533979","PI533912","PI533937","PI656105","PI534108","PI576333",
      "PI595702","PI576339","PI533985","PI656107","PI533991","PI533949","PI533915","PI597945","PI576359","PI576347",
      "PI656025","PI656026","PI561472","PI607931","PI655975"
    ),
    sorghum_PCIL_recurrent_parents = c("PI655981","PI656031","PI656050","PI565121"),
    sorghum_de_novo_references = c(
      "PI154844","PI156178","PI180348","PI276816","PI276837","PI329301","PI329501","PI513676","PI533766",
      "PI534133","PI564163","PI569459","PI570071","PI576434","PI585966","PI597980","PI655988","PI656015",
      "PI656023","PI656027","PI656031","PI656044","PI656057","PI656111","PI660557","PI660563","PI660565"
    ),
    verbose = TRUE
){
  
  ## 1) Load inventory
  if (is.data.frame(inventory_path)) {
    bms_inventory <- inventory_path
  } else if (grepl("\\.xlsx$", inventory_path, ignore.case = TRUE)) {
    bms_inventory <- readxl::read_excel(inventory_path, col_types = "text")
  } else {
    bms_inventory <- utils::read.csv(inventory_path, stringsAsFactors = FALSE)
  }

  
  ## 2) Filter to PI and known inventory: Keep only samples that have a PI number and are physically present in our lab inventory
  variant_summary_pi <- variant_summary$accessions_variants[!is.na(variant_summary$accessions_variants$PInumber), ]
  variant_summary_inventory <- variant_summary_pi[variant_summary_pi$PInumber %in% bms_inventory$PI_NAME, ]
  
  ## 4) Identify variant columns
  variant_cols <- names(variant_summary_inventory)[!(names(variant_summary_inventory) %in% metadata_columns)]
  
  ## 5) Loop per variant: Categorize available inventory into Reference vs. Alternate alleles for each priority tier
  by_variant <- list()
  
  for (v in 1:length(variant_cols)) {
    temp <- variant_summary_inventory[c(metadata_columns, variant_cols[v])]
    names(temp)[4] <- "Allele"
    
    reference_in_inventory <- temp[temp$Allele == "Reference", ]
    alternate_in_inventory <- temp[temp$Allele == "Alternate", ]
    
    reference_PCIL_recurrent <- reference_in_inventory[reference_in_inventory$PInumber %in% sorghum_PCIL_recurrent_parents, ]
    alternate_PCIL_recurrent <- alternate_in_inventory[alternate_in_inventory$PInumber %in% sorghum_PCIL_recurrent_parents, ]
    
    reference_PCIL_donor <- reference_in_inventory[reference_in_inventory$PInumber %in% sorghum_PCIL_donors, ]
    alternate_PCIL_donor <- alternate_in_inventory[alternate_in_inventory$PInumber %in% sorghum_PCIL_donors, ]
    
    reference_de_novo <- reference_in_inventory[reference_in_inventory$PInumber %in% sorghum_de_novo_references, ]
    alternate_de_novo <- alternate_in_inventory[alternate_in_inventory$PInumber %in% sorghum_de_novo_references, ]
    
    by_variant[[ variant_cols[v] ]] <- list(
      variant = variant_cols[v],
      reference_in_inventory = reference_in_inventory,
      alternate_in_inventory = alternate_in_inventory,
      reference_recurrent = reference_PCIL_recurrent,
      alternate_recurrent = alternate_PCIL_recurrent,
      reference_donor = reference_PCIL_donor,
      alternate_donor = alternate_PCIL_donor,
      reference_de_novo = reference_de_novo,
      alternate_de_novo = alternate_de_novo
    )
    
    if (verbose) {
      message(sprintf(
        "\nVariant: %s\n  %-5s | %-9s | %-6s | %-7s | %-9s\n  %-5s | %-9d | %-6d | %-7d | %-9d\n  %-5s | %-9d | %-6d | %-7d | %-9d\n",
        variant_cols[v],
        "Class", "Recurrent", "Donors", "De Novo", "Inventory",
        "ALT", nrow(alternate_PCIL_recurrent), nrow(alternate_PCIL_donor), nrow(alternate_de_novo), nrow(alternate_in_inventory),
        "REF", nrow(reference_PCIL_recurrent), nrow(reference_PCIL_donor), nrow(reference_de_novo), nrow(reference_in_inventory)
      ))
    }
  }

  
  ## 6) Select min_req per variant based on tier priority
  selected_per_variant <- list()
  
  for (v in names(by_variant)) {
    vobj <- by_variant[[v]]
    
    # Priority selection: Stack available samples in strict priority order (Recurrent > Donor > De Novo > General Inventory). 
    # Using head() naturally picks the top 'min_req' samples, automatically favoring the highest priority tier available.
    ref_stack <- rbind(vobj$reference_recurrent, vobj$reference_donor, vobj$reference_de_novo, vobj$reference_in_inventory)
    ref_stack <- ref_stack[!duplicated(ref_stack$PInumber), ]
    ref_sel <- head(ref_stack, min_req)
    if (nrow(ref_sel) > 0) ref_sel$Class <- "Reference"
    
    alt_stack <- rbind(vobj$alternate_recurrent, vobj$alternate_donor, vobj$alternate_de_novo, vobj$alternate_in_inventory)
    alt_stack <- alt_stack[!duplicated(alt_stack$PInumber), ]
    alt_sel <- head(alt_stack, min_req)
    if (nrow(alt_sel) > 0) alt_sel$Class <- "Alternate"
    
    # ADD THIS (per-variant checks)
    if (nrow(alt_sel) < min_req) {
      warning(
        paste0(
          "Variant ", v, ": minimum requirement not met for Alternate allele (needed ",
          min_req, ", found ", nrow(alt_sel), "). Consider removing this variant for validation."
        ),
        call. = FALSE
      )
    }
    
    if (nrow(ref_sel) < min_req) {
      warning(
        paste0(
          "Variant ", v, ": minimum requirement not met for Reference allele (needed ",
          min_req, ", found ", nrow(ref_sel), "). Consider removing this variant for validation."
        ),
        call. = FALSE
      )
    }  
    
    sel_v <- rbind(ref_sel, alt_sel)
    if (nrow(sel_v) > 0) sel_v$Variant <- v
    
    selected_per_variant[[v]] <- sel_v
  }
  
  sel_df <- do.call(rbind, selected_per_variant)
  rownames(sel_df) <- NULL
  
  ## Collapse variant list per PI: If an accession is chosen for multiple variants, merge them into a single row to avoid duplicate physical samples.
  unique_pinumbers_sel_df <- unique(sel_df$PInumber)
  sel_df_merged_list <- vector("list", length(unique_pinumbers_sel_df))

  for (i in seq_along(unique_pinumbers_sel_df)) {
    current_pinumber <- unique_pinumbers_sel_df[i]
    subset_data <- sel_df[sel_df$PInumber == current_pinumber, ]

    # Get the first LIB, Sample, and Class (assuming consistency within PInumber group)
    first_lib <- subset_data$LIB[1]
    first_sample <- subset_data$Sample[1]
    first_class <- subset_data$Class[1]

    # Get unique variants and collapse them
    collapsed_variants <- paste(unique(subset_data$Variant), collapse = ";")

    sel_df_merged_list[[i]] <- data.frame(
      PInumber = current_pinumber,
      LIB = first_lib,
      Sample = first_sample,
      Class = first_class,
      Variants = collapsed_variants,
      stringsAsFactors = FALSE
    )
  }
  sel_df_merged <- do.call(rbind, sel_df_merged_list)
  
  ## Artificial hets: Generate Artificial Heterozygotes by pairing one Reference and one Alternate accession for each variant.
  ref_par_subset <- sel_df[sel_df$Class == "Reference", ]
  ref_par <- ref_par_subset[order(ref_par_subset$PInumber), ]

  alt_par_subset <- sel_df[sel_df$Class == "Alternate", ]
  alt_par <- alt_par_subset[order(alt_par_subset$PInumber), ]
  
  art_het <- data.frame()
  var <- unique(sel_df$Variant)
  
  for (v in 1:length(var)) {
    tmp_r <- ref_par[ref_par$Variant==var[v],]
    tmp_a <- alt_par[alt_par$Variant==var[v],]
    # Ensure tmp_r and tmp_a are not empty before proceeding
    if (nrow(tmp_r) == 0 || nrow(tmp_a) == 0) {
      next # Skip if no reference or alternate accession found for this variant
    }
    tmp_art <- paste(tmp_r$PInumber, tmp_a$PInumber, sep = "|")
    tmp_df <- data.frame(PInumber=tmp_art, Class="Heterozygous", Variant=var[v], unique_name=NA)
    art_het <- rbind(art_het, tmp_df)
  }
  

  unique_pinumbers_art_het <- unique(art_het$PInumber)
  art_het_merged_list <- vector("list", length(unique_pinumbers_art_het))

  for (i in seq_along(unique_pinumbers_art_het)) {
    current_pinumber <- unique_pinumbers_art_het[i]
    subset_data <- art_het[art_het$PInumber == current_pinumber, ]

    # Assuming 'Class' is consistent within a PInumber group, take the first one
    first_class <- subset_data$Class[1]

    # Get unique variants and collapse them into a single string
    collapsed_variants <- paste(unique(subset_data$Variant), collapse = ";")

    art_het_merged_list[[i]] <- data.frame(
      PInumber = current_pinumber,
      unique_name = NA, # As specified in the original dplyr code
      Class = first_class,
      Variants = collapsed_variants,
      stringsAsFactors = FALSE
    )
  }
  art_het_merged <- do.call(rbind, art_het_merged_list)
  
  ## Expand the minimal list by adding the required biological and technical replicates
  sample_list_no_hets <- do.call(rbind, replicate(min_req, sel_df_merged, simplify = FALSE)) 
  sample_list_no_hets$Rep <- rep(1:min_req, each = nrow(sel_df_merged))
  sample_list_no_hets$unique_name <- paste0(sample_list_no_hets$PInumber, "_R", sample_list_no_hets$Rep)
  sample_list_no_hets <- sample_list_no_hets[c("PInumber", "Class", "Variants", "unique_name")]
  
  # Adding technical replicates
  sample_list_no_hets_techreps <- do.call(rbind, replicate(tech_rep, sample_list_no_hets, simplify = FALSE)) 
  tmp_count <- rep(1:tech_rep, each = nrow(sample_list_no_hets))
  sample_list_no_hets_techreps$unique_name <- paste0(sample_list_no_hets_techreps$unique_name, "_T", tmp_count)
  sample_list_no_hets <- sample_list_no_hets_techreps  
  
  # Replicating artificial heterozygotes
  art_het_trip <- do.call(rbind, replicate(het_tech_rep, art_het_merged, simplify = FALSE)) 
  art_het_trip$Rep <- rep(1:het_tech_rep, each = nrow(art_het_merged))
  art_het_trip$unique_name <- paste0(art_het_trip$PInumber, "_R", art_het_trip$Rep)
  art_het_trip <- art_het_trip[c("PInumber", "Class", "Variants", "unique_name")]
  
  # Generating Empty control wells
  empty_controls <- data.frame(PInumber = "Empty", Class = "NegativeControl", Variants = NA)
  empty_controls_list <- do.call(rbind, replicate(blank_reps, empty_controls, simplify = FALSE)) 
  empty_controls_list$Rep <- rep(1:blank_reps) 
  empty_controls_list$unique_name <- paste0(empty_controls_list$PInumber, "_R", empty_controls_list$Rep)
  empty_controls_list <- empty_controls_list[c("PInumber", "Class", "Variants", "unique_name")]
  
  
  sample_list_min <- rbind(sample_list_no_hets, art_het_trip, empty_controls_list) 
  
  
  return(list(
    distribution_by_variant = by_variant,
    selected_per_variant = selected_per_variant,
    minimal_sample_list = sample_list_min
  ))
}


#' Create empty plate design with blanks
#' This function generates an empty plate design for KASP and/or DArT assays, incorporating blank controls in specified positions. It calculates the number of plates needed based on the sample list and assay type, and produces a plate layout with unique identifiers for each well.
#' @param sample_list A character vector of sample names to be included in the plate design.
#' @param assay_type A character string specifying the assay type ("kasp", "dart", or "both") to  determine the number of controls per plate.
#' @param controls_in_samples A logical value indicating whether controls are included in the sample list (default: FALSE). If TRUE, the function will not add additional controls and will ignore the `controls` argument.
#' @param controls A character vector of control names to be included in the plate design if `controls_in_samples` is FALSE (default: NULL). The function will replicate these controls across plates as needed.
#' @param reps An integer specifying the number of replicates for each control when `controls_in_samples` is FALSE (default: 1). This argument is ignored if `controls_in_samples` is TRUE.
#' @return A list containing the plate layout as  a data frame with unique well identifiers, the plate matrix with control names, and a summary message about the design.
#'
#'
create_empty_plate_with_blanks <- function(
  sample_list,
  assay_type,
  controls_in_samples = FALSE,
  controls = NULL,
  reps = 1
) {
  if (controls_in_samples && !is.null(controls)) {
    warning(
      "You set `controls_in_samples = TRUE` but also provided a `controls` vector. The `controls` argument will be ignored."
    )
  }

  if (controls_in_samples && reps != 1) {
    warning(
      "You set `controls_in_samples = TRUE`, so `reps = ",
      reps,
      "` will be ignored."
    )
  }

  if (any(duplicated(sample_list))) {
    warning(
      "Sample list contains duplicate names. Plate design will proceed, but Intertek requires unique identifiers."
    )
  }

  sample_number <- length(sample_list)

  # Initialize a standard 96-well plate layout (8 rows x 12 columns)
  row_names <- LETTERS[1:8]
  col_names <- 1:12
  empty_plate <- as.data.frame(matrix(nrow = 8, ncol = 12))
  names(empty_plate) <- sprintf("%02d", col_names) # padded column names
  row.names(empty_plate) <- row_names

  # Create a long-format 1D tracking list for the wells
  empty_plate_list <- data.frame(matrix(nrow = 96, ncol = 3))
  names(empty_plate_list) <- c("Position", "Plate", "Sample")

  # Use standard A01–H12 alphanumeric well formatting
  pos <- character(96)
  counter <- 1
  for (c in 1:12) {
    for (r in 1:8) {
      col_padded <- sprintf("%02d", col_names[c])
      pos[counter] <- paste0(row_names[r], col_padded)
      counter <- counter + 1
    }
  }
  empty_plate_list$Position <- pos

  # Set mandatory control well counts based on the assay vendor specifications
  assay_control_map <- c("kasp" = 2, "dart" = 2, "both" = 3)
  assay_controls <- assay_control_map[assay_type]
  if (is.na(assay_controls)) {
    stop("Invalid 'assay_type'. Must be 'kasp', 'dart', or 'both'.")
  }

  # Calculate extra sample slots if controls need to be generated and are not in sample_list
  if (!controls_in_samples) {
    if (is.null(controls)) {
      stop(
        "You must provide a 'controls' vector if controls_in_samples is FALSE."
      )
    }
    if (!is.numeric(reps) || reps < 1) {
      stop(
        "'reps' must be a positive integer when controls are to be replicated."
      )
    }
    replicated_controls <- paste0(
      rep(controls, each = reps),
      "_R",
      rep(1:reps, times = length(controls))
    )
    sample_number <- sample_number + length(replicated_controls)
  }

  # Determine total plates required. If adding controls pushes us into a new plate, recalculate the controls needed.
  number_of_plates <- ceiling(sample_number / 96)
  total_controls <- number_of_plates * assay_controls
  number_samples_with_controls <- total_controls + sample_number
  n_plates_with_control <- ceiling(number_samples_with_controls / 96)

  if (n_plates_with_control > number_of_plates) {
    number_of_plates <- n_plates_with_control
    total_controls <- number_of_plates * assay_controls
  }

  number_samples_with_controls <- total_controls + sample_number
  empty_wells <- number_of_plates * 96 - number_samples_with_controls

  control_id <- 1:total_controls
  control_names <- paste("Blank_", control_id, sep = "")

  plate_matrix_list <- vector("list", number_of_plates)
  all_plate_lists <- vector("list", number_of_plates) # To collect temp_plate_list for efficient rbind
  control_counter <- 1

  # Build each plate iteratively, hardcoding the required blank positions (e.g., H11, H12)
  for (plates in 1:number_of_plates) {
    temp_plate <- empty_plate
    temp_plate_list <- empty_plate_list

    if (assay_type == "kasp") {
      temp_plate_list$Sample[grepl(
        "H11",
        temp_plate_list$Position
      )] <- control_names[control_counter]
      temp_plate["H", "11"] <- control_names[control_counter]
      control_counter <- control_counter + 1
      temp_plate_list$Sample[grepl(
        "H12",
        temp_plate_list$Position
      )] <- control_names[control_counter]
      temp_plate["H", "12"] <- control_names[control_counter]
      control_counter <- control_counter + 1
    }

    if (assay_type == "dart") {
      temp_plate_list$Sample[grepl(
        "G12",
        temp_plate_list$Position
      )] <- control_names[control_counter]
      temp_plate["G", "12"] <- control_names[control_counter]
      control_counter <- control_counter + 1
      temp_plate_list$Sample[grepl(
        "H12",
        temp_plate_list$Position
      )] <- control_names[control_counter]
      temp_plate["H", "12"] <- control_names[control_counter]
      control_counter <- control_counter + 1
    }

    if (assay_type == "both") {
      temp_plate_list$Sample[grepl(
        "H11",
        temp_plate_list$Position
      )] <- control_names[control_counter]
      temp_plate["H", "11"] <- control_names[control_counter]
      control_counter <- control_counter + 1
      temp_plate_list$Sample[grepl(
        "G12",
        temp_plate_list$Position
      )] <- control_names[control_counter]
      temp_plate["G", "12"] <- control_names[control_counter]
      control_counter <- control_counter + 1
      temp_plate_list$Sample[grepl(
        "H12",
        temp_plate_list$Position
      )] <- control_names[control_counter]
      temp_plate["H", "12"] <- control_names[control_counter]
      control_counter <- control_counter + 1
    }

    plate_matrix_list[[plates]] <- as.data.frame(temp_plate)
    names(plate_matrix_list)[[plates]] <- paste("Plate_", plates, sep = "")
    temp_plate_list$Plate <- paste("Plate_", plates, sep = "")
    all_plate_lists[[plates]] <- temp_plate_list # Store in list
  }

  plate_list <- do.call(rbind, all_plate_lists) # Combine all plate lists once

  # Generate unique 4-character alphabetical IDs for each well to avoid naming collisions
  n <- nrow(plate_list)
  unique_ids <- character(0)
  while (length(unique_ids) < n) {
    new_ids <- replicate(
      n - length(unique_ids),
      paste0(sample(LETTERS, 4, replace = TRUE), collapse = "")
    )
    unique_ids <- unique(c(unique_ids, new_ids))
  }

  plate_list$Unique_Well_ID <- unique_ids

  if (assay_type == "both") {
    assay_print <- "DarT and KASP"
  } else {
    assay_print <- assay_type
  }

  # Create 2D visual matrix representations of each plate for easy visual inspection of the layout
  plate_matrix_unique_id <- list()
  for (plate_name in unique(plate_list$Plate)) {
    tmp_df <- plate_list[plate_list$Plate == plate_name, ]
    tmp_df$Row <- gsub("[0-9]+", "", tmp_df$Position)
    tmp_df$Col <- sprintf(
      "%02d",
      as.numeric(gsub("^[A-Z]", "", tmp_df$Position))
    )

    mat <- matrix(
      NA,
      nrow = 8,
      ncol = 12,
      dimnames = list(LETTERS[1:8], sprintf("%02d", 1:12))
    )

    for (i in seq_len(nrow(tmp_df))) {
      mat[tmp_df$Row[i], tmp_df$Col[i]] <- tmp_df$Unique_Well_ID[i]
    }
    plate_matrix_unique_id[[plate_name]] <- as.data.frame(mat)
  }

  message(sprintf(
    "\n----------------------------------------\nPlate Design Summary\n----------------------------------------\n%-20s : %s\n%-20s : %d\n%-20s : %d\n%-20s : %d\n%-20s : %d\n%-20s : %d\n----------------------------------------\n* Add %d samples to your sample list to fill all wells.\n",
    "Assay Type",
    assay_print,
    "Samples",
    sample_number,
    "Plates Required",
    number_of_plates,
    "Total Blanks",
    total_controls,
    "Blanks per Plate",
    assay_controls,
    "Empty Wells",
    empty_wells,
    empty_wells
  ))

  return(list(
    plate_list = plate_list,
    plate_matrix_list = plate_matrix_list,
    plate_matrix_unique_id = plate_matrix_unique_id
  ))
}


#' Fill plate design with samples
#' This function takes the empty plate design generated by `create_empty_plate_with_blanks` and fills it with the provided sample list. It can randomize the placement of samples across the plates or fill them in a sequential manner, depending on the user's choice. The function updates the plate layout and provides a summary of the design after filling the samples.
#' @param plate_lists A list object containing the plate layout components (plate_list, plate_matrix_list, plate_matrix_unique_id) as generated by `create_empty_plate_with_blanks`.
#' @param sample_list A character vector of sample names to be filled into the plate design. The function will check for the number of samples against the available empty wells and will stop with an error if there are not enough empty wells to accommodate all samples.
#' @param randomization A logical value indicating whether to randomize the placement of samples across the plates (default: TRUE). If FALSE, samples will be filled sequentially in the order they appear in the sample list.
#' @return A list containing the updated `plate_list` (long format) and `plate_matrix_list` (2D matrix format) with filled sample information.
#'
#' 
fill_plate_samples <- function(
    plate_lists,              
    sample_list,             
    randomization = TRUE     
) {
  ### Filling the empty layout plates with the finalized samples
  # Isolate the list.
  plate_list <- plate_lists[['plate_list']]

  if (randomization) {
    # Randomly assign samples to empty wells to avoid batch/edge effects during plate reading
    empty_plate_list <- plate_list[is.na(plate_list$Sample), ]
    
    if (length(sample_list) > nrow(empty_plate_list)) {
      stop(paste("Not enough empty wells:",
                 "you have", length(sample_list), "samples but only",
                 nrow(empty_plate_list), "empty wells"))
    }
    
    sample_wells <- empty_plate_list[sample(nrow(empty_plate_list), size = length(sample_list)), ]
    sample_wells$Sample <- sample(sample_list, length(sample_list), replace = FALSE)
    
    plate_list$Sample[plate_list$Unique_Well_ID %in% sample_wells$Unique_Well_ID] <-
      sample_wells$Sample[match(
        plate_list$Unique_Well_ID[plate_list$Unique_Well_ID %in% sample_wells$Unique_Well_ID],
        sample_wells$Unique_Well_ID
      )]
    
  } else {
    # Fill the plate sequentially (A01, B01, etc.) for easier manual pipetting/tracking
    empty_plate_list <- plate_list[is.na(plate_list$Sample), ]
    sample_wells <- empty_plate_list[1:length(sample_list), ]
    sample_wells$Sample <- sample_list
    
    plate_list$Sample[plate_list$Unique_Well_ID %in% sample_wells$Unique_Well_ID] <-
      sample_wells$Sample[match(
        plate_list$Unique_Well_ID[plate_list$Unique_Well_ID %in% sample_wells$Unique_Well_ID],
        sample_wells$Unique_Well_ID
      )]
  }
  
  # Push the 1D updated assignments back into the 2D visual plate matrix
  plate_matrix_list <- plate_lists[['plate_matrix_list']] # Re-extract as it's modified
  for (plates in 1:length(plate_matrix_list)) {
    temp_name <- names(plate_matrix_list)[plates]
    tmp_plate_matrix <- plate_matrix_list[[temp_name]]
    names(tmp_plate_matrix) <- sprintf("%02d", 1:12) # Ensure consistent column naming

    tmp_plate_all <- plate_list[plate_list$Plate == temp_name, ] # Filter plate_list for the current plate

    # Get all filled wells (non-NA samples) for the current plate
    tmp_plate_samples <- tmp_plate_all[!is.na(tmp_plate_all$Sample), ]
    
    # Extract and format column names to match the matrix's column names (e.g., "01", "02")
    tmp_plate_samples$Column <- sprintf("%02d", as.numeric(gsub(
      "^[A-Z]",
      "",
      tmp_plate_samples$Position
    )))
    tmp_plate_samples$Row <- gsub("[0-9]+", "", tmp_plate_samples$Position)
    
    # Convert to matrix for efficient vectorized assignment, then back to data.frame
    tmp_plate_matrix_as_mat <- as.matrix(tmp_plate_matrix)
    tmp_plate_matrix_as_mat[cbind(tmp_plate_samples$Row, tmp_plate_samples$Column)] <- tmp_plate_samples$Sample
    plate_matrix_list[[temp_name]] <- as.data.frame(tmp_plate_matrix_as_mat)
  }
  
  
  empty_wells <- sum(is.na(plate_list$Sample))
  message(paste("Your testing plate has", length(sample_list), "samples.",
                "This design has", empty_wells, "empty wells.",
                "If you haven't yet filled the controls yet, make sure they fit in the plate.",
                "If you have already, add", empty_wells, "samples to fill all wells."))
  
  return(list(
    plate_list = plate_list,
    plate_matrix_list = plate_matrix_list
  ))
}



#' @param selection A character vector of variant IDs (e.g., "SNP_Chr03_79037855") to summarize.
#' @param output_dir The directory where output files will be written (default: current directory).
#' @param vcf A character string or vector specifying the path(s) to the VCF file(s).
#' @param write_files A logical value indicating whether to write output tables to files (default: TRUE).
#' @param metadata_csv A character string specifying the path to a CSV file containing sample metadata, including 'LIB' and 'PInumber' columns.
#' @return A list containing `accessions_variants` (a data frame of accessions and their allele classifications)
#'   and `variant_summary` (a data frame summarizing allele frequencies and PCV annotations).
#'
#' @importFrom utils read.csv
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR rbind2
#'
summarize_variant_alleles <- function(
  selection,
  output_dir = ".",
  vcf = "", # directory where the VCFs live
  write_files = TRUE,
  metadata_csv
) {
  # --- Inputs & paths ---
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # --- Load metadata ---
  metadat = read.csv(metadata_csv)
  # --- Load VCFs (Supports single or multiple overlapping VCF files) ---
  if (length(vcf) == 1) {
    # Case 1: only one VCF provided
    snp_vcf <- vcfR::read.vcfR(vcf)
  } else {
    # Case 2: multiple VCFs provided → try merging
    vcfs <- lapply(vcf, vcfR::read.vcfR)

    # Check if all sample columns (excluding FORMAT) are identical
    sample_cols <- lapply(vcfs, function(v) colnames(v@gt)[-1])
    same_samples <- all(vapply(
      sample_cols[-1],
      function(x) identical(x, sample_cols[[1]]),
      logical(1)
    ))

    # If samples differ → stop with a clear message
    if (!same_samples) {
      stop(
        "The VCF files cannot be merged because they do not have identical sample sets or sample order.\n",
        "Please ensure all VCFs contain the same samples in the same order before merging, ",
        "or process each variant separately with its corresponding VCF."
      )
    }

    # If samples match → merge them
    snp_vcf <- Reduce(vcfR::rbind2, vcfs)

    # Optional: warn if duplicate IDs exist
    fix_ids <- as.data.frame(snp_vcf@fix, stringsAsFactors = FALSE)$ID
    dups <- unique(fix_ids[duplicated(fix_ids)])
    if (length(dups) > 0) {
      warning(sprintf(
        "Merged VCFs contain duplicate variant IDs (n=%d). First occurrence will be used downstream.",
        length(dups)
      ))
    }
  }

  # --- Convert standard VCF genotype strings into numeric classes (0=Ref, 1=Het, 2=Alt, NA=Missing) ---
  dat <- snp_vcf@gt
  for (i in seq_len(ncol(dat))) {
    dat[dat[, i] %in% c("1|1", "1/1"), i] <- 2
    dat[dat[, i] %in% c("0|1", "1|0", "0/1", "1/0"), i] <- 1
    dat[dat[, i] %in% c("0|0", "0/0"), i] <- 0
    dat[dat[, i] %in% c(".|.", "./."), i] <- NA
    dat[, i] <- suppressWarnings(as.numeric(dat[, i]))
  }
  dat[, 1] <- "GT" # keep GT header
  snp_vcf@gt <- dat

  # Helper function to process the subsetting and aggregation for a single variant ID
  process_one <- function(sel_id) {
    # --- Build genotype matrix for the selected variant ---
    geno <- as.data.frame(snp_vcf@gt, stringsAsFactors = FALSE)
    geno <- geno[, -1, drop = FALSE] # drop GT header col
    fix <- as.data.frame(snp_vcf@fix, stringsAsFactors = FALSE)

    # Locate the specific variant in the VCF by its ID (avoiding partial matches)
    idx <- which(fix$ID == sel_id)
    if (length(idx) == 0) {
      stop(sprintf("Selected variant '%s' not found in VCF.", sel_id))
    }
    if (length(idx) > 1) {
      warning(sprintf(
        "Selected variant '%s' found %d times; using first.",
        sel_id,
        length(idx)
      ))
    }
    idx <- idx[1]

    # REF/ALT from VCF (handle multi-ALT by taking first)
    REF <- fix$REF[idx]
    ALT <- fix$ALT[idx]
    if (grepl(",", ALT, fixed = TRUE)) {
      ALT_split <- strsplit(ALT, ",", fixed = TRUE)[[1]]
      warning(sprintf(
        "Multi-ALT site detected (%s). Using first ALT: %s",
        ALT,
        ALT_split[1]
      ))
      ALT <- ALT_split[1]
    }

    selection_geno <- geno[idx, , drop = FALSE]
    t_selection_geno <- as.data.frame(
      t(selection_geno),
      stringsAsFactors = FALSE
    )

    # --- Haplotype summary: Tally up the frequencies of Reference, Heterozygous, and Alternate alleles ---
    allele_values <- t_selection_geno[, 1]
    allele_levels <- c("0", "1", "2")
    allele_labels <- c("Reference", "Heterozygous", "Alternate")

    haplotype_summary <- table(factor(allele_values, levels = allele_levels))
    haplotype_summary <- as.data.frame(haplotype_summary)
    names(haplotype_summary) <- c("Allele", "Frequency")
    haplotype_summary$Allele <- factor(
      haplotype_summary$Allele,
      levels = allele_levels,
      labels = allele_labels
    )

    # Rename Frequency to be unique when merging
    names(haplotype_summary)[names(haplotype_summary) == "Frequency"] <- paste(
      "Frequency",
      sel_id,
      sep = "_"
    )

    # PCV annotation
    haplotype_summary$PCV <- NA
    haplotype_summary$PCV[haplotype_summary$Allele == "Reference"] <- REF
    haplotype_summary$PCV[haplotype_summary$Allele == "Heterozygous"] <- paste(
      REF,
      "/",
      ALT,
      sep = ""
    )
    haplotype_summary$PCV[haplotype_summary$Allele == "Alternate"] <- ALT

    #  Rename "PCV" column to the variant ID
    names(haplotype_summary)[names(haplotype_summary) == "PCV"] <- sel_id

    # --- Accessions (allele class + LIB + Sample) ---
    accessions_variants <- t_selection_geno
    names(accessions_variants)[1] <- sel_id
    accessions_variants$LIB <- row.names(accessions_variants)
    accessions_variants[, 1] <- factor(
      accessions_variants[, 1],
      levels = allele_levels,
      labels = allele_labels
    )

    # Merge metadata to get "Sample" (and keep LIB)
    accessions_variants <- merge(
      accessions_variants,
      metadat,
      by.x = "LIB",
      by.y = "LIB",
      all.x = TRUE
    )
    # Keep: [AlleleClass, LIB, Sample]
    keep_cols <- c(sel_id, "LIB", "Sample", "PInumber")
    keep_cols <- keep_cols[keep_cols %in% colnames(accessions_variants)]
    accessions_variants <- accessions_variants[, keep_cols, drop = FALSE]

    # --- Filenames & write ---
    file_accession_hap <- file.path(
      output_dir,
      paste(sel_id, "sample_hap_summary_table.txt", sep = "_")
    )
    file_haplotype_freq <- file.path(
      output_dir,
      paste(sel_id, "haplo_freq_summary_table.txt", sep = "_")
    )

    if (isTRUE(write_files)) {
      write.table(
        accessions_variants,
        file_accession_hap,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      write.table(
        haplotype_summary,
        file_haplotype_freq,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }

    # --- Return both tables ---
    return(list(
      accessions_variants = accessions_variants,
      variant_summary = haplotype_summary
    ))
  }

  # Process all selected variant IDs using the defined helper function
  selections <- as.character(selection)
  if (length(selections) == 1) {
    return(process_one(selections[1]))
  } else {
    out_list <- lapply(selections, process_one)
    names(out_list) <- selections

    # --- Combine haplotype summaries (using renamed PCV column) ---
    haplo_combined <- NULL
    for (sel in selections) {
      freq_col <- paste("Frequency", sel, sep = "_")
      tmp <- out_list[[sel]]$variant_summary[, c("Allele", freq_col, sel)]
      if (is.null(haplo_combined)) {
        haplo_combined <- tmp
      } else {
        haplo_combined <- merge(
          haplo_combined,
          tmp,
          by = "Allele",
          all.x = TRUE
        )
      }
    }

    # --- Combine accession alleles (allelic class) and keep metadata ---
    first <- selections[1]
    # Start with LIB + Sample + PInumber from the first variant
    accessions_combined <- out_list[[first]]$accessions_variants[,
      c("LIB", "Sample", "PInumber"),
      drop = FALSE
    ]
    # Add each variant’s allele classification
    for (sel in selections) {
      tmp <- out_list[[sel]]$accessions_variants[, c("LIB", sel), drop = FALSE]
      accessions_combined <- merge(
        accessions_combined,
        tmp,
        by = "LIB",
        all.x = TRUE
      )
    }

    return(list(
      accessions_variants = accessions_combined,
      variant_summary = haplo_combined
    ))
  }
}


#' Optimize minimal sample list parameters for perfect plate design
#' 
#' This function automates the search for the best combination of parameters 
#' (`min_req`, `tech_rep`, `het_tech_rep`, and `blank_reps`) to perfectly fill 
#' a testing plate with exactly 0 empty wells, avoiding manual guesswork.
#' 
#' @param variant_summary A list containing variant summary data from `summarize_variant_alleles`.
#' @param inventory_path A URL or file path pointing to the sample inventory, or a loaded data frame.
#' @param assay_type A character string specifying the assay type ("kasp", "dart", or "both").
#' @param min_req_range Numeric vector of `min_req` values to explore (default: 3:5).
#' @param tech_rep_range Numeric vector of `tech_rep` values to explore (default: 1:3).
#' @param het_tech_rep_range Numeric vector of `het_tech_rep` values to explore (default: 3:6).
#' @return A list identical to the output of `build_min_sample_list`, but optimized to perfectly fill plates.
#' @importFrom readxl read_excel
#' 
optimize_plate_params <- function(
    variant_summary,
    inventory_path,
    assay_type = "kasp",
    min_req_range = 3:5,
    tech_rep_range = 1:3,
    het_tech_rep_range = 3:6
) {
  
  # Determine mandatory control well counts based on the assay vendor specifications
  assay_control_map <- c("kasp" = 2, "dart" = 2, "both" = 3)
  assay_controls <- assay_control_map[assay_type]
  if (is.na(assay_controls)) {
    stop("Invalid 'assay_type'. Must be 'kasp', 'dart', or 'both'.")
  }
  usable_per_plate <- 96 - assay_controls
  
  # Load inventory (similar logic to build_min_sample_list)
  if (is.data.frame(inventory_path)) {
    bms_inventory <- inventory_path
  } else if (grepl("\\.xlsx$", inventory_path, ignore.case = TRUE)) {
    bms_inventory <- readxl::read_excel(inventory_path, col_types = "text")
  } else {
    bms_inventory <- utils::read.csv(inventory_path, stringsAsFactors = FALSE)
  }
  
  # Generate all possible combinations of the provided parameters to test
  grid <- expand.grid(min_req = min_req_range, tech_rep = tech_rep_range, het_tech_rep = het_tech_rep_range)
  
  best_diff <- Inf
  best_params <- NULL
  best_blank_reps <- NULL
  
  message("Scanning parameter combinations to find the optimal plate fill...\n")
  
  # Iterate through every combination, quietly building the sample list to evaluate how many wells it utilizes
  for (i in 1:nrow(grid)) {
    p <- grid[i, ]
    temp_rec <- suppressMessages(suppressWarnings(build_min_sample_list(
      variant_summary, bms_inventory, min_req = p$min_req, tech_rep = p$tech_rep, het_tech_rep = p$het_tech_rep, blank_reps = 0
    )))
    
    base_samples <- nrow(temp_rec$minimal_sample_list)
    # Calculate how many wells would be left empty on the final plate
    remainder <- base_samples %% usable_per_plate
    required_blanks <- if (remainder == 0) 0 else usable_per_plate - remainder
    
    # Keep the parameters that produce the least amount of wasted plate space
    if (required_blanks < best_diff) {
      best_diff <- required_blanks
      best_params <- p
      best_blank_reps <- required_blanks
    }
    if (required_blanks == 0) break
  }
  
  # Re-run and return the list built with the discovered optimal parameters
  final_rec <- suppressMessages(suppressWarnings(build_min_sample_list(
    variant_summary, bms_inventory, min_req = best_params$min_req, tech_rep = best_params$tech_rep, het_tech_rep = best_params$het_tech_rep, blank_reps = best_blank_reps
  )))
  
  message(paste0(
    "\n",
    "==================================================\n",
    "  OPTIMIZATION COMPLETE\n",
    "--------------------------------------------------\n",
    sprintf("%-25s | %s\n", "Parameter", "Optimal Value"),
    "--------------------------|-----------------------\n",
    sprintf("%-25s | %d\n", "Min. Requirements", best_params$min_req),
    sprintf("%-25s | %d\n", "Tech. Replicates", best_params$tech_rep),
    sprintf("%-25s | %d\n", "Het. Tech. Replicates", best_params$het_tech_rep),
    sprintf("%-25s | %d\n", "Blank Replicates", best_blank_reps),
    "--------------------------|-----------------------\n",
    sprintf("%-25s | %s\n", "Plate Efficiency", "100% (0 empty wells)"),
    "==================================================\n\n"
  ))
  
  return(final_rec)
}



#' Export Plate Layouts to CSV Files
#'
#' This function takes the filled plate design and sample recommendation objects and exports the plate list (long format)
#' and the first plate's matrix (2D format) to CSV files. The filenames are dynamically generated based on the variants
#' included in the `sample_recommendation_obj`.
#' @param plates_filled_obj A list object containing the filled plate layouts, typically the output from `fill_plate_samples()`.
#' @param sample_recommendation_obj A list object containing sample recommendations, typically the output from `build_min_sample_list()`.
#' @return This function does not return an R object. It writes two CSV files to the specified `output_dir`:
#'   - `[base_name_file]_plate_design_list.csv`: Contains the plate layout in a long (list) format.
#'   - `[base_name_file]_plate_design_matrix.csv`: Contains the layout of the first plate in a 2D matrix format.
#' 
export_plate_layouts_to_csv <- function(plates_filled_obj, sample_recommendation_obj) {
  # Extract variant names from the sample recommendation for dynamic file naming.
  variants_in_plate <- names(sample_recommendation_obj[["selected_per_variant"]])
  variants_tag <- paste(variants_in_plate, collapse = "_")

  # Construct a base filename using the variant tag.
  base_name_file <- paste0(
    "Validation_plate_for_",
    variants_tag,
    "_for_KASP_development"
  )

  # Define the full filenames for the plate list and plate matrix CSVs.
  list_name <- paste0(base_name_file, "_plate_design_list.csv")
  matrix_name <- paste0(base_name_file, "_plate_design_matrix.csv")

  # Write the plate list (long format) to a CSV file.
  write.csv(
    plates_filled_obj[["plate_list"]],
    file = list_name,
    quote = FALSE,
    row.names = FALSE
  )

  # Write the  plate's matrix (2D format) to a CSV file.
  # Get the names of all the plates generated (e.g., "Plate_1", "Plate_2")
plate_names <- names(plates_filled[["plate_matrix_list"]])

# Loop through each plate and save it as a distinct CSV
for (p_name in plate_names) {
  # Create a unique file name for each plate
  current_file_name <- paste0(base_name_file, "_", p_name, ".csv")

  # Write the specific plate matrix to the file
   write.csv(
    plates_filled_obj[["plate_matrix_list"]][[p_name]],
     file = current_file_name,
     quote = FALSE,
    row.names = TRUE
   )
}
}
