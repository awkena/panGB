# ==============================================================================
# KASP Marker Validation Plate Design Workflow
#
# Description:
# This script processes variant data (from VCFs) and physical seed inventory
# to design 96-well validation plates for KASP markers. It identifies REF/ALT
# alleles, selects optimal validation accessions, and generates ready-to-use
# randomized plate layouts.
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1: Environment Setup
# ------------------------------------------------------------------------------

# It is highly recommended to use RStudio Projects or the 'here' package
# instead of absolute paths for better reproducibility across different computers.
# setwd('/Users/israeltawiahtetteh/Desktop/PanGB pipeline scripts')

# Load custom plate design and variant summarization functions
source("Plate_design_scripts.R")


# ==============================================================================
# PART A: SINGLE VARIANT WORKFLOW
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 2: Extract Variant Data & Select Core Samples
# ------------------------------------------------------------------------------

# 2.1 Extract variant data and classify alleles
# Parse the VCF file for a specific target SNP to identify which accessions
# carry the Reference (REF) or Alternate (ALT) alleles.
variant_summary <- summarize_variant_alleles(
  selection = "SNP_Chr03_79037855",
  vcf = "Sobic.003G421300_SNP_INDEL_snpeff_combined.vcf.gz",
  output_dir = "./",
  write_files = FALSE,
  install_missing = TRUE
)

# 2.2 Select core validation samples
# Build a minimal list of samples from the lab's seed inventory that ensures
# both REF and ALT alleles are adequately represented.
sample_recommendation <- build_min_sample_list(
  variant_summary = variant_summary,
  inventory_path = "Seed Inventory for BMS_Clara.xlsx"
)

# --- Optional Inspection Block ---
# View(variant_summary$accessions_variants)
# View(variant_summary$variant_summary)
# View(sample_recommendation$distribution_by_variant$SNP_Chr03_79037855$reference_in_inventory)
# View(sample_recommendation$selected_per_variant$SNP_Chr03_79037855)
# View(sample_recommendation$minimal_sample_list)
# -------------------------------

# ------------------------------------------------------------------------------
# STEP 3: Initial Plate Scaffolding & Export (Single Variant)
# ------------------------------------------------------------------------------

# 3.1 Initial plate scaffolding
# Generate an empty 96-well plate design containing mandatory vendor blanks
# (e.g., H11, H12 for KASP) and allocate slots for our samples.
my_sample_list <- sample_recommendation$minimal_sample_list$unique_name

plate_design <- create_empty_plate_with_blanks(
  sample_list = my_sample_list,
  assay_type = "kasp",
  controls_in_samples = TRUE
)

plates_filled <- fill_plate_samples(
  plate_lists = plate_design,
  sample_list = my_sample_list,
  randomization = TRUE
)

# 3.2 Export the final plate layouts to CSV files for Intertek
export_plate_layouts_to_csv(
  plates_filled_obj = plates_filled,
  sample_recommendation_obj = sample_recommendation
 )

# NOTE: End here if you do not need to fill all 96 wells on the plate.

# ------------------------------------------------------------------------------
# STEP 4: Plate Optimization & Manual Adjustments (Fill all 96 wells)
# ------------------------------------------------------------------------------

# Option A: Auto-optimize parameters to achieve exactly 0 empty wells
sample_recommendation_opt <- optimize_plate_params(
  variant_summary = variant_summary,
  inventory_path = "Seed Inventory for BMS_Clara.xlsx",
  assay_type = "kasp"
)

my_sample_list_opt <- sample_recommendation_opt$minimal_sample_list$unique_name

plate_design_opt <- create_empty_plate_with_blanks(
  sample_list = my_sample_list_opt,
  assay_type = "kasp",
  controls_in_samples = TRUE
)

plates_filled_opt <- fill_plate_samples(
  plate_lists = plate_design_opt,
  sample_list = my_sample_list_opt,
  randomization = TRUE
)


# Export plates as csv files.
export_plate_layouts_to_csv(
  plates_filled_obj = plates_filled_opt,
  sample_recommendation_obj = sample_recommendation_opt )


# Option B: Manually force specific accessions into the plate
# Format strings for manual additions (e.g., PI514353_R1_T1, PI514353_R1_T2)
tech_rep <- 2
bio_rep <- 3

manual_sample <- rep("PI514353", bio_rep)
manual_sample <- paste0(manual_sample, "_R", seq_len(bio_rep))
manual_sample <- rep(manual_sample, tech_rep)
manual_sample <- paste0(
  manual_sample,
  "_T",
  rep(seq_len(tech_rep), each = bio_rep)
)

# Append the manual samples to the newly optimized list with extra blanks
sample_recommendation_manual <- build_min_sample_list(
  variant_summary = variant_summary,
  inventory_path = "Seed Inventory for BMS_Clara.xlsx",
  tech_rep = 2,
  min_req = 4,
  het_tech_rep = 5,
  blank_reps = 4
)

my_sample_list_manual <- c(
  sample_recommendation_manual$minimal_sample_list$unique_name,
  manual_sample
)

plates_design_manual_complete <- create_empty_plate_with_blanks(
  sample_list = my_sample_list_manual,
  assay_type = "kasp",
  controls_in_samples = TRUE
)

plates_filled_manual_complete <- fill_plate_samples(
  plate_lists = plates_design_manual_complete,
  sample_list = my_sample_list_manual,
  randomization = TRUE
)

# Export result as an excel file.
export_plate_layouts_to_csv(
  plates_filled_obj = ,
  sample_recommendation_obj = )


# ==============================================================================
# PART B: MULTI-VARIANT (MULTIPLEX) WORKFLOW
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1: Summarize & Filter Multiple Variants
# ------------------------------------------------------------------------------

# 1.1 Extract and summarize allele data for multiple variants simultaneously
variant_summary_mult <- summarize_variant_alleles(
  selection = c(
    "INDEL_Chr03_79037889",
    "SNP_Chr03_79037855",
    "SNP_Chr03_79037944"
  ),
  vcf = "Sobic.003G421300_SNP_INDEL_snpeff_combined.vcf.gz",
  output_dir = "./",
  write_files = FALSE,
  install_missing = TRUE
)

# 1.2 Build a single validation sample list that satisfies ALL variants
sample_recommendation_mult <- build_min_sample_list(
  variant_summary = variant_summary_mult,
  inventory_path = "Seed Inventory for BMS_Clara.xlsx"
)


# ------------------------------------------------------------------------------
# STEP 1.3: Multi-Variant Plate Scaffolding & Optimization
# ------------------------------------------------------------------------------

# Reduce technical replicates to force the design into a single plate
sample_recommendation_mult_opt1 <- optimize_plate_params(
  variant_summary = variant_summary_mult,
  inventory_path = "Seed Inventory for BMS_Clara.xlsx",
  assay_type = "kasp"
)

sample_list_mult_opt1 <- sample_recommendation_mult_opt1$minimal_sample_list$unique_name

plate_design_mult_opt1 <- create_empty_plate_with_blanks(
  sample_list = sample_list_mult_opt1,
  assay_type = "kasp",
  controls_in_samples = TRUE
)

plate_filled_mult_opt1 <- fill_plate_samples(
  plate_lists = plate_design_mult_opt1,
  sample_list = sample_list_mult_opt1,
  randomization = TRUE
)

# Export Multi-Variant Plate Layouts
export_plate_layouts_to_csv(
  plates_filled_obj = plate_filled_mult_opt1,
  sample_recommendation_obj = sample_recommendation_mult_opt1
)
