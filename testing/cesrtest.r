library(cancereffectsizeR)
library(data.table)

# Download TCGA lung adenocarcinoma (LUAD) somatic variant data.
tcga_maf_file <- "TCGA-LUAD.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "LUAD", filename = "TCGA-LUAD.maf.gz")
}

# Prepare data
maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")

# Create cancereffectsizeR analysis and load data
cesa <- CESAnalysis(refset = "ces.refset.hg38")
cesa <- load_maf(cesa = cesa, maf = maf)

# Infer trinculeotide-context-specific relative rates of SNV mutation from
# a mutational signature analysis (leaving out signatures not found in LUAD)
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "LUAD", treatment_naive = TRUE)
cesa <- trinuc_mutation_rates(
  cesa = cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
  signature_exclusions = signature_exclusions
)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$lung)

# Infer scaled selection coefficients under the default model of clonal selection.
# By default, inference is restricted to recurrent mutations.
cesa <- ces_variant(cesa, run_name = "example")

# Visualize top-effect variants.
plot_effects(effects = cesa$selection$example, color_by = "#DB382D", topn = 20)

# Attribute effects to mutational signatures
mut_effects <- mutational_signature_effects(cesa, cesa$selection$example)

# Plot a comparison of how signatures contribute to mutation vs. selection
plot_signature_effects(mut_effects, viridis_option = "F", num_sig_groups = 5)

# See the full tutorial for more details and a broader view of functionality!