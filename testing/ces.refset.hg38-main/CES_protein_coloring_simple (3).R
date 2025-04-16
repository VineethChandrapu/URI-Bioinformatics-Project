remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR", dependencies = TRUE, repos = BiocManager::repositories(),force = T)
remotes::install_github("nvelden/NGLVieweR")
install.packages("bio3d", dependencies=TRUE)
library(cancereffectsizeR)
library(data.table)
library(NGLVieweR)
library(bio3d)

####THIS CHUNK JUST RUNS CES TUTORIAL ON LOCAL DATA####
tcga_maf_file <- "TCGA-BRCA.maf.gz"
if (!file.exists(tcga_maf_file2)) {get_TCGA_project_MAF(project = "BRCA", filename = tcga_maf_file)}
tcga_clinical <- fread(system.file("tutorial/TCGA_BRCA_clinical.txt", package = "cancereffectsizeR"))
setnames(tcga_clinical, "patient_id", "Unique_Patient_Identifier")
tcga_maf_file<-read.csv("~/../Downloads/brca_tcga_pan_can_atlas_2018/data_mutations.txt",header = T,comment.char = "#",sep = "\t")
tcga_maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38",chain_file = "~/../Desktop/hg19ToHg38.over.chain")
tgs_maf_file <- system.file("tutorial/metastatic_breast_2021_hg38.maf", package = "cancereffectsizeR")
tgs_maf <- preload_maf(maf = tgs_maf_file, refset = "ces.refset.hg38")
cesa <- CESAnalysis(refset = "ces.refset.hg38")
cesa <- load_maf(cesa = cesa, maf = tcga_maf_file, maf_name = "BRCA")
top_tgs_genes <- c("TP53", "PIK3CA", "ESR1", "CDH1", "GATA3", "KMT2C","MAP3K1", "AKT1", "ARID1A", "FOXA1", "TBX3", "PTEN")
tgs_coverage <- ces.refset.hg38$gr_genes[ces.refset.hg38$gr_genes$gene %in% top_tgs_genes]
tgs_maf$pM <- "M1"
cesa <- load_maf(cesa,maf = tgs_maf, sample_data_cols = "pM", maf_name = "MBC",coverage = "targeted", covered_regions = tgs_coverage,covered_regions_name = "top_genes", covered_regions_padding = 10)


signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "BRCA", treatment_naive = TRUE)

cesa <- trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
                              signature_exclusions = signature_exclusions
)
cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$breast)
cesa <- ces_variant(cesa = cesa, run_name = "recurrents")
###### END TUTORIAL RUNNING 


#### YOU NOW HAVE A CESA OBJECT CALLED cesa.
#### AT THE MOMENT YOU ALSO NEED TO MANUALLY DOWNLOAD A PDB FOR YOUR GENE
#### THIS IS BEING WORKED ON BUT ISN'T YET DONE


#####INPUTS YOU CAN EDIT HERE
#PDB <= hopefully will have code that does this soon
this_pdb_file<-"~/../Downloads/AF-P03372-F1-model_v4.pdb"
#GENE YOU ARE PLOTTING AS NAMED IN cesa
gene_to_plot<-"ESR1"
#NOT YET IMPLEMENTED; WOULD GIVE OPTION OF RENDERING EITHER NORMAL OR MUTANT FORM
#norm_or_mut<-"mut"
#INCLUDE LABELS?
doLabels<-TRUE
#COLOR BY LOG INTENSITY OR NO?
logScale<-T
#COLOR OF NON-INTERESTING RESIDUES
default_col<-c(211, 211, 211) #grey
#COLORS IN RGB FOR MIN->MID-MAX
#YOU CAN SET mid_ramp_color<-null if you don't want to specify a mid color
max_ramp_color<-c(220, 20, 60) #dark red
min_ramp_color<-c(255, 245, 238) #v light pink
mid_ramp_color<-c(255,140,0) #orange

####END OF INPUT
#### AFTER RUNNING THE CODE BELOW
#### YOU'LL GET AN _NGL.R FILE WRITTEN TO YOUR WORKING DIR
#### OPEN IT UP AND RUN IT TO RENDER THE PROTEIN

this_pdb<-read.pdb(this_pdb_file)
fseq<-(tolower(as.character(this_pdb$seqres)))
fseq<-paste(toupper(substr(fseq, 1, 1)), substr(fseq, 2, nchar(fseq)), sep="")
fseq<-seqinr::a(fseq)

col_table<-data.frame(fseq)
col_table$mut_to<-"none"
col_table$sel<-0
col_table$color<-"default"
all_hits<-cesa$selection$recurrents[cesa$selection$recurrents$gene==gene_to_plot,]
all_hits<-all_hits[order(all_hits$selection_intensity),] #small to biggest, so small gets overwritten purposely

for(i in 1:nrow(all_hits)){
  this_hit<-all_hits[i,]
  short_name<-gsub(pattern = paste0(this_hit$gene," "),replacement = "",this_hit$variant_name)
  print(short_name)
  this_from<-substr(short_name, 1, 1)
  this_to<-substr(short_name,nchar(short_name),nchar(short_name))
  this_pos<-as.numeric(substr(short_name,2,nchar(short_name)-1))
  this_sel<-this_hit$selection_intensity
  
  if(col_table$fseq[this_pos]!=this_from){
    warning("CES AA does not match structure AA at this position")
    warning("if this happens often, the protein and CES are not using the same reference")
  }
  col_table$mut_to[this_pos]<-this_to
  col_table$sel[this_pos]<-this_sel
}

col_table$log_sel<-col_table$sel
col_table$log_sel[col_table$log_sel>0]<-log(col_table$log_sel[col_table$log_sel>0])

maxHex <- rgb(max_ramp_color[1], max_ramp_color[2], max_ramp_color[3], maxColorValue = 255)
minHex <- rgb(min_ramp_color[1], min_ramp_color[2], min_ramp_color[3], maxColorValue = 255)

if (!is.null(mid_ramp_color)) {
  midHex <-
    rgb(mid_ramp_color[1], mid_ramp_color[2], mid_ramp_color[3], maxColorValue = 255)
}

defHex <- rgb(default_col[1], default_col[2], default_col[1], maxColorValue = 255)


numberOfColors <- sum(col_table$sel>0)
rel_ents<-col_table[col_table$sel>0,]


#get color, as normal
values <- rel_ents$sel
if(logScale){values<-rel_ents$log_sel}
ii <-
  cut(
    values,
    breaks = seq(min(values), max(values), len = numberOfColors),
    include.lowest = TRUE
  )
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
if (is.null(midHex)) {
  myColors <-
    colorRampPalette(c(minHex, maxHex))((numberOfColors - 1))[ii]
}
if (!is.null(midHex)) {
  myColors <-
   colorRampPalette(c(minHex, midHex, maxHex))((numberOfColors - 1))[ii]
}


rel_ents$color <- myColors
if(doLabels){rel_ents$label<-paste0(rel_ents$fseq,rownames(rel_ents),rel_ents$mut_to)}

###Build NGL_R file for viewing

NGL_file<-c()
NGL_file<-c(NGL_file,"NGLVieweR(this_pdb_file) %>%")
NGL_file<-c(NGL_file,'addRepresentation("cartoon", param = list(color = defHex))%>%')
###loop to add the colors
color_template<-'addRepresentation("ball+stick", param = list(
    colorScheme = "element",
    colorValue = "HEXCOLOR",
    sele = "POSITION"
  ))%>%'
for(i in 1:nrow(rel_ents)){
  thisEntry<-color_template
  thisEntry<-gsub(pattern = "HEXCOLOR",replacement = rel_ents$color[i],thisEntry)
  thisEntry<-gsub(pattern = "POSITION",replacement = rownames(rel_ents)[i],thisEntry)
  NGL_file<-c(NGL_file,thisEntry)
}


###loop to add the labels
label_template<-'addRepresentation("label",
                    param = list(
                      sele = "POSITION",
                      labelType = "format",
                      labelGrouping = "residue", # or "atom" (eg. sele = "20:A.CB")
                      labelFormat="MYLABEL",
                      color = "white",
                      fontFamiliy = "sans-serif",
                      xOffset = 1,
                      yOffset = 0,
                      zOffset = 0,
                      fixedSize = TRUE,
                      radiusType = 1,
                      radiusSize = 1.5, # Label size
                      showBackground = FALSE
                    )
  )%>%'


###add labels if and only if doLabels is true
if(doLabels){
  for(i in 1:nrow(rel_ents)){
    thisEntry<-label_template
    thisEntry<-gsub(pattern = "MYLABEL",replacement = rel_ents$label[i],thisEntry)
    thisEntry<-gsub(pattern = "POSITION",replacement = rownames(rel_ents)[i],thisEntry)
    #need to get rid of hanging pipe thing
    if(i==nrow(rel_ents)){
      thisEntry<-gsub(pattern = "%>%","",thisEntry,fixed = T)
    }
    NGL_file<-c(NGL_file,thisEntry)
  }
}

NGL_file_name<-paste0(gene_to_plot,"_NGL.R")
writeLines(text = NGL_file,con = NGL_file_name)

#sourcing doesn't work for some reason
#you have to open up the file in R
#run it from a new window
#source(NGL_file_name)




rel_ents


#swapaa ARG #0:79.A

color_one_ces_res<-function(input_pdb,input_nuc_fasta){
  
  
}



plot_effects(effects = cesa$selection$recurrents)
