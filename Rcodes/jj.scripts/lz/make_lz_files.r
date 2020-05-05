# Make the data file (--metal option) & the batch mode specification file (--hitspec option)
# for LocusZoom from the input file containing the SNPs (RSIDs) and the pvalues.

# Remove all variables from the workspace
rm(list=ls())

# input files - QC'ed MAGWAS files
file_list <- list.files(path="/home3/fsakhtar/1000g_cell_line/magwas_qced", full.names=TRUE)

# 1000g significant hits at -log(p) >= 6
sig_hits <- read.table('/home3/fsakhtar/1000g_cell_line/1000gs.top.genes.csv',
    header=FALSE, stringsAsFactors=FALSE)
colnames(sig_hits) <- c('drug','ChrNum', 'RSID','basepairPosition',
    'nAA', 'nAa', 'naa','nAA', 'nAa', 'naa','pvalue','logp','gene')

# output dir for LocusZoom output
lz_dir <- c("/home3/fsakhtar/1000g_cell_line/lz/input_files/")

for (f in file_list) {
    print(paste0("Reading file: ", f))
    file_name <- gsub(".*/", "", f)
    df.input <- read.table(f, header=TRUE, stringsAsFactors=FALSE)

    ## Create METAL file for each drug
    # Get required columns for --metal file
    metal_cols <- c('RSID','pvalue_Pillai')
    df.metal <- df.input[, metal_cols]
    colnames(df.metal) <- c('MarkerName', 'P-value')

    # checking that there are no dups
    dups <- subset(df.metal, duplicated(MarkerName))
    if (nrow(dups) > 0) {
        print(paste0('ERROR: dups in ', file_name))
        print(dups)
    }

    metal_file <- c(paste0(lz_dir, file_name, '_metal'))
    write.table(df.metal, metal_file, row.names=FALSE, quote=FALSE, sep='\t')


    ## Create hitspec file for each drug
    drug_sig_hits <- subset(sig_hits, drug==file_name)
    df.hitspec <- data.frame(snp=factor(), chr=numeric(), start=numeric(),
        end=numeric(), flank=factor(), run=factor(), m2zargs=factor())

    if (nrow(drug_sig_hits) > 0) {
        # 
        # sig_snps <- drug_sig_hits$RSID
        # ref_cmd <- rep('add-refsnps=',length(sig_snps))
        # ref_str <- paste0(ref_cmd,sig_snps)
        # locuszoom ... --denote-markers-file <your file>

        # Create an entry in the hitspec file for each chromosome that has a
        # significant SNP
        num_chr <- unique(drug_sig_hits$ChrNum)
        for (chr in num_chr) {
            # All significant entries for this chromosome number
            snps_per_chr <- subset(drug_sig_hits, ChrNum==chr)

            # Plot +/-500kb flanking the region from the smallest and the largest
            # significant SNP on this chromosome
            flank <- 500000
            min_bp <- max((min(snps_per_chr$basepairPosition) - flank), 1)
            max_bp <- max(snps_per_chr$basepairPosition) + flank

            hit_entry <- data.frame(snp='NA', chr=chr, start=min_bp, end=max_bp,
                flank='NA', run='yes', m2zargs='showAnnot=T')
            df.hitspec <- rbind(df.hitspec, hit_entry)
        }
        hitspec_file <- c(paste0(lz_dir, file_name, '_hitspec'))
        write.table(df.hitspec, hitspec_file, row.names=FALSE, quote=FALSE, sep='\t')
    }
}