# ISSAAC-seq by Droplet
This section contains information about ISSAAC-seq using the droplet workflow. Here, the [10x Genomics Single Cell ATAC](https://www.10xgenomics.com/products/single-cell-atac) kit was used, but any droplet systems with a Nextera capture sequence should work, such as [Bio-Rad](https://www.bio-rad.com/en-us/product/surecell-atac-seq-library-prep-kit?ID=PEXSR1MC1ORV) and [HyDrop](https://hydrop.aertslab.org/).

## About the library and files

|          |     Read     |                                         Description                                        |
|:--------:|:------------:|:------------------------------------------------------------------------------------------:|
| __ATAC__ |  Read 1 (R1) |                                          gDNA read                                         |
|          |  Read 2 (R2) |                                          gDNA read                                         |
|          | Index 1 (I1) |                           Sample index (not needed for analysis)                           |
|          | Index 2 (I2) |                                  i5 index (cell barcodes)                                  |
|----------|--------------|--------------------------------------------------------------------------------------------|
|  __RNA__ |  Read 1 (R1) | cDNA sequence, might have some adaptor contamination depending on how long you sequence it |
|          |  Read 2 (R2) |            The first 10 bp are UMIs, the rest are ignored as they are mostly dT            |
|          | Index 1 (I1) |                           Sample index (not needed for analysis)                           |
|          | Index 2 (I2) |                                  i5 index (cell barcodes)                                  |

Check [this page](https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html#Droplet) to see how a step-by-step guide of how the libraries are generated. In this workflow, single nuclei are caputred using the droplet microfluidics after in situ treatment. Since we use the `10x Genomics Single Cell ATAC` kit for both RNA and ATAC libraries, the sequencing configuration is the same for both modalities:

```
> 50 cycles for Read 1 (R1)
> 50 cycles for Read 2 (R2)
8-10 cycles for I1 (i7) <-- This is the sample index
16 cycles for I2 (i5) <-- This is the cell barcode
```

That means you sequence the libraries as if they are 10x scATAC-seq libraries. Only `R1`, `R2` and `I2` are needed for the analysis. The mouse cortex data (two replicates) can be downloaded from here:

- __Rep1 ATAC__

  - [mCortex_rep1_Droplet_ATAC_S1_L001_R1_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_Droplet_ATAC_S1_L001_R2_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_Droplet_ATAC_S1_L001_I2_001.fastq.gz](ENA link to be activated)


- __Rep1 RNA__

  - [mCortex_rep1_Droplet_RNA_S1_L001_R1_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_Droplet_RNA_S1_L001_R2_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_Droplet_RNA_S1_L001_I2_001.fastq.gz](ENA link to be activated)


- __Rep2 ATAC__

  - [mCortex_rep2_Droplet_ATAC_S1_L001_R1_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep2_Droplet_ATAC_S1_L001_R2_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep2_Droplet_ATAC_S1_L001_I2_001.fastq.gz](ENA link to be activated)


- __Rep2 RNA__

  - [mCortex_rep2_Droplet_RNA_S1_L001_R1_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep2_Droplet_RNA_S1_L001_R2_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep2_Droplet_RNA_S1_L001_I2_001.fastq.gz](ENA link to be activated)

## To get the gene expression matrix

1. Collect the cell barcodes and UMI information from `I2` and `R2` read files. Use `collect_cb_umi_issaac_droplet.py` from the `scripts_data` directory and do:

```
# rep1
python collect_cb_umi_issaac_droplet.py \
    mCortex_rep1_Droplet_RNA_S1_L001_R2_001.fastq.gz \
    mCortex_rep1_Droplet_RNA_S1_L001_I2_001.fastq.gz | \
    gzip > mCortex_rep1_Droplet_RNA_CB_UMI.fastq.gz

# rep2
python collect_cb_umi_issaac_droplet.py \
    mCortex_rep2_Droplet_RNA_S1_L001_R2_001.fastq.gz \
    mCortex_rep2_Droplet_RNA_S1_L001_I2_001.fastq.gz | \
    gzip > mCortex_rep2_Droplet_RNA_CB_UMI.fastq.gz
```

2. Simply run [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) (we used `v2.7.9a`):


```
# rep1
STAR --genomeDir /path/to/mouse/star_index \
     --readFilesCommand zcat \
     --readFilesIn mCortex_rep1_Droplet_RNA_S1_L001_R1_001.fastq.gz mCortex_rep1_Droplet_RNA_CB_UMI.fastq.gz \
     --clip3pNbases 116 \
     --soloCBstart 1 --soloCBlen 8 \
     --soloUMIstart 9 --soloUMIlen 10 \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist 737K-cratac-v1_rc.txt \
     --runThreadN 40 \
     --outSAMattributes CB UB --outSAMtype BAM SortedByCoordinate

# rep2
STAR --genomeDir /path/to/mouse/star_index \
     --readFilesCommand zcat \
     --readFilesIn mCortex_rep2_Droplet_RNA_S1_L001_R1_001.fastq.gz mCortex_rep2_Droplet_RNA_CB_UMI.fastq.gz \
     --clip3pNbases 116 \
     --soloCBstart 1 --soloCBlen 8 \
     --soloUMIstart 9 --soloUMIlen 10 \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist 737K-cratac-v1_rc.txt \
     --runThreadN 40 \
     --outSAMattributes CB UB --outSAMtype BAM SortedByCoordinate
```

Note that we used 151 bp for the cDNA reads, which was unnecessarily too long. Since we sequenced the libraries with other people, we don't really have a choice. There will be some adaptor sequence at the 3' end of the the cDNA read, so we only used the first 35 bp for the mapping. That's what the `--clip3pNbases 116` flag is about. For the `737K-cratac-v1_rc.txt` file, it is in the `scripts_data` directory, and this is basically the reverse complementary to the cell barcode sequences in the `737K-cratac-v1.txt.gz` file from [Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/overview/welcome). The original `737K-cratac-v1.txt.gz` file is located in (using v2.0.0 as an example) `cellranger-atac-2.0.0/lib/python/atac/barcodes`. You can get the reverse complements by running:

```
zcat 737K-cratac-v1.txt.gz | rev | tr 'ACGT' 'TGCA' > 737K-cratac-v1_rc.txt
```

If you are still using the old `HiSeq2500` for sequencing, you should use the original `737K-cratac-v1.txt` as the whitelist. If you are using NextSeq and NovaSeq v1.5, you should use `737K-cratac-v1_rc.txt` as the whitelist. You probably need the `737K-cratac-v1_rc.txt` file.

## To get the peak count matrix

This is even simpler. You can just run [Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/overview/welcome) (we used v2.0.0) to get the output. In order to run the program, you need to rename the fastq files according to their convention, which is very strict:

|                   original file                   |                    renamed file                   |
|:-------------------------------------------------:|:-------------------------------------------------:|
| mCortex_rep1_Droplet_ATAC_S1_L001_R1_001.fastq.gz | mCortex_rep1_Droplet_ATAC_S1_L001_R1_001.fastq.gz |
| mCortex_rep1_Droplet_ATAC_S1_L001_R2_001.fastq.gz | mCortex_rep1_Droplet_ATAC_S1_L001_R3_001.fastq.gz |
| mCortex_rep1_Droplet_ATAC_S1_L001_I2_001.fastq.gz | mCortex_rep1_Droplet_ATAC_S1_L001_R2_001.fastq.gz |
| mCortex_rep2_Droplet_ATAC_S1_L001_R1_001.fastq.gz | mCortex_rep2_Droplet_ATAC_S1_L001_R1_001.fastq.gz |
| mCortex_rep2_Droplet_ATAC_S1_L001_R2_001.fastq.gz | mCortex_rep2_Droplet_ATAC_S1_L001_R3_001.fastq.gz |
| mCortex_rep2_Droplet_ATAC_S1_L001_I2_001.fastq.gz | mCortex_rep2_Droplet_ATAC_S1_L001_R2_001.fastq.gz |

Put files associated with each replicate in one directory, and simply run:

```
# rep1
cellranger-atac-2.0.0/cellranger-atac count \
    --id mCortex_ATAC_rep1 \
    --reference /path/to/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
    --fastqs /path/to/rep1/fastq \
    --localcores 40
    
# rep1
cellranger-atac-2.0.0/cellranger-atac count \
    --id mCortex_ATAC_rep2 \
    --reference /path/to/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
    --fastqs /path/to/rep2/fastq \
    --localcores 40
```

Then, aggregate the two replicates using the `aggr` function (check [the manual](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/aggr) for more info):

```
cellranger-atac-2.0.0/cellranger-atac aggr \
    --id mCortex_ATAC_aggr \
    --csv /path/to/sample_info.csv \
    --normalize depth \
    --reference /path/to/refdata-cellranger-arc-mm10-2020-A-2.0.0
```

## Data analysis

First, change the cell barcodes of the RNA part to match the ATAC cell barcodes. The cell barcodes are in the `Solo.out/Gene/filtered/barcodes.tsv`. Modify that file by simply reverse complement the sequence and add `-1` at the end to match the 10x convention:

```
for i in $(cat barcodes.tsv); do
    echo "1-${i}" | rev | tr 'ACGT' 'TGCA'
done > tmp

mv tmp barcodes.tsv
```

Then use the `R` script `mCortex_analysis.R` from the `scripts_data` directory to analyse the data. The `R` script will output the following files that are needed to reproduce the results:

[mCortex_all_metadata.csv](https://github.com/dbrg77/ISSAAC-seq/blob/main/scripts_data/mCortex_all_metadata.csv) # metadata for all cells  
[mCortex_all_expression_matrix.csv.gz](https://github.com/dbrg77/ISSAAC-seq/blob/main/scripts_data/mCortex_all_expression_matrix.csv.gz) # RNA count matrix (UMI)  
[mCortex_all_RNA_UMAP_coordinates.csv](https://github.com/dbrg77/ISSAAC-seq/blob/main/scripts_data/mCortex_all_RNA_UMAP_coordinates.csv) # RNA UMAP coordinates (2D)  
[mCortex_all_ATAC_UMAP_coordinates.csv](https://github.com/dbrg77/ISSAAC-seq/blob/main/scripts_data/mCortex_all_ATAC_UMAP_coordinates.csv) # ATAC UMAP coordinates (2D)  
[mCortex_all_OPC_to_Oligo_peudotime_information.csv](https://github.com/dbrg77/ISSAAC-seq/blob/main/scripts_data/mCortex_all_OPC_to_Oligo_peudotime_information.csv) # Pseudotime information of OPCs and Oligos  
[mCortex_all_OPC_to_Oligo_cor_gene_expression_pseudotime.csv](https://github.com/dbrg77/ISSAAC-seq/blob/main/scripts_data/mCortex_all_OPC_to_Oligo_cor_gene_expression_pseudotime.csv) # Smoothed gene expression along the pseudotime trajectory  
[mCortex_all_OPC_to_Oligo_cor_gene_activity_score_pseudotime.csv](https://github.com/dbrg77/ISSAAC-seq/blob/main/scripts_data/mCortex_all_OPC_to_Oligo_cor_gene_activity_score_pseudotime.csv) # Smoothed gene activity score along the pseudotime trajectory  

Check the [mCortex_cell_types.ipynb](https://nbviewer.org/github/dbrg77/ISSAAC-seq/blob/master/main/mCortex_cell_types.ipynb) and [mCortex_pseudotime.ipynb](https://nbviewer.org/github/dbrg77/ISSAAC-seq/blob/master/main/mCortex_pseudotime.ipynb) notebook files for details about reproducing the figures.
