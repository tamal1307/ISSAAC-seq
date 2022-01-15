# ISSAAC-seq by FACS
This section contains information about ISSAAC-seq using a plate-based workflow. FACS is commonly used to sort single nuclei into wells, but handpick should also work.

## About the libraries and files

|          |     Read     |                                         Description                                        |
|:--------:|:------------:|:------------------------------------------------------------------------------------------:|
| __ATAC__ |  Read 1 (R1) |                                          gDNA read                                         |
|          |  Read 2 (R2) |                                          gDNA read                                         |
|          | Index 1 (I1) |                                          i7 index                                          |
|          | Index 2 (I2) |                                          i5 index                                          |
|----------|--------------|--------------------------------------------------------------------------------------------|
|  __RNA__ |  Read 1 (R1) |            The first 10 bp are UMIs, the rest are ignored as they are mostly dT            |
|          |  Read 2 (R2) | cDNA sequence, might have some adaptor contamination depending on how long you sequence it |
|          | Index 1 (I1) |                                          i7 index                                          |
|          | Index 2 (I2) |                                          i5 index                                          |

Check [this page](https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html#FACS) to see how a step-by-step guide of how the libraries are generated. In this workflow, single nuclei are sorted into invidividual wells containing index primers. The library preparation is done individually. In this case, the well barcode is the cell barcode. The combination of `Index 1 (i7)` and `Index 2 (i5)` defines a cell. There are two common scenarios: 

1. You equence your libraries from a core facility. In this case, you probably need to provide the index information `Index 1 (i7) + Index 2 (i5)`, and the core will demultiplex for you. In this case, you have two fastq files for each cell per modality: `R1.fastq.gz` and `R2.fastq.gz`.

2. You sequence by yourself, and run `bcl2fastq` from `Illumina` to generate fastq files. You can certainly put the `Index 1 (i7) + Index 2 (i5)` information in the `SampleSheet.csv`, and each cell will get demultiplexed by `bcl2fastq`. A simpler way is simply put `NNNNNNNN + NNNNNNNN` in the `SampleSheet.csv`. Then run `bcl2fastq` like this:

```
bcl2fastq -r 4 -p 4 -w 4 --create-fastq-for-index-reads --no-lane-splitting -o /path/to/output/dir
```

Due to a known bug in `bcl2fastq`, the program look for the literal `N` in the index. Therefore, you will get the following four files:

```
Undetermined_S0_L001_R1_001.fastq.gz
Undetermined_S0_L001_R2_001.fastq.gz
Undetermined_S0_L001_I1_001.fastq.gz
Undetermined_S0_L001_I2_001.fastq.gz
```

You can start from there, and rename them to something meaningful. The mouse cortex data using the FACS workflow from the paper can be downloaded from here:

- __ATAC__

  - [mCortex_rep1_FACS_ATAC_S1_R1_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_FACS_ATAC_S1_R2_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_FACS_ATAC_S1_I1_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_FACS_ATAC_S1_I2_001.fastq.gz](ENA link to be activated)


- __RNA__

  - [mCortex_rep1_FACS_RNA_S2_R1_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_FACS_RNA_S2_R2_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_FACS_RNA_S2_I1_001.fastq.gz](ENA link to be activated)
  - [mCortex_rep1_FACS_RNA_S2_I2_001.fastq.gz](ENA link to be activated)

## To get the gene expression matrix

1. Collect the cell barcodes and UMI information from `I1`, `I2` and `R1` read files. Use `collect_cb_umi_issaac_plate.py` from the `scripts_data` directory and do:

```
python collect_cb_umi_issaac_plate.py \
    mCortex_rep1_FACS_RNA_S2_I1_001.fastq.gz \
    mCortex_rep1_FACS_RNA_S2_I2_001.fastq.gz \
    mCortex_rep1_FACS_RNA_S2_R1_001.fastq.gz | \
    gzip > mCortex_rep1_FACS_RNA_CB_UMI.fastsq.gz
```

2. Simply run [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) (we used `v2.7.9a`):


```
STAR --genomeDir /path/to/mouse/star_index \
     --readFilesCommand zcat \
     --readFilesIn mCortex_rep1_FACS_RNA_S2_R2_001.fastq.gz mCortex_rep1_FACS_RNA_CB_UMI.fastsq.gz \
     --clip3pNbases 116 \
     --soloCBstart 1 --soloCBlen 8 \
     --soloUMIstart 9 --soloUMIlen 10 \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist mcortex_rep1_facs_rna_whitelist.txt \
     --runThreadN 40 \
     --outSAMattributes CB UB --outSAMtype BAM SortedByCoordinate
```

Note that we used 151 bp for the cDNA reads, which was unnecessarily too long. Since we sequenced the libraries with other people, we don't really have a choice. There will be some adaptor sequence at the 3' end of the the cDNA read, so we only used the first 35 bp for the mapping. That's what the `--clip3pNbases 116` flag is about. For the `mcortex_rep1_facs_rna_whitelist.txt`, this is basically the combination of your `Index 1 (i7)` and `Index 2 (i5)` index sequences. The one used for the mouse cortex FACS experiment can be found under the `scripts_data` directory.

## To get the peak count matrix

A simple and elegant way is to use the recently developed [chromap](https://github.com/haowenz/chromap) aligning method. First, we need to get the cell barcodes into one fastq file. You can use `concat_i7i5.py` from the `scripts_data` directory and do:

```
python concat_i7i5.py \
    mCortex_rep1_FACS_ATAC_S1_I1_001.fastq.gz \
    mCortex_rep1_FACS_ATAC_S1_I2_001.fastq.gz | \
    gzip > mCortex_rep1_FACS_ATAC_CB.fastq.gz
```

Then, simply run `chromap` by:

```
chromap --preset atac \
        -x index \
        -r ref.fa \
        -1 mCortex_rep1_FACS_ATAC_S1_R1_001.fastq.gz \
        -2 mCortex_rep1_FACS_ATAC_S1_R2_001.fastq.gz \
        -o aln.bed \
        -b mCortex_rep1_FACS_ATAC_CB.fastq.gz
        --barcode-whitelist mcortex_rep1_facs_atac_whitelist.txt
```

Again, the `mcortex_rep1_facs_atac_whitelist.txt` is basically the combination of your `Index 1 (i7)` and `Index 2 (i5)` index sequences. The one used for the mouse cortex FACS experiment can be found under the `scripts_data` directory.

However, during the preparation of the manuscript, we were using the old way. That is, we demultiplexed each single cells into individual `fastq` files, and used the [scATAC_snakemake](https://github.com/dbrg77/scATAC_snakemake) pipeline for the data processing. If you like this way, you can use [deML](https://github.com/grenaud/deml) for the demultiplexing.