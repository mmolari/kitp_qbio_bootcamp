# Bioinformatic bootcamp notes - read alignment to reference

This small repository contains notes for the bioinformatics bootcamp of the Qbio course, in occasion of the 2024 KITP program [Horizontal Gene Transfer and Mobile Elements in Microbial Ecology and Evolution](https://online.kitp.ucsb.edu/online/hgt24/directory.html)

They cover the topic of aligning reads to a reference. We consider the example of a dataset obtained by Nanopore sequencing a phage isolate. The tutorial is divided in three sections:
- [a look into the fastq file](note1.md)
- [aligning to a reference](note2.md)
- [inspecting alignments with IGV](note3.md)

## setup

Get the required software by:
- installing `samtools` and `minimap2` using conda or mamba:
    ```sh
    conda create -n aln_env -c bioconda -c conda-forge minimap2 samtools
    ```
    remeber to activate the environment before using the tools!
- download [IGV](https://igv.org/doc/desktop/#DownloadPage/).

you should have in your folder:
- the fast file `reads.fq` containing the reads
- the reference sequence `bas54.fa`, corresponding to [phage Bas54](https://www.ncbi.nlm.nih.gov/nuccore/2071745857) from the Basel collection.