# exploring alignments with IGV

Parsing the sam file can be complicated. There are however convenient tools to visually inspect the aligned reads, such as [IGV](https://igv.org/).

## sorting and indexing mapped reads

As a first step to prepare the data for IGV, we need to
- convert the read aligment from `sam` to `bam` format. The two format are equivalent but the latter is binarized, and thus smaller in size.
- sort the `bam` file, i.e. order the mappings according to which comes first on the reference.
- index the `bam` file, i.e. create an auxiliary file that allows for quick access to the mappings.

```sh
# converting sam to bam file
samtools view -b map.sam > map.bam

# sorting bam file
samtools sort map.bam > map.sorted.bam

# indexing bam file (produces map.sorted.bam.bai)
samtools index map.sorted.bam
```

## visualizing alignments with IGV

To visualize the alignments, open IGV and load the reference genome and the sorted and indexed bam file. You can navigate the genome and zoom in/out to inspect the reads. You can also filter the reads by mapping quality, read length, etc.

## if you're bored...

- why is there a peak in the coverage pattern?
- is there a region with a high number of mismatches? What could be the reason?
- is there a region with many clipped reads?
___

$\leftarrow$ [previous part](note2.md)