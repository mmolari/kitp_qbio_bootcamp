# exploring alignments with IGV

```sh
# converting sam to bam file
samtools view -b map.sam > map.bam

# sorting bam file
samtools sort map.bam > map.sorted.bam

# indexing bam file (produces map.sorted.bam.bai)
samtools index map.sorted.bam
```
___

[$\leftarrow$ to the previous part](note2.md)