# aligning reads to a reference

Most sequencers deliver only short DNA fragments and/or make mistakes. Hence the output of the sequencer is quite far from what you'd like to have in the end and substantial post-processing is necessary.

The two main strategies to analyze sequencing data further are (i) mapping to reference genomes or (ii) assembly of all the little reads into longer parts.
- Mapping is the method of choice if a close, reliable reference sequence exists and the emphasis is on single nucleotide changes or small indels.
- Assembly requires longer reads, more coverage, and more computing power but if successful gives a genome where there was none before. This is an important 

In our case we know the reference genome for our sequences: the [Bas54 pahge](https://www.ncbi.nlm.nih.gov/nuccore/2071745857) Aligning the reads to the assembled genome is still however an important sanity check, to spot possible inconsistencies between the reads and the assembly.


## Mapping the reads

In mapping, one piles up sequencing read on their homologous positions of a reference genome.
Mapping programs typically generate an index of the reference genome first and then use this index to rapidly find aligment of the reads on the reference genome.
The positions can be found in a time that logarithmic in the size of the genome and large data sets can be mapped against large genomes very rapidly.
But many-fold repeated regions can slow the mapping process and mapping doesn't allow a large number of differences between the query and the reference.

In our case we will use [minimap2](https://github.com/lh3/minimap2), which is well-adapted for mapping long reads:

```sh
minimap2 -a -x map-ont bas54.fa reads.fq > map.sam
```

With the `-x map-ont` flag minimap adapts its parameters to map ONT reads.


### mapping output

Mapping program typically produce output files in `sam` format, which is compact way to specify where and how each read aligned (mapped) to the reference.
Samfiles typically start with a header, composed of lines starting with the `@` character. The header contains various information, including the names and lengths of the different pieces (chromosomes, contigs, etc) the reference came in, and the command used to map the reads.

In our case for example it is:
```
@HD	VN:1.6	SO:unsorted	GO:query
@SQ	SN:MZ501093.1	LN:136346
@PG	ID:minimap2	PN:minimap2	VN:2.28-r1209	CL:minimap2 -a -x map-ont -t 8 bas54.fa reads.fq
```
Each subsequent line contains the alignment of an individual read:
```
4b9d52af-2178-4e3c-82c5-f0cec1b173ce	16	MZ501093.1	20014	60	5S155M1I83M1I39M1D61M2I24M1D23M1I72M3D23M1D14M1D133M	*	0	0	CTAGAAGCAAGAATGTATTTGCTCATGACCCAGATACTTAGTGCTACTACACGTGGAACGGGGATTATCCTGAAATAATCCTCCCAGACTCCACTGTAGAGGGTGCTGGGGGCGTGTCTGCAAATGCCTGGTCTGTTTTTGGTGAGCTTGCTGCAACTTCTTTCTGGCAGTATAGTTGATTACGGATTTATCGGTGGTCAGCTTGATATGGATCTTGAAGTTGCGGATACCTTTAAGGTTCGACGCAACTTCAAATACCACAATCTCTTTTGAAAATCAAACAGAGGTCTTGAAGGGGTGGCGAGAACTATAACCGTATGCATTACTCAGACCTCTGGCGGCAACAAAAAGTCTATTGGCCTGGTAATGTTAATGGTCCTACGGAAGAGATCCTGATCCTAACATTTACAGCTGGTGAAACAGATATTTTCAAGTTGGAAACCTACGACAACGGGTTAACTTGGTACGCACTGATTATCGCAGGAGCTATGAAGATGCGTATTCGCATAACATTGATAACGTATTTCAGATGATTGAGGGGCACAGGAAGTGCCTAGATAACAACACAGGTATTACCAACACCTCTACACAAAATAAGCTAGTAAACCCTGATGGTCTTATAGCTGCAGGGTGGGCA	%%&&)*)(+'%%%%%&%&&&'('(-..2/,%%%%%&''')((******6221*'((*/(047<@>>@.+*'&*),.3345<<=@???@ACBCKIDDAA@A@76:@A/./2888==HKMHHNMPJIE???<777BKJG>AE@><==<=<==<......1.*.(&&(-***+,69999:MEE999900298::9::;EEKGDAAAA@AB@86666667CKHHFHGGKGMGJNME00:>0//1*))*--212223ACCEOUPQJJQLKOT{LAI38E?KB>><<20+*1:;;:=>FGIHK{{RLLRNLNIIKMPGGHHJHFHICDA???@CJIJJGI=....--,10**&&'**++,>?AIDFCDEEDFIFD65=<843,,++0+*+/.-.--,,./197779AHDDFHEIIIHD@*))+/125;=LMKKMLJEEEADJLGFNFG76667HFFEEIG@8731123++('''(+++DSGH{JHI>54444('''(0111111222,,,**++;53334JGGFFGT{JKJKJLKMPHH000728,++++5,)))))((((4446;;==@DJJJNHLOLLJIM{NJIHHFEEDD;13553388BBEGIHKQLKDCCCHHHJZQKIQLODCAB43:=;<?@@?D	NM:i:26	ms:i:1113	AS:i:1110	nn:i:0	tp:A:P	cm:i:74	s1:i:490	s2:i:0	de:f:0.0362	rl:i:0
```
This is a tab-separated format. Let's break it down in some relevant fields:

| col | value                                  | description                                                                   |
| --- | -------------------------------------- | ----------------------------------------------------------------------------- |
| 1   | `4b9d52af-2178-4e3c-82c5-f0cec1b173ce` | read id                                                                       |
| 2   | `16` ( = `b000010000000`)              | [bit-packed flag](https://broadinstitute.github.io/picard/explain-flags.html) |
| 3   | `MZ501093.1`                           | reference id                                                                  |
| 4   | `20014`                                | read length                                                                   |
| 5   | `60`                                   | read mapping quality                                                          |
| 6   | `5S 155M 1I 83M 1I 39M 1D 61M ...`     | [CIGAR string](https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/)  |
| ... | ...                                    | ...                                                                           |
| 10  | `CTAGAAGCAAGAATGTATTTGCTCATGA...`      | read sequence                                                                 |
| 11  | `%%&&)*)(+'%%%%%&%&&&'('(-..`          | read quality score                                                            |

The first entry is the read ID, the second is a [complicated flag with various information on the mapping](https://broadinstitute.github.io/picard/explain-flags.html), the third is the reference sequence id. Then comes the mapping position (base 20014). The fifth is the so called [CIGAR string](https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/) that specifies the aligment to the reference sequence. In this case, after 6 soft-clipped bases, it specifies a 155 base match (155M), a one base insertion (1I), and another 83 base match (83M) etc. Fields 10 and 11 contain the sequence and the quality score.
The full format specification can be found [here](https://samtools.github.io/hts-specs/SAMv1.pdf)
Note that the samfile contains all the information from the fastq file and information on where the read mapped.

We won't go into the details of how the mapping and alignment is done algorithmically.
Most programs initally generate an index of the reference genome (a structure that allows rapid look-up of the positions of kmers), find seeds in reads using the index, and do a local pairwise aligment.

If you want to explore the sam file more, check out [this notebook](scripts/exploring_sam_file.ipynb).

## if you're bored...

- what is the fraction of forward / reverse mappings?
- how many secondary/supplementary alignments are there? Where do most supplementary alignments start/end?
- are there reads that do not map? What do they contain?
- is there a correlation between the quality score of a site and whether the site differs from the reference?
- what are clipped sequences?
- what is the distribution of in/del sizes?

___

$\leftarrow$ [previous part](note1.md) \
$\rightarrow$ [next part](note3.md)