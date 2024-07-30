# what you get from sequencing

Sequencing technologies read some type of signal (usually light or current) from which they infer the sequence.
In our case we consider reads obtained with [Nanopore sequencing](https://www.nature.com/articles/s41587-021-01108-x/figures/1). This technology generates fairly long reads. But irrespective of the sequencing method, results will most of the time be obtained in the form of a [fastq file](https://en.wikipedia.org/wiki/FASTQ_format).

If you look at the `reads.fq` file, the first four lines should look like this:
```
@a6b08eb2-38bc-4ee8-bd44-f9ef896b4f76 runid=18083a38f24669f60331f9c85cdc711b77a76b5d sampleid=no_sample read=28793 ch=426 start_time=2023-07-14T18:04:07Z model_version_id=dna_r10.4.1_e8.2_sup@v3.5.1 barcode=barcode02
ATCAAGTGCACAAAAAGATGATGCTGAACTTTGCAATGGGTATCTACATGGACGGGTATAAGGGAGGTGATCCGGTCCTCAAGCTTCGCTGCCCGACACTGTGGGAAAGAGTGGTTGACAGTAAGGACGGCACTGCCCTGCTCGAACTCCACACCGCAGGTCGCGATTTGGCGGTATGTAGACGGATACCCTCTGTGGGCATGTGGGCCTGGGGAATTTGTTGACAGACCTAGCGCGGTAGTGGACTTGCTATTACAGTTTGAATAAAAGGTGAGTTATAAGGTATTCATGCACCTACCAACAACAATACAACGGAGAAAGTGGAATGACGAAAGAAGAGGCCGACGTATTAGATTTTGCACTAGAGGCTACGCAGGAGGCTTTGGTGTTAGCTTGGAAAATGAACCAGCGTGCTGCAGAATGCAACCGCCGCATCAATGAGATTCTGGGCATGTCTATGTATGATTTTACCATGCATGCCGGAGAGATTCCTGAATTAACTGCAGAGATGAATGATCTTGACAATAAAATTAATATCATTCTGAAGACAGCAGAAACAACTCACCAAACCGCCCGTAAATTTATCAACGAGAAGGTGTGAATATGGACTACAGCTACTGCGGGGAAGAAGTTCCTTTTAACTACAGCGATCTTACCACGAGCAATCCATATCGAGAACTCGCCACGGCACGCACATCATCGCCATTGACTTTGATGGAACGTTTAATGCGAGCCCTTTCATCTTCTCCCGCCTGATCAGGAATATAAAGAGGGGTGGACAGTCTACCTTGTGACTTTCCGATTCATTACTCAAGACAATGTAGACATCGAGTATTGGCGCAAGGAGCTCAATGTTGGTGTGTTTTATACTGGTGGTGTTCAGAAGGCAGCCTACCTGGCGCAGTTTGATATCTACCCGGGATCTA
+
*+-//--)(()*6:@MR{LMHIJFF6//0896*+**+,5DFGH>?>>HF@AAACHJFKC544434DF=,,,,,...EJRLGFIFJHHHJKHGJHGIIKLOA??DGKK{PVKKGFDDBDDBGFDCEEFJKNKIFHEGLGE444@??@CCAA83336459854.-,++,-112;?@FBDDC@@@@AFIEEFEEEEFCB?>@>:864**)**-159;HHIKNHGHEEEGIHGFI>?AADDCEDFFFJIEFEECDDFEGA@@CDBBIHIBCEEHKIG:5554520///127777EHDCBDCCC@@@@AE7776<===+)))*))(((((()),--143385214723?>?>>?AABBBEECBABEEOFEDEE@CA<1=44489EDILMMUWbNNGE>>722<=CCA?ABADCCDEEL{PKOIIHIHRKQIJIIFGFC:++++66?@ALGIIGIGIIHEKFKJIFIGIIFM{MPMRM{OJLIHIHKJI{M{IFFF<::8=HHT{MLMGHEFFHIJDKIKEDDDF<:9:>JIIFdMNSILJLP@?????PL{PSNMINOLKJHNILLIGJJMDCCCCNVMPJJHJLI{YJP00003>@@A@?BNMGHLGGE?BBABCDFMOJHHGEDA?CE?OHF<8546666;419:?<=<<>.....@;;?BEFEHJDFKFCHFFE32210213488DDDFDFBD;;;;:>>>?NLJIMIJMNAG]NLPTJMLJJIPKQGMMEDCCDG{LLT@::@?@@==<>?9EGN<<<;:8;(('''&&&')..5>?AEFFGGIHFHFH{NKKMNFHYH00000?@7===<>>>@GIIGJLIIA=====+***,BCC>??AAMJKLPLDJFHGGD{LMIJKIGEM{KMNGFEDDFDFFEEHCCCDF{LGHEDEIB@@AAHFAA@@AE@:::43-,+*+,(+,,+,+'
```

Fastq files are organized in groups of four line, with each line representing:
```
@read_id
sequence
+
quality_string
```

The quality string specifies `confidence` of the base call, that a measure of the error probability of this base. Quality scores are encoded on a log-scale and represented as ASCII characters.
The most commonly used scale is Q = char(33+q) where
$$
q = -10\log_{10} p
$$
and $p$ is the probability of miscall.
Hence a base that the sequencer thinks is 99.9% correct will have a $q=30$ and be encoded as character `?`.
Depending on the sequencing technology, these quality score can be fairly reliable measures of sequencing quality, or only rough indicators.

## if you're bored...

- how long is the reference genome?
- how many reads does the file contain?
- what is the distribution of read length?
- what is the expected coverage if every read maps?
- what is the distribution of read quality?

___

$\rightarrow$ [next part](note2.md)