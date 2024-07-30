# aligning reads to a reference


```sh
minimap2 -a -x map-ont bas54.fa reads.fq > map.sam
```



## if you're bored...

- what is the fraction of forward / reverse mappings?
- are there reads that do not map? What do they contain?
- is there a correlation between the quality score of a site and whether the site differs from the reference?
- what are clipped sequences?

___

[$\leftarrow$ to the previous part](note1.md) \
[$\rightarrow$ to the next part](note3.md)