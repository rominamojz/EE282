# Homework 3 - Pipeline
In Homework 3, we analyze the Drosophila melanogaster genome by processing its FASTA and GTF files from FlyBase. Using tools like faSize and bioawk, we calculate genome metrics and summarize annotation features, including gene counts per chromosome arm. File integrity is verified with checksums, and all code and results are documented directly in the current file.
## Summarizing Genome Assembly 
1.  We download our file using using the following command
```
wget  ftp://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz
```
2. Now we check the integrity of the file using md5sum. The process of checking the integrity of the file goes at follows
```
md5sum dmel-all-chromosome-r6.48.fasta.gz > md5check.txt
md5sum -c md5check.txt
```
3. Using faSize to analyze the data from the file that we downloaded.
```
faSize dmel-all-chromosome-r6.48.fasta.gz 
```
| Total number of nucleotides | Total number of N's | Total number of sequences |
|-----------------------------|---------------------|---------------------------|
| 143726002                  | 1152978            | 1870                      |
## Summarizing Annotation File
1. We download the annotation file from the flybase using this command:
```
wget ftp://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/gtf/dmel-all-r6.60.gtf.gz
```
2. We checked the integrity of the file using similar command as we used for genome assembly section:
```
md5sum dmel-all-r6.60.gtf.gz > md5heck.txt
md5sum -c md5check.txt
```
3. We analyzed the data on the file using bioawk. 
a. The following command is to get the total number of features of each type, sorted from the most common to the least common:
```
bioawk -c gff '{ print $feature }' dmel-all-r6.60.gtf.gz | sort | uniq -c | sort -nr
```
Here's the result of the command:
| Type         | Count    |
|--------------|----------|
| exon         | 190038   |
| CDS          | 163253   |
| 5UTR         | 46806    |
| 3UTR         | 33741    |
| start_codon  | 30888    |
| stop_codon   | 30828    |
| mRNA         | 30802    |
| gene         | 17872    |
| ncRNA        | 3059     |
| miRNA        | 485      |
| pseudogene   | 365      |
| tRNA         | 312      |
| snoRNA       | 270      |
| pre_miRNA    | 262      |
| rRNA         | 115      |
| snRNA        | 32       |

b. The following command is to get the total number of genes per chromosome arm (X, Y, 2L, 2R, 3L, 3R, 4)
``` 
bioawk -c gff '$3 == "gene" { print $1 }' dmel-all-r6.60.gtf.gz | sort | uniq -c | grep -E " (X|Y|2L|2R|3L|3R|4)$"
```

| Gene         | Count    |
|--------------|----------|
|  2L          |  3508    |
| 2R           | 3649     |
| 3L           | 3481     |
| 3R           | 4226     |
| 4            | 114      |
| X            | 2704     |
| Y            | 113      |
