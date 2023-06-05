# Introduction
Below the assigned cases and their segregation pattern 

* 1677 Autosomal Dominant 
* 1689 Autosomal Dominant
* 1832 Autosomal Dominant
* 1743 Autosomal Recessive
* 1796 Autosomal Dominant

Each number represent a triplets - father, mother and child. Each exome sequence is saved in a server in .fq.gz file format.


# Alignment and indexing 
For each locus to which a sequence can be mapped, determine the optimal base by base alignment of the query sequence to the reference sequence. Output of bowtie2 is a .sam alignment file, with SAMtools we obtain a sorted .bam file
```
bowtie2 -U case1677_father.fq.gz -p 8 -x /home/BCG2023_genomics_exam/uni --rg-id 
'F' --rg "SM:father" | samtools view -Sb | samtools sort -o case1677_father.bam

bowtie2 -U case1677_mother.fq.gz -p 8 -x /home/BCG2023_genomics_exam/uni --rg-id 
'M' --rg "SM:mother" | samtools view -Sb | samtools sort -o case1677_mother.bam

bowtie2 -U case1677_child.fq.gz -p 8 -x /home/BCG2023_genomics_exam/uni --rg-id 
'C' --rg "SM:child" | samtools view -Sb | samtools sort -o case1677_child.bam
```
Parameters: 
* -U → Comma-separated list of files containing unpaired reads to be aligned, so the FASTAQ file
* -p → when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes. Set to 8
* -x → The basename of the index for the reference genome, in this case files uni in the exam directory
* —rg-id → Set the read group ID to `<text>`. This causes the SAM `@RG` header line to be printed, with `<text>` as the value associated with the ID: tag. It also causes the`RG:Z:`
extra field to be attached to each SAM output record, with value set to `text`. In this case: M mother, F father and C child
* —rg →Add `<text>` as a field on the @RG header line. Note: in order for the @RG line to appear, --rg-id must also be specified.
* samtools view -Sb → specifies that the input is in SAM format and the output will be be BAM format
samtools sort -o → sets the name of the output file

    >This informations are taken from the offical documentation 

This sorted BAM alignment file can now be indexed: 
```
samtools index case1677_father.bam
samtools index case1677_mother.bam
samtools index case1677_child.bam
```

# Quality check 
Before starting any type of analysis is best practice to perform some quality check to esure that data are of high quality. In this project i use 
fastQC, Qualimap and MultiQC. In the final report you can find some screenshot of the reports.
 ## FastQC
 This is one of the most popular tools for quality check and it aims to provide a simple way to do some quality control checks on raw sequence data. In order to perform a complete analysis I will merge those file with MultiQC.
 ```
 fastqc case1677_child.fq.gz
fastqc case1677_mother.fq.gz
fastqc case1677_father.fq.gz
```

## Qualimap
>Qualimap examines sequencing alignment data in SAM/BAM files according to the features of the mapped reads and provides an overall view of the data that helps to the detect biases in the sequencing and/or mapping of the data and eases decision making for further analysis - see official documentation.

I apply this to perform genomewide coverage analysis 
```
qualimap bamqc -bam case1677_father.bam -gff /home/BCG2023_genomics_exam/exons16Padded_sorted.bed -outdir case1677_father
qualimap bamqc -bam case1677_mother.bam -gff /home/BCG2023_genomics_exam/exons16Padded_sorted.bed -outdir case1677_mother
qualimap bamqc -bam case1677_child.bam -gff /home/BCG2023_genomics_exam/exons16Padded_sorted.bed -outdir case1677_child
```
Note that I'm using a target region file provided by the professor. 

## MultiQC 
Let's put together the information of the three samples. MultiQC help to summarize experiments containing multiple samples and multiple analysis steps, such as FastQC and Qualimap. 

In the final report you can find the interpretation of the resultant quality check. 

# Coverage with bedtools 
The next step is computing coverage track with bedtools.

>`bedtools genomecov` computes histograms (default), per-base reports (`-d`) and BEDGRAPH (`-bg`) summaries of feature coverage (e.g., aligned sequences) for a given genome.

```
bedtools genomecov -ibam case1677_father.bam -bg -trackline -trackopts 'name="father"' -max 100 >  case1677_fatherCov.bg

bedtools genomecov -ibam case1677_mother.bam -bg  -trackline -trackopts 'name="mother"' -max 100 >  case1677_motherCov.bg

bedtools genomecov -ibam case1677_child.bam -bg  -trackline -trackopts 'name="child"' -max 100 >  case1677_childCov.bg
```
Parameters: 
* -ibam →BAM file as input for coverage. 
* -bg → Report depth in BedGraph format, a concise representation since consecutive positions with the same coverage are reported as a single output line describing the start and end coordinate of the interval having the coverage level. 
* -trackline → Adds a UCSC/Genome-Browser track line definition in the first line of the output. This parameter is needed by BedGraph
* -trackopts →Writes additional track line definition parameters in the first line, such as the name for BedGraph
* max → Combine all positions with a depth >= max into a single bin in the histogram. (like the x50 for qualimp)

> These description are taken from the official documentation

# Variant calling - Freebayes

> Freebayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment. - from Offical documentation

```
freebayes -f /home/BCG2023_genomics_exam/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 case1677_child.bam case1677_mother.bam case1677_father.bam > case1677.vcf
```

Parameters:
* -f → fasta reference, so the universe.fasta file in the exam directory. It contains the sequence of chr16
*- -t → targets FILE. Limit analysis to targets listed in the BED-format FILE.
- -m n → -min-mapping-quality, exclude alignments from analysis if they have a mapping quality less than n.
- -C → Require at least 5 supporting observations to consider a variant, standard value is 5
- -Q n → -mismatch-base-quality-threshold, count mismatches toward --read-mismatch-limit if the base quality of the mismatch is >= n.
- --min-coverage → Require at least this coverage to process a site.

> These description are taken from the official documentation

The output of this a VCF file.

# Selecting segregation pattern
To do this I need to chech the samples order in the VCF file
```
less -L case1677.vcf #father, mother and child
```
Then I check the possible genotype for mother, father and child to filter the "impossible" combinations 
```
grep -v "#" case1677.vcf | cut -f 10 | cut -d ":" -f 1 | sort | uniq -c
grep -v "#" case1677.vcf | cut -f 11 | cut -d ":" -f 1 | sort | uniq -c
grep -v "#" case1677.vcf | cut -f 12 | cut -d ":" -f 1 | sort | uniq -c
```

## Autosomal dominant case
This case was reported to be affected by an autosomal dominant disease: both parents are healthy and so they are homozygous while the child must be heterozygous → de novo mutation so I need to select case where only the child is mutated (and so different from the parents):

```
grep "#" case1677.vcf > candilist1677.vcf
grep "0/0.*0/0.*0/1" case1677.vcf >> candilist1677.vcf
grep "0/0.*0/0.*0/2" case1677.vcf >> candilist1677.vcf
grep "1/1.*1/1.*1/2" case1677.vcf >> candilist1677.vcf
grep "0/0.*0/0.*0/3" case1677.vcf >> candilist1677.vcf
grep "2/2.*2/2.*1/2" case1677.vcf >> candilist1677.vcf
```

This result in a new VCF file with the usual header and just the correct case for my analysis. 

## Autosomal recessive case
Assuming healthy parents means that we select cases with heterozygous mother and father,
while the affected child must be homozygous with the correct allele inherited from parents.
```
grep "#" case1743.vcf > candilist1743.vcf
grep "0/1.*0/1.*1/1" case1743.vcf >> candilist1743.vcf
grep "1/2.*1/2.*2/2" case1743.vcf >> candilist1743.vcf
grep "1/2.*1/2.*1/1" case1743.vcf >> candilist1743.vcf
grep "0/2.*1/2.*2/2" case1743.vcf >> candilist1743.vcf
grep "1/2.*0/2.*2/2" case1743.vcf >> candilist1743.vcf
grep "0/1.*1/2.*1/1" case1743.vcf >> candilist1743.vcf
grep "1/2.*0/1.*1/1" case1743.vcf >> candilist1743.vcf
grep "0/2.*0/2.*2/2" case1743.vcf >> candilist1743.vcf
grep "1/3.*0/1.*1/1" case1743.vcf >> candilist1743.vcf
grep "1/3.*1/2.*1/1" case1743.vcf >> candilist1743.vcf
grep "2/3.*1/2.*2/2" case1743.vcf >> candilist1743.vcf
grep "2/3.*0/2.*2/2" case1743.vcf >> candilist1743.vcf
```

## Intersect with bedtools
Use bedtools intersect to filter for “target regions”

```
bedtools intersect -a candilist1677.vcf -b /home/BCG2023_genomics_exam /exons16Padded_sorted.bed -u > 1677candilistTG.vcf
```

# Conclusions
I obtained the final VCF files for each samples assigned. In the final analysis report you can find additional informations about quality chech and the diagnosis made for each samples

