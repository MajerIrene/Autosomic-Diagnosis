# Molecular diagnosis of rare Mendelian autosomic disease 
This is a simple workflor for diagnosis of autosomic disease from exon sequencing data, a genomic technique for sequencing all of the protein-coding regions of genes in a genome. It's a cheaper alternative to other technologies such as WGS, although is less precise in capturing complex variants. 
In this project the aim is to identify the correct diagnosis, having at disposal the exome sequencing data of an individual and of his parents.

## General Information
The analysis was performed over data provided during Genomics course at University of Milan. 
The data were studied taking into account some assumptions:
- Disease-causing variant are located on chromosome 16
- A set of the possible diseases was given a priori
- The parents were supposed to be heathy
- The reference genome was the hg19.37
- High-quality sequencing data is available – probabilistic modelling of SNPs is not necessary

Here you can find the [Workflow](https://github.com/MajerIrene/Autosomic-Diagnosis/blob/main/Workflow.md) used fo the analysis.

Here you can find the final [Report](https://github.com/MajerIrene/Autosomic-Diagnosis/blob/main/Genomics_report.md).

## Technologies Used
- bowtie2
- SAMtools 
- FastQC
- Qualimap
- Bedtools
- Freebayes
- VEP
- UCSC genome browser
