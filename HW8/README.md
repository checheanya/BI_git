## Project #8. “Immune repertoire annotation: a RepSeq data analysis tutorial“ (Immune Repertoire Sequencing data analysis)

> Anna Chechenina, Kristina Zheltova

> April 2023, Bioinformatics Institute

In this repo you can find:
* task file with goals and described steps
* final table with phenotypes and cell types for all samples
* tutorial in rmarkdown format ([from this repo](https://github.com/mikeraiko/repseq-annotation-tutorial)) and knitted file with data analysis and plots in pdf format

The aim of the project was to perform basic Immune Repertoire Sequencing (RepSeq) data analysis on T-cell receptor (TCR) repertoires. In the data folder you can find txt files for 16 samples of 10000 random reads from two donors from [Qi et al.](http://www.pnas.org/content/111/36/13139.short) study, without sample labels and TCR nucleotide sequences.

We reproduced the code in the rmarkdown file from the tutorial, previously installed needed requirements:

<p><code>install.packages(c("data.table","dplyr","reshape2","ggplot2","NMF","scales","forcats","parallel","stringr"))</code>.</p>

Based on the performed analysis we build a table (.tsv file) identifying cell types for each sample.
As we know, memory CD8 T-cells have lower diversity, whereas CD4 memory cells were the second less diverse subjects. The reason of such a low level of diversity is that these cells are targeting specific antigen. To determine CMV status we compared epitope recognition pattern.

For the patients recognition we used HLA and antigen recognition values, which were the same for the same patient. Also we can distinguish two donors by having or not CMV. Naive T-cells were assigned to patients by comparison of the frequencies of different genes usege and the difference between them and memory cells. Their subsets were defined by analyzing diversity and similarity with the matched memory cells samples.
