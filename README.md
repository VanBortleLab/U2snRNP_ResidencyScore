
![Banner](https://github.com/user-attachments/assets/251045aa-aa6a-4e12-8735-a7193b2b7892)
# Overview
We applied a uniform framework for scoring U2 snRNP residency at the 3' end of all introns using ENCODE eCLIP data ([ENCODE project](https://www.encodeproject.org/)) performed in HepG2 and K562 cell line, facilitating the understanding of U2 snRNP binding patterns across introns. This repository specifically includes details and links to processed files related to our manuscript:
https://www.biorxiv.org/content/10.1101/2024.07.02.601767v1
# Data processing
We retrieved eCLIP-seq data targeting U2 snRNPs (SF3A3, SF3B4) from ENCODE(see underlying data below). 
Next, we applied a Poisson framework for establishing a global U2 occupancy score at the 3' intron-exon junction. Briefly, 500 bp upstream and downstream of every 3’ intron-exon junctions using hg38 as reference genome are binned into 100 10bp binds. Reads mapping to each bin for all of the junctions are extracted from eCLIP-seq using bedtools. Enrichment scores for 5 bins into the introns of the 3’ junction were determined using a Poisson framework. Briefly, for each bin, the expected signal within 500bp upstream (λup), 500bp downstream (λdown), and an average across the entire 1000bp window for all transcript (λwhole) is computed. The λ  value  is  thereafter determined  as  the  maximum  lambda  value  of  the  distance set. The probability density function (PDF) utilized is given by `P(X = x) = (e^(-λ) * λ^x) / x!`. If k represents the observed signal for the bin, the p-value tied to a specific bin enrichment is calculated using P(X > k). Finally, a U2 snRNP residency score is computed as the negative logarithm (base 10).
![Fig1](https://github.com/user-attachments/assets/3be09ee4-2005-461b-8c4c-2be414fdbdfc)
## Sample file

