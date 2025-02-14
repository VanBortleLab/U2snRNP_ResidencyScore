
Banner![Screenshot 2025-01-13 at 6 17 51 PM](https://github.com/user-attachments/assets/657259ee-a6b8-4c6f-9ce8-997a32f6fd55)

# Overview
We applied a uniform framework for scoring U2 snRNP residency at the 3' end of all introns using ENCODE eCLIP data ([ENCODE project](https://www.encodeproject.org/)) performed in HepG2 and K562 cell line, facilitating the understanding of U2 snRNP binding patterns across introns. This repository specifically includes details and links to processed files related to our manuscript:
https://www.biorxiv.org/content/10.1101/2024.07.02.601767v1
# Data processing
We retrieved eCLIP-seq data targeting U2 snRNPs (SF3A3, SF3B4) from ENCODE (see underlying data below).
Next, we applied a Poisson framework for establishing a global U2 occupancy score at the 3' intron-exon junction. Briefly, 500 bp upstream and downstream of every 3’ intron-exon junctions using hg38 as reference genome are binned into 100 10bp binds. Reads mapping to each bin for all of the junctions are extracted from eCLIP-seq using bedtools. Enrichment scores for 5 bins into the introns of the 3’ junction were determined using a Poisson framework. Briefly, for each bin, the expected signal within 500bp upstream (λup), 500bp downstream (λdown), and an average across the entire 1000bp window for all transcript (λwhole) is computed. The λ  value  is  thereafter determined  as  the  maximum  lambda  value  of  the  distance set. The probability density function (PDF) utilized is given by `P(X = x) = (e^(-λ) * λ^x) / x!`. If k represents the observed signal for the bin, the p-value tied to a specific bin enrichment is calculated using P(X > k). Finally, a U2 snRNP residency score is computed as the negative logarithm (base 10).
![Fig1](https://github.com/user-attachments/assets/3be09ee4-2005-461b-8c4c-2be414fdbdfc)
## Sample data
The [sample data](https://github.com/VanBortleLab/U2snRNP_ResidencyScore/blob/main/sample.bed) contains counts extracted from 100 windows with 10bp resolution from U2 snRNP eCLIP-seq.
```r
df <- read.delim("sample.bed")
head(df)
```
| Chr  | Start     | End       | Index                                             | Null | Strand | GSM2423237_eCLIP_SF3B4_K562 | GSM2423238_eCLIP_SF3B4_K562 | GSM2423259_eCLIP_SF3B4_HepG2 | GSM2423260_eCLIP_SF3B4_HepG2 |
|------|-----------|-----------|--------------------------------------------------|------|--------|----------------------------|----------------------------|----------------------------|----------------------------|
| chr1 | 100010864 | 100010873 | chr1_100007156_100011364_ENSG00000283761/clean  | 0    | +      | 0                          | 0                          | 0                          | 0                          |
| chr1 | 100010874 | 100010883 | chr1_100007156_100011364_ENSG00000283761/clean  | 0    | +      | 0                          | 0                          | 0                          | 0                          |
| chr1 | 100010884 | 100010893 | chr1_100007156_100011364_ENSG00000283761/clean  | 0    | +      | 0                          | 0                          | 0                          | 0                          |
| chr1 | 100010894 | 100010903 | chr1_100007156_100011364_ENSG00000283761/clean  | 0    | +      | 0                          | 0                          | 0                          | 0                          |
| chr1 | 100010904 | 100010913 | chr1_100007156_100011364_ENSG00000283761/clean  | 0    | +      | 0                          | 0                          | 0                          | 0                          |

Run `df` through the example [code](https://github.com/VanBortleLab/U2snRNP_ResidencyScore/blob/main/code.R) to get the U2 snRNP residency p.value  
# Data deposit
Our U2 residency p.value over all introns is available for download via [Zenodo](https://zenodo.org/records/13760839) 
# Underlying data
| Protein | Cell Line | Type   | Read 1         | Read 2         |
|---------|-----------|--------|----------------|----------------|
| SF3A3   | HepG2     | PAIRED | SRR5111429_1   | SRR5111429_2   |
| SF3A3   | HepG2     | PAIRED | SRR5111430_1   | SRR5111430_2   |
| SF3B4   | K562      | PAIRED | SRR5111338_1   | SRR5111338_2   |
| SF3B4   | K562      | PAIRED | SRR5111339_1   | SRR5111339_2   |
| SF3B4   | HepG2     | PAIRED | SRR5111367_1   | SRR5111367_2   |
| SF3B4   | HepG2     | PAIRED | SRR5111368_1   | SRR5111368_2   |


