---
title: Exploring Genetic Structures with gggenes
author: Tom Stanton
date: '2022-03-14'
slug: gggenes-tutorial
categories: []
tags: []
subtitle: ''
summary: ''
authors: []
lastmod: '2022-03-14T15:57:22Z'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: []
---



## Introduction

Why would we want to plot genes? Well there can be many reasons why we would want to do this, not least to spice up any publication or presentation. Plotting genes gives us a visual representation of the genetic organisation of anything of interest in a genome, such as a bio-synthetic gene cluster or operon.

The advantage of doing this in R is that we can filter huge genome data-sets to plot specific things of interest, and we easily grab any public data we want using Entrez. Additionally, any filtering and plotting can be turned into a function which can be automated for lots of repetitive tasks or even interactive use on a website.

## Setting up gggenes, genbankr and rentrez


```r
install.packages("gggenes")
BiocManager::install("genbankr")
install.packages("rentrez")
```

## Import required packages


```r
# Reporting
library(knitr)
# Wrangling
library(tidyverse)
library(stringr)
# Plotting
library(ggplot2)
library(gggenes)
# Genbank fetching and parsing
library(genbankr)
library(rentrez)
```

## Plotting siderophore BGCs on virulence plasmids

#### Search plasmids of interest with rentrez

We can search for a particular virulence plasmid with a high quality sequence using rentrez and a specific search term.


```r
plasmid_search <- entrez_search(
  db="nuccore",
  term='SGH10 AND "Klebsiella pneumoniae"[porgn] AND 
  (refseq[filter] AND plasmid[filter])')

plasmid_search$ids
```

```
## [1] "1302480780"
```

This only returns 1 ID which is exactly the plasmid we're looking for. We can use this ID to fetch the accessions, which returns a new-line separated string, which we can coerce into a vector with `strsplit`.


```r
plasmid_accessions <- entrez_fetch(db="nuccore", id=plasmid_search$ids, rettype="acc")
plasmid_accessions <- strsplit(plasmid_accessions, '\n')[[1]]
plasmid_accessions
```

```
## [1] "NZ_CP025081.1"
```

#### Parse GenBank using plasmid accession

We can use the accession as a `GBAccession` object with `readGenBank`, which will fetch and parse the GenBank file for us, how convenient!


```r
plasmid_accession <- GBAccession(plasmid_accessions[1])
plasmid <- readGenBank(plasmid_accession)
```

#### Filter data-set and format for gene-plotting

As before, we can turn our GenBank coding sequences into a data-frame and filter for siderophore BGCs using a regex. All gene names are present and we don't need to flip the gene coordinates, but we do remove the translation column. Finally, we are also adding a new column called `siderophore` to we can create facets for each siderophore cluster of interest. This is done with `mutate()`, `ifelse()` and a regex using `grepl()`.


```r
siderophores <- as(cds(plasmid), "data.frame") %>% 
  filter(grepl('Iro|Iut|Iuc|aerobactin|salmochelin', product)) %>%
  mutate(middle = (start + end) / 2,
         siderophore=ifelse(
           grepl('Iro|salmochelin', product), 'Salmochelin', 'Aerobactin'))  %>% 
  select(!translation)
```

#### Plot genes

Now we can plot the genes as before where the aesthetics are:

-   `xmin` and `xmax` as the start and end of each gene respectively

-   `y` as the name of each siderophore locus

-   `x` as the midpoint of our gene (for adding the gene labels)

-   `fill` means that each gene will be colored by it's receptive protein product

-   `label` as the name of our gene

We can add gene labels with `geom_text()` and plot each one as a facet with `facet_wrap()`.


```r
siderophores %>%
ggplot(aes(xmin = start, xmax = end, y = siderophore, x=middle,
               fill = product, label = gene)) +
    geom_gene_arrow(arrowhead_height = unit(5, "mm"), 
                    arrow_body_height=unit(5, "mm"), 
                    arrowhead_width = unit(1, "mm"), 
                    show.legend = T) +
    geom_text(angle = 45, hjust = 0.2, nudge_y=0.2, size = 2) +
    facet_wrap(~ siderophore, scales = "free", ncol = 1) +
    labs(y=NULL, x='Length (basepairs)') +
    theme_genes() + theme(legend.position='bottom')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="672" />

## Plotting T3SS genes on *Shigella flexneri* plasmids

Lets try a different genes of interest on plasmids from a different species.

#### Search plasmids of interest with rentrez

Again, we can use rentrez to search for our species of interest, and filter for known plasmid sequences that are in RefSeq.


```r
plasmid_search <- entrez_search(
  db="nuccore",
  term='"Shigella flexneri"[porgn] AND 
  (refseq[filter] AND plasmid[filter])')

plasmid_accessions <- entrez_fetch(db="nuccore", id=plasmid_search$ids, rettype="acc")
plasmid_accessions <- strsplit(plasmid_accessions, '\n')[[1]]
plasmid_accessions
```

```
##  [1] "NZ_CP012139.1"     "NZ_CP012138.1"     "NZ_CP012736.1"    
##  [4] "NZ_CP012734.1"     "NZ_CP012733.1"     "NZ_CP012732.1"    
##  [7] "NZ_CP026773.1"     "NZ_CP024472.1"     "NZ_CP024471.1"    
## [10] "NZ_CP026800.1"     "NZ_LR878367.1"     "NZ_LR878366.1"    
## [13] "NC_024996.1"       "NZ_CP037924.1"     "NZ_WACK01000004.1"
## [16] "NZ_WACK01000003.1" "NZ_WACK01000002.1" "NZ_LR213457.1"    
## [19] "NZ_LR213456.1"     "NZ_LR213454.1"
```

#### Parse GenBank using plasmid accession

Lets pick a plasmid at random and parse it as before.


```r
plasmid_accession <- GBAccession(plasmid_accessions[8])
plasmid <- readGenBank(plasmid_accession)
```

#### Filter data-set and format for gene-plotting

As before, we can turn our GenBank coding sequences into a data-frame and filter for T3SS using a regex with `grepl()`.


```r
t3ss <- as(cds(plasmid), "data.frame") %>% 
  filter(grepl('type 3|type iii|secretion', product, ignore.case = T)) %>%
  mutate(middle = (start + end) / 2,
         gene = coalesce(gene,gene_id)) %>%
  mutate(x = if_else(strand=='-', end, start),
         y = if_else(strand=='-', start, end)) %>%
  mutate(start = x, end = y) %>% select(!c(x, y)) %>% 
  select(!translation)
```

#### Plot genes


```r
t3ss %>%
  ggplot(aes(xmin = start, xmax = end, y = seqnames, x=middle,
               fill = product, label = gene)) +
    geom_gene_arrow(arrowhead_height = unit(5, "mm"), 
                    arrow_body_height=unit(5, "mm"), 
                    arrowhead_width = unit(1, "mm"), 
                    show.legend = T) +
    geom_text(angle = 45, hjust = 0.2, nudge_y=0.2, size = 2) +
    labs(y=NULL, x='Length (basepairs)', title = plasmid_accession,
         subtitle = 'Shigella flexneri plasmid', 
         caption = 'Type III Secretion System Genes') +
    theme_genes() + theme(legend.position='bottom')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" width="672" />
