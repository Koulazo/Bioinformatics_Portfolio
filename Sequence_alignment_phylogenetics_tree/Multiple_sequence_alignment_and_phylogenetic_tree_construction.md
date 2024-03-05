# Multiple sequence alignment and phylogenetic tree construction
## This project aims to determine the origin of the COVID-19 virus using multiple sequence alignment and phylogenetic tree construction
## To achieve this, we will *detail steps here*
### This project's completion is attributed to my successful completion of the BIS180L course during my attendance at UC Davis. The course materials and resources were provided by Dr. Julin Maloof and can be found on the BIS180L course GitHub page hosted at https://jnmaloof.github.io/BIS180L_web/.


First we acquire the sequences to complete the analysis.

```sh
wget https://bis180ldata.s3.amazonaws.com/downloads/Assignment4/selected_viral_seqs_195v2.fa -P output/
```
## Create Alignments
For the multiple sequence alignment we will use MAFFT.

```sh
#MAFFT installation within Linux kernel
sudo apt install mafft
```

```sh
#Run MAFFT
time mafft --maxiterate 100 --thread 3 --reorder --op 0.5 selected_viral_seqs_195v2.fa  > mafft_maxiter100_195_op.5.fa
```

To view the genome alignments created by MAFFT we will use the in-browser viewer [Wasabi](was.bi).
![Sequence Alignment](images/Sequence_Alignment.PNG)


## Subset sequences in R-Studio



```{r,eval = True}
library(tidyverse)
library(Biostrings)
```



```{r}
inpath <- "output/mafft_maxiter100_195v2_op.5.fa"
outpath <- "output/mafft_maxiter100_195v2_op.5_trimmed_75pct.fa"
alignment <- readDNAMultipleAlignment(inpath)
alignment
```


```{r}
#Trim the beginning and ends with incomplete sequences
alignment <- DNAMultipleAlignment(alignment,start=1000,end=48449)

#Mask sites that contain more than 25% gaps
alignment <- maskGaps(alignment, min.fraction=0.25, min.block.width=1)
maskedratio(alignment) #what proportion got masked? (rows and columns)
```

```{r}
#Change alignments into a stringset to more easily manipulate the names
alignment <- alignment %>% as("DNAStringSet") 
```

```{r}
newnames <- names(alignment) %>% 
  tibble(name=.) %>%
  mutate(name=str_replace_all(name," ","_")) %>% #replace " " with "_" because some programs truncate name at " "
  separate(name, 
           into=c("acc", "isolate", "complete", "name", "country", "host"),
           sep="_?\\|_?") %>%
  mutate(name=str_replace(name, "Middle_East_respiratory_syndrome-related","MERS"),   # abbreviate
         name=str_replace(name, "Severe_acute_respiratory_syndrome-related", "SARS"), # abbreviate
         newname=paste(acc,name,country,host,sep="|")) %>% # select columns for newname
  pull(newname) #return newname
```

```{r}
head(newnames)
```

```{r}
names(alignment) <- newnames
alignment %>% writeXStringSet(outpath)
```


```{bash}
#Create a phylogenetic tree from the subsetted sequences
FastTree -nt -gtr -gamma -out output/mafft_maxiter100_195v2_op.5_trimmed.fasttree.tre output/mafft_maxiter100_195v2_op.5_trimmed_75pct.fa
```

The phylogenetic tree can be viewed at the in-browser tool [iTol](https://itol.embl.de/)

![Coronavirus_genetic_tree](images/coronavirus_tree120.png)

# Results

Seq H has SARS Coronavirus within its taxon and Bat Coronavirus and as sister taxa, putting it well within range for the genetic classification of a Coronavirus. All of the hosts of these nearby viruses are bats in the Rhinolophus family, making it unlikely that Seq_H evolved from any other host.
