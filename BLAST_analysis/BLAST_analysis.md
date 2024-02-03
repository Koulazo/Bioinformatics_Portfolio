# Blast Analysis
## This project aims to simulate a genomic investigation to identify the probable cause of a respiratory illness. To achieve this, we will conduct BLAST analyses on patient samples, comparing them against a database of known viral genomes.
### This project's completion is attributed to my successful completion of the BIS180L course during my attendance at UC Davis. The course materials and resources were provided by Dr. Julin Maloof and can be found on the BIS180L course GitHub page hosted at https://jnmaloof.github.io/BIS180L_web/.


To begin the BLAST I downloaded the Refseq files from the host server and unzipped their contents

```sh
wget https://bis180ldata.s3.amazonaws.com/downloads/ncbi_virus_110119_2.txt.gz

wget https://bis180ldata.s3.amazonaws.com/downloads/patient_viral.txt.gz

gunzip -c ncbi_virus_110119_2.txt.gz > ncbi_virus_110119_2.txt

gunzip -c Data/patient_viral.txt.gz > Data/patient_viral.txt
```

To search the entries of the ncbi_virus_110119_2.txt using BLAST, a database is required.

```sh
sudo apt install ncbi-blast+

makeblastdb -in ncbi_virus_110119_2.txt -dbtype nucl
```

We BLAST this database against the patient_viral data to filter for significant results(<1e-3). Then create a .tsv of the output to organize the results.

```sh
time blastn -db ../input/ncbi_virus_110119_2.txt -query ../input/patient_viral.txt -evalue 1e-3 > blastout.default

time blastn  -db ../input/ncbi_virus_110119_2.txt -query ../input/patient_viral.txt -evalue 1e-3 -outfmt '7 std stitle' > blastout.default.tsv
```

## First 10 BLAST Results

The information contained within this BLAST indicates that Seq_A, sampled from human patients, is the complete genome of the Hepatitis B virus. However this is only the BLAST for one of the sample sequences and does not point to any respiratory infection. If we limit the result outputs it becomes easier to look at all the samples and be more confident in our results.


| query acc.ver | subject acc.ver | % identity | alignment length | mismatches | gap opens | q. start | q. end | s. start | s. end | evalue | bit score | subject title | definition | genome completion | origin | country | species sampled |
|---------------|------------------|------------|-------------------|------------|-----------|----------|--------|----------|--------|--------|------------|----------------|-------------|---------------------|--------|---------|-------------------|
| Seq_A         | HQ700546        | 99.067     | 3215             | 30         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5775      | HQ700546 |Hepatitis B virus isolate BC_NC_AM263| complete genome|Hepatitis B virus||                        |
| Seq_A         | MF674513        | 98.974     | 3215             | 33         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5755      | MF674513 |Hepatitis B virus isolate VN05T0253| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens      |
| Seq_A         | KP341010        | 98.818     | 3215             | 38         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5727      | KP341010 |Hepatitis B virus isolate bba9| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens           |
| Seq_A         | MF674488        | 98.694     | 3216             | 40         | 2         | 1        | 3215   | 1        | 3215   | 0.0    | 5705      | MF674488 |Hepatitis B virus isolate VN05A0034| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens      |
| Seq_A         | MF674384        | 98.600     | 3215             | 45         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5688      | MF674384 |Hepatitis B virus isolate HBV040008| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens      |
| Seq_A         | MF674453        | 98.569     | 3215             | 46         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5683      | MF674453 |Hepatitis B virus isolate HBV050074| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens      |
| Seq_A         | MF674424        | 98.569     | 3215             | 46         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5683      | MF674424 |Hepatitis B virus isolate HBV040080| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens      |
| Seq_A         | MF674470        | 98.538     | 3215             | 47         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5677      | MF674470 |Hepatitis B virus isolate HBVM1014| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens       |
| Seq_A         | MF674452        | 98.476     | 3215             | 49         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5666      | MF674452 |Hepatitis B virus isolate HBV050071| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens      |
| Seq_A         | MF674510        | 98.476     | 3215             | 49         | 0         | 1        | 3215   | 1        | 3215   | 0.0    | 5666      | MF674510 |Hepatitis B virus isolate VN05T0230| complete genome|Hepatitis B virus|Viet Nam|Homo sapiens      |



# test limit sequences
```sh
time blastn  -max_target_seqs 5 -db ../input/ncbi_virus_110119_2.txt -query ../input/patient_viral.txt -evalue 1e-3 -outfmt '7 std stitle' > limited_query_blastout.default.tsv

time blastn  -max_target_seqs 5 -perc_identity 98 -db ../input/ncbi_virus_110119_2.txt -query ../input/patient_viral.txt -evalue 1e-3 -outfmt '7 std stitle' > limited_query_perc_ID_blastout.default.tsv
```


## Best Hits
### These filtered BLAST queries allows us to find the best hits on each sequence from the patient.

Sequence G corresponds with Staphylococcus aureus and sequence F corresponds with Streptococcus pneumoniae. Both of these sequences have low evalues, but the bit score for the Streptococcus pneumoniae is significantly higher, indicating that the match is more significant and is not due to random chance. While there is still a chance of the respiratory infection being attributed to Staphylococcus aureus, the statistical evidence favors Streptococcus pneumoniae and the cause of infection.

| query acc.ver | subject acc.ver | % identity | alignment length | mismatches | gap opens | q. start | q. end | s. start | s. end | evalue   | bit score | subject title | definition                                                | genome completion  | origin                        | country sequenced | species sampled          |
|---------------|-----------------|------------|------------------|------------|-----------|----------|--------|----------|--------|----------|-----------|---------------|-----------------------------------------------------------|--------------------|-------------------------------|-------------------|--------------------------|
| Seq_A         | HQ700546        | 99.067     | 3215             | 30         | 0         | 1        | 3215   | 1        | 3215   | 0.0      | 5775      | HQ700546      | Hepatitis B virus isolate BC_NC_AM263                     |  complete genome   | Hepatitis B virus             |                   |                          |
| Seq_B         | MG926379        | 90.651     | 1690             | 135        | 21        | 1        | 1683   | 3        | 1676   | 0.0      | 2226      | MG926379      | UNVERIFIED: Hepatitis delta virus isolate 2000074_243110  |  complete genome   | Hepatitis delta virus         | Israel            | Homo sapiens             |
| Seq_C         | BK000394        | 100.000    | 166              | 0          | 0         | 233668   | 233833 | 192345   | 192180 | 8.57e-79 | 307       | BK000394      | TPA_inf: Human herpesvirus 5 strain AD169 substrain varUK |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | BK000394        | 100.000    | 33               | 0          | 0         | 190815   | 190847 | 190375   | 190407 | 7.37e-05 | 62.1      | BK000394      | TPA_inf: Human herpesvirus 5 strain AD169 substrain varUK |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | BK000394        | 100.000    | 33               | 0          | 0         | 76       | 108    | 190407   | 190375 | 7.37e-05 | 62.1      | BK000394      | TPA_inf: Human herpesvirus 5 strain AD169 substrain varUK |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | BK000394        | 100.000    | 33               | 0          | 0         | 76       | 108    | 229705   | 229737 | 7.37e-05 | 62.1      | BK000394      | TPA_inf: Human herpesvirus 5 strain AD169 substrain varUK |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | BK000394        | 100.000    | 33               | 0          | 0         | 190815   | 190847 | 229737   | 229705 | 7.37e-05 | 62.1      | BK000394      | TPA_inf: Human herpesvirus 5 strain AD169 substrain varUK |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 166              | 0          | 0         | 233668   | 233833 | 193850   | 193685 | 8.57e-79 | 307       | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 33               | 0          | 0         | 191386   | 191418 | 603      | 571    | 7.37e-05 | 62.1      | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 33               | 0          | 0         | 191386   | 191418 | 191309   | 191341 | 7.37e-05 | 62.1      | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 33               | 0          | 0         | 190815   | 190847 | 191880   | 191912 | 7.37e-05 | 62.1      | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 33               | 0          | 0         | 76       | 108    | 191912   | 191880 | 7.37e-05 | 62.1      | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 33               | 0          | 0         | 76       | 108    | 231210   | 231242 | 7.37e-05 | 62.1      | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 33               | 0          | 0         | 190815   | 190847 | 231242   | 231210 | 7.37e-05 | 62.1      | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 32               | 0          | 0         | 77       | 108    | 1        | 32     | 2.65e-04 | 60.2      | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | FJ527563        | 100.000    | 32               | 0          | 0         | 190815   | 190846 | 32       | 1      | 2.65e-04 | 60.2      | FJ527563      | Human herpesvirus 5 strain AD169                          |  complete genome   | Human betaherpesvirus 5       | USA               | Homo sapiens             |
| Seq_C         | AC146999        | 100.000    | 137              | 0          | 0         | 195338   | 195474 | 231013   | 231149 | 1.13e-62 | 254       | AC146999      | Human Herpesvirus 5 AD169-BAC isolate                     |  complete sequence | Human betaherpesvirus 5       |                   |                          |
| Seq_C         | AC146999        | 100.000    | 33               | 0          | 0         | 76       | 108    | 30700    | 30732  | 7.37e-05 | 62.1      | AC146999      | Human Herpesvirus 5 AD169-BAC isolate                     |  complete sequence | Human betaherpesvirus 5       |                   |                          |
| Seq_C         | AC146999        | 100.000    | 33               | 0          | 0         | 190815   | 190847 | 30732    | 30700  | 7.37e-05 | 62.1      | AC146999      | Human Herpesvirus 5 AD169-BAC isolate                     |  complete sequence | Human betaherpesvirus 5       |                   |                          |
| Seq_C         | AC146999        | 100.000    | 33               | 0          | 0         | 191386   | 191418 | 31349    | 31317  | 7.37e-05 | 62.1      | AC146999      | Human Herpesvirus 5 AD169-BAC isolate                     |  complete sequence | Human betaherpesvirus 5       |                   |                          |
| Seq_C         | AC146999        | 100.000    | 33               | 0          | 0         | 191386   | 191418 | 31759    | 31727  | 7.37e-05 | 62.1      | AC146999      | Human Herpesvirus 5 AD169-BAC isolate                     |  complete sequence | Human betaherpesvirus 5       |                   |                          |
| Seq_C         | AC146999        | 100.000    | 33               | 0          | 0         | 191386   | 191418 | 221698   | 221730 | 7.37e-05 | 62.1      | AC146999      | Human Herpesvirus 5 AD169-BAC isolate                     |  complete sequence | Human betaherpesvirus 5       |                   |                          |
| Seq_C         | AC146999        | 100.000    | 33               | 0          | 0         | 190815   | 190847 | 222338   | 222370 | 7.37e-05 | 62.1      | AC146999      | Human Herpesvirus 5 AD169-BAC isolate                     |  complete sequence | Human betaherpesvirus 5       |                   |                          |
| Seq_C         | AC146999        | 100.000    | 33               | 0          | 0         | 76       | 108    | 222370   | 222338 | 7.37e-05 | 62.1      | AC146999      | Human Herpesvirus 5 AD169-BAC isolate                     |  complete sequence | Human betaherpesvirus 5       |                   |                          |
| Seq_C         | KP745687        | 100.000    | 35               | 0          | 0         | 191110   | 191144 | 194501   | 194535 | 5.70e-06 | 65.8      | KP745687      | Human herpesvirus 5 strain BE/36/2011                     |  complete genome   | Human betaherpesvirus 5       | Belgium           | Homo sapiens             |
| Seq_C         | KP745658        | 100.000    | 35               | 0          | 0         | 191110   | 191144 | 194506   | 194540 | 5.70e-06 | 65.8      | KP745658      | Human herpesvirus 5 strain BE/14/2012                     |  complete genome   | Human betaherpesvirus 5       | Belgium           | Homo sapiens             |
| Seq_D         | KP745638        | 100.000    | 32               | 0          | 0         | 236045   | 236076 | 194951   | 194920 | 2.68e-04 | 60.2      | KP745638      | Human herpesvirus 5 strain BE/15/2010                     |  complete genome   | Human betaherpesvirus 5       | Belgium           | Homo sapiens             |
| Seq_D         | KP745638        | 100.000    | 32               | 0          | 0         | 195494   | 195525 | 194920   | 194951 | 2.68e-04 | 60.2      | KP745638      | Human herpesvirus 5 strain BE/15/2010                     |  complete genome   | Human betaherpesvirus 5       | Belgium           | Homo sapiens             |
| Seq_D         | KP745638        | 100.000    | 32               | 0          | 0         | 361      | 392    | 194951   | 194920 | 2.68e-04 | 60.2      | KP745638      | Human herpesvirus 5 strain BE/15/2010                     |  complete genome   | Human betaherpesvirus 5       | Belgium           | Homo sapiens             |
| Seq_F         | MK044829        | 100.000    | 54               | 0          | 0         | 12053    | 12106  | 16505    | 16558  | 2.51e-17 | 100       | MK044829      | Streptococcus phage 109751                                |  complete genome   | Streptococcus phage 109751    | Thailand          | Streptococcus            |
| Seq_F         | NC_024361       | 100.000    | 46               | 0          | 0         | 4221     | 4266   | 32402    | 32357  | 7.03e-13 | 86.1      | NC_024361     | Streptococcus phage DCC1738                               |  complete sequence | Streptococcus phage DCC1738   |                   | Streptococcus pneumoniae |
| Seq_G         | EU136189        | 100.000    | 31               | 0          | 0         | 9901     | 9931   | 9914     | 9884   | 7.34e-05 | 58.4      | EU136189      | Staphylococcus phage SAP-2                                |  complete genome   | Staphylococcus phage SAP-2    | South Korea       | Staphylococcus aureus    |
| Seq_G         | KU992911        | 100.000    | 31               | 0          | 0         | 9901     | 9931   | 9692     | 9662   | 7.34e-05 | 58.4      | KU992911      | Staphylococcus phage SLPW                                 |  complete genome   | Staphylococcus phage SLPW     | China             | Staphylococcus aureus    |
| Seq_G         | MN098325        | 100.000    | 31               | 0          | 0         | 9901     | 9931   | 9882     | 9852   | 7.34e-05 | 58.4      | MN098325      | Staphylococcus phage Portland                             |  complete genome   | Staphylococcus phage Portland | USA               | Staphylococcus aureus    |
| Seq_G         | NC_031008       | 100.000    | 31               | 0          | 0         | 9901     | 9931   | 9692     | 9662   | 7.34e-05 | 58.4      | NC_031008     | Staphylococcus phage SLPW                                 |  complete genome   | Staphylococcus phage SLPW     | China             | Staphylococcus aureus    |



