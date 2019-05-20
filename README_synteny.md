# Comparative macrosynteny analysis using SynChro

***This tutorial takes EMBL/Genbank genome files or genome FASTAs and GFF annotations, extracts protein FASTA records from FASTAs and GFFs, and runs all genomes through SynChro.*** *Here we use this protocol to get genome files from* Gallus gallus *and* Anolis carolinensis, *and genome FASTA and GFF annotation from* Gopherus agassizii. *We pull protein FASTA records from* Gopherus agassizii *and isolate proteins from the top 26 scaffolds with the most genes. We input these files into the CHROnicle directory and run SynChro to create macrosynteny figures.*



### Step 1: Clone in GitHub repository at https://github.com/mmoral31/GopAga2.0-Genome.git

*The scripts will run straight out of the box (note—these were only tested on unix and linux platforms)*

```bash
> git clone https://github.com/mmoral31/GopAga2.0-Genome.git
```



### Step 2: Gather EMBL or Genbank (.dat) genome files for species of interest from ftp.ensembl.org. If no genome files exist for your species of interest, gather genome FASTAs and GFF annotations.

some species have curated EMBL or Genbank genome files on Ensembl that can be easily converted for use in SynChro. however, if the species of interest does not have these files, the genome FASTA and GFF annotation can be curated in steps 3 and 4 for use in SynChro.



###Step 3: For species without EMBL or Genbank genome files, use gff2fasta.pl to pull protein FASTA records from genome FASTAs and GFF annotations.

> Requires BioPerl

to use genome FASTAs and GFF annotations in SynChro, protein FASTA records must be extracted from both files and supplemented with data on gene position and scaffold. this modified gff2fasta.pl script will perform this task.

gff2fasta.pl takes three parameters in this order:

- genome FASTA for species of interest (e.g., Cmyd_genome.fa)
- GFF annotation for species of interest (e.g., Cmyd_annotation.gff)
- shortened species name for use in the prefix of output files (e.g., Cmyd)

```bash
> python gff2fasta.pl Cmyd_genome.fa Cmyd_annotation.gff Cmyd
```

only modifications made to the gff2fasta.pl script at https://github.com/ISUgenomics/common_scripts/blob/master/gff2fasta.pl is the addition of gene start, gene end, strand, and scaffold data for each gene in the ID  of its respective protein.

output file for further curation will be named {shortened species name}.pep.fasta (e.g., Cmyd.pep.fasta)



### Step 4: To reduce memory requirements in SynChro, use pull_id_match.py to isolate proteins from a reduced number of scaffolds.

> Requires Biopython

preliminary analyses using SynChro with protein FASTAs showed that genomes with genes spread across too many scaffolds will quickly reach memory limitations and will cause the program to stall. to combat this, only include proteins from genes in scaffolds that contain the greatest number of genes. this will allow for sufficient coverage without memory errors.

first, parse the protein FASTA of the species of interest to determine scaffolds with the most genes.

```bash
grep '>' Cmyd_proteins.fa | cut -d " " -f5 | sort | uniq -c | sort
```

example output:

```bash
 364 ScCC6lQ_38465_HRSCAF_64658
 376 ScCC6lQ_161303_HRSCAF_192837
 395 ScCC6lQ_161354_HRSCAF_192888
 401 ScCC6lQ_3500_HRSCAF_17357
 406 ScCC6lQ_97_HRSCAF_1492
 415 ScCC6lQ_28138_HRSCAF_52803
 420 ScCC6lQ_3789_HRSCAF_18069
 429 ScCC6lQ_161370_HRSCAF_192904
 446 ScCC6lQ_14958_HRSCAF_36502
 450 ScCC6lQ_161315_HRSCAF_192849
 454 ScCC6lQ_3764_HRSCAF_18001
 455 ScCC6lQ_40749_HRSCAF_67263
 503 ScCC6lQ_161394_HRSCAF_192928
 508 ScCC6lQ_24292_HRSCAF_48220
 528 ScCC6lQ_140_HRSCAF_2235
 580 ScCC6lQ_1923_HRSCAF_13295
 715 ScCC6lQ_16796_HRSCAF_38896
 718 ScCC6lQ_161421_HRSCAF_192955
 861 ScCC6lQ_161298_HRSCAF_192832
```

then, create a text file with the names of desired scaffolds; one scaffold name per line.

```bash
ScCC6lQ_14958_HRSCAF_36502
ScCC6lQ_161315_HRSCAF_192849
ScCC6lQ_3764_HRSCAF_18001
ScCC6lQ_40749_HRSCAF_67263
ScCC6lQ_161394_HRSCAF_192928
ScCC6lQ_24292_HRSCAF_48220
ScCC6lQ_140_HRSCAF_2235
ScCC6lQ_1923_HRSCAF_13295
ScCC6lQ_16796_HRSCAF_38896
ScCC6lQ_161421_HRSCAF_192955
ScCC6lQ_161298_HRSCAF_192832
```

finally, run pull_id_match.py isolate proteins from desired scaffolds.

pull_id_match.py takes three parameters in this order:

- name of protein FASTA file to pull from (e.g., Cmyd_proteins.fa)
- name of text file with list of desired scaffolds (e.g., Cmyd_scaffolds.txt)
- name of output FASTA file (e.g., Cmyd_protreduct.fa)

```bash
> python pull_id_match.py Cmyd_proteins.fa Cmyd_scaffolds.txt Cmyd_protreduct.fasta
```



### Step 5: Convert genome files and protein FASTAs to SynChro format using ConvertGbk.py, ConvertEMBL.py, or ConvertFasta.py for Genbank, EMBL, and protein FASTA files respectively.

> Requires CHROnicle 

genome files and protein FASTAs must now be converted to SynChro format for downstream analyses. **from this point forward, all scripts are run from inside the CHROnicle folder** that can be downloaded from [http://www.lcqb.upmc.fr/CHROnicle/SynChro.html](http://www.lcqb.upmc.fr/CHROnicle/SynChro.html). 

create a directory for the current macrosynteny analysis (e.g. Reptiles). inside this directory, create a directory "00rawGenom". in this directory, create a directory for each species in this analysis (e.g., CheloniaMydas, GopherusAgassizii, etc.). in each of these directories, place the respective genome files or protein FASTAs. 

from the top level of the CHROnicle directory, move into the "Programs/0Convert2InputF" directory. from here, run one of three scripts (*ConvertGbk.py*, *ConvertEMBL.py*, or *ConvertFasta.py*) for each species.

all scripts require at least three parameters in this order:

- name of the directory for the current macrosynteny analysis (e.g. Reptiles)
- name of the directory for the target species (e.g. CheloniaMydas)
- shortened species name for output files (e.g. CMYD)

if the target species has Genbank genome files, run ConvertGbk.py with no additional parameters.

```bash
> ./ConvertGbk.py Reptiles CheloniaMydas CMYD
```

if the target species has EMBL genome files, run ConvertEMBL.py with no additional parameters.

```bash
> ./ConvertEMBL.py Reptiles CheloniaMydas CMYD
```

if the target species has a protein FASTA, run ConvertFasta.py with additional parameters specifying the columns of each protein FASTA record ID that contain information on gene name, scaffold, start, end, and strand. the last parameter specifies whether the protein FASTA is sorted or not. **the example parameterization is sufficient for output from *gff2fasta.pl***.

*note: ConvertFasta.py will only work if protein FASTA has the suffix ".fasta" and will otherwise throw an InputF error*

```bash
> ./ConvertFasta.py Reptiles CheloniaMydas CMYD 1 5 2 3 4 0
```



###Step 6: Run SynChro.py to perform macrosynteny analyses

now, SynChro.py can be run to perform macrosynteny analyses.

from the top level of the CHROnicle directory, move into the "Programs/1SynChro" directory. from here, run SynChro.py once.

SynChro.py takes three parameters in this order:

- name of the directory for the current macrosynteny analysis (e.g. Reptiles)
- parameter set to 0 if all species in the directory are to be analyzed or 1 if only a select few are to be analyzed
- delta value which specifies how many orthologous genes need to be syntenic in both genomes to be considered a synteny block (for more information, visit: [http://www.lcqb.upmc.fr/CHROnicle/SynChro.html](http://www.lcqb.upmc.fr/CHROnicle/SynChro.html))

```bash
> ./SynChro.py Reptiles 0 2
```

