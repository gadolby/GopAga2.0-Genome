# Automating gene family evolution analyses

**_This tutorial takes a list of homologous genes of interest, BLASTs to obtain these sequences from NCBI, combines the outputs, and aligns them with ClustalW._** _Here we use this protocol to get 10 Toll-Like Receptor genes from tetrapods with available genomes on NCBI (N=22). We divide the sequences into TLR gene subfamilies (known a priori), and use the auto-alignments to generate gene family trees_



### Step 1: Clone in GitHub repository at https://github.com/mmoral31/GopAga2.0-Genome.git

_The scripts will run straight out of the box (note—these were only tested on unix and linux platforms)_

```bash
> git clone https://github.com/mmoral31/GopAga2.0-Genome.git
```



### Step 2: Gather protein FASTA to be searched against and create a BLAST database

> Requires BLAST+

if you want to query against whole proteomes, search for the target species on NCBI and pull protein FASTAs (e.g., *Chelonia mydas*)

```bash
> makeblastdb -in GCF_000344595.1_CheMyd_1.0_protein.faa -dbtype prot -parse_seqids -out CmydPro
```

**_repeat for each species you want to have these genes for_**



### Step 3: Set up query sequences nd input file

**Need to:** 

**1)** combine query reference sequences in a single FASTA, and 

**2)** set up input text file with gene names

determine the query sequences to be used and place into a single FASTA file labeled as {gene family}_ref.fasta where {gene family} is the name given to the group of genes used as a query (e.g., TLR_ref.fasta)

```
>NP_003256.1 toll-like receptor 3 precursor [Homo sapiens]
MRQTLPCIYFWGGLLPFGMLCASSTTKCTVSHEVADCSHLKLTQVPDDLPTNITVLNLTHNQLRRLPAAN
FTRYSQLTSLDVGFNTISKLEPELCQKLPMLKVLNLQHNELSQLSDKTFAFCTNLTELHLMSNSIQKIKN
NPFVKQKNLITLDLSHNGLSSTKLGTQVQLENLQELLLSNNKIQALKSEELDIFANSSLKKLELSSNQIK
EFSPGCFHAIGRLFGLFLNNVQLGPSLTEKLCLELANTSIRNLSLSNSQLSTTSNTTFLGLKWTNLTMLD
LSYNNLNVVGNDSFAWLPQLEYFFLEYNNIQHLFSHSLHGLFNVRYLNLKRSFTKQSISLASLPKIDDFS
FQWLKCLEHLNMEDNDIPGIKSNMFTGLINLKYLSLSNSFTSLRTLTNETFVSLAHSPLHILNLTKNKIS
KIESDAFSWLGHLEVLDLGLNEIGQELTGQEWRGLENIFEIYLSYNKYLQLTRNSFALVPSLQRLMLRRV
ALKNVDSSPSPFQPLRNLTILDLSNNNIANINDDMLEGLEKLEILDLQHNNLARLWKHANPGGPIYFLKG
LSHLHILNLESNGFDEIPVEVFKDLFELKIIDLGLNNLNTLPASVFNNQVSLKSLNLQKNLITSVEKKVF
GPAFRNLTELDMRFNPFDCTCESIAWFVNWINETHTNIPELSSHYLCNTPPHYHGFPVRLFDTSSCKDSA
PFELFFMINTSILLIFIFIVLLIHFEGWRISFYWNVSVHRVLGFKEIDRQTEQFEYAAYIIHAYKDKDWV
WEHFSSMEKEDQSLKFCLEERDFEAGVFELEAIVNSIKRSRKIIFVITHHLLKDPLCKRFKVHHAVQQAI
EQNLDSIILVFLEEIPDYKLNHALCLRRGMFKSHCILNWPVQKERIGAFRHKLQVALGSKNSVH
>NP_003257.1 toll-like receptor 4 isoform C [Homo sapiens]
MELNFYKIPDNLPFSTKNLDLSFNPLRHLGSYSFFSFPELQVLDLSRCEIQTIEDGAYQSLSHLSTLILT
GNPIQSLALGAFSGLSSLQKLVAVETNLASLENFPIGHLKTLKELNVAHNLIQSFKLPEYFSNLTNLEHL
DLSSNKIQSIYCTDLRVLHQMPLLNLSLDLSLNPMNFIQPGAFKEIRLHKLTLRNNFDSLNVMKTCIQGL
AGLEVHRLVLGEFRNEGNLEKFDKSALEGLCNLTIEEFRLAYLDYYLDDIIDLFNCLTNVSSFSLVSVTI
ERVKDFSYNFGWQHLELVNCKFGQFPTLKLKSLKRLTFTSNKGGNAFSEVDLPSLEFLDLSRNGLSFKGC
CSQSDFGTTSLKYLDLSFNGVITMSSNFLGLEQLEHLDFQHSNLKQMSEFSVFLSLRNLIYLDISHTHTR
VAFNGIFNGLSSLEVLKMAGNSFQENFLPDIFTELRNLTFLDLSQCQLEQLSPTAFNSLSSLQVLNMSHN
NFFSLDTFPYKCLNSLQVLDYSLNHIMTSKKQELQHFPSSLAFLNLTQNDFACTCEHQSFLQWIKDQRQL
LVEVERMECATPSDKQGMPVLSLNITCQMNKTIIGVSVLSVLVVSVVAVLVYKFYFHLMLLAGCIKYGRG
ENIYDAFVIYSSQDEDWVRNELVKNLEEGVPPFQLCLHYRDFIPGVAIAANIIHEGFHKSRKVIVVVSQH
FIQSRWCIFEYEIAQTWQFLSSRAGIIFIVLQKVEKTLLRQQVELYRLLSRNTYLEWEDSVLGRHIFWRR
LRKALLDGKSWNPEGTVGTGCNWQEATSI
>NP_057646.1 toll-like receptor 7 precursor [Homo sapiens]
MVFPMWTLKRQILILFNIILISKLLGARWFPKTLPCDVTLDVPKNHVIVDCTDKHLTEIPGGIPTNTTNL
TLTINHIPDISPASFHRLDHLVEIDFRCNCVPIPLGSKNNMCIKRLQIKPRSFSGLTYLKSLYLDGNQLL
EIPQGLPPSLQLLSLEANNIFSIRKENLTELANIEILYLGQNCYYRNPCYVSYSIEKDAFLNLTKLKVLS
LKDNNVTAVPTVLPSTLTELYLYNNMIAKIQEDDFNNLNQLQILDLSGNCPRCYNAPFPCAPCKNNSPLQ
IPVNAFDALTELKVLRLHSNSLQHVPPRWFKNINKLQELDLSQNFLAKEIGDAKFLHFLPSLIQLDLSFN
FELQVYRASMNLSQAFSSLKSLKILRIRGYVFKELKSFNLSPLHNLQNLEVLDLGTNFIKIANLSMFKQF
KRLKVIDLSVNKISPSGDSSEVGFCSNARTSVESYEPQVLEQLHYFRYDKYARSCRFKNKEASFMSVNES
CYKYGQTLDLSKNSIFFVKSSDFQHLSFLKCLNLSGNLISQTLNGSEFQPLAELRYLDFSNNRLDLLHST
AFEELHKLEVLDISSNSHYFQSEGITHMLNFTKNLKVLQKLMMNDNDISSSTSRTMESESLRTLEFRGNH
LDVLWREGDNRYLQLFKNLLKLEELDISKNSLSFLPSGVFDGMPPNLKNLSLAKNGLKSFSWKKLQCLKN
LETLDLSHNQLTTVPERLSNCSRSLKNLILKNNQIRSLTKYFLQDAFQLRYLDLSSNKIQMIQKTSFPEN
VLNNLKMLLLHHNRFLCTCDAVWFVWWVNHTEVTIPYLATDVTCVGPGAHKGQSVISLDLYTCELDLTNL
ILFSLSISVSLFLMVMMTASHLYFWDVWYIYHFCKAKIKGYQRLISPDCCYDAFIVYDTKDPAVTEWVLA
ELVAKLEDPREKHFNLCLEERDWLPGQPVLENLSQSIQLSKKTVFVMTDKYAKTENFKIAFYLSHQRLMD
EKVDVIILIFLEKPFQKSKFLQLRKRLCGSSVLEWPTNPQAHPYFWQCLKNALATDNHVAYSQVFKETV
```

create a text file with each line containing the name of one gene family member in the order of the reference FASTA (e.g., TLR.txt)

```
TLR3
TLR4
TLR7
```



### Step 4: Use _blastp.py_ to extract genes from the target species

> Requires BLAST+ and Biopython

blastp.py takes three parameters in this order:

- name of the created BLAST database without suffix (e.g., CmydPro)
- name of the gene family for query (e.g., tlr). this is important for distinguishing the reference FASTA (e.g., TLR_ref.fasta)
- name of the gene family text file (e.g., TLR.txt)

```bash
> python blastp.py CmydPro tlr TLR.txt
```

this will output one file per gene family member to record standard output BLAST results and one aggregate file with all pulled sequences from the target species.

- {gene family member}_{db name} for all of the individual records

- {gene family}_{db name}.fa for the aggregate FASTA file

  

### Step 5: Use _pull_id_match.py_ to separate output into subfamilies for analysis

> Requires Biopython

if the gene family of interest has separate subfamilies (as in TLRs), use pull_seqids.py to separate and extract a subfamily of interest

pull_id_match.py takes three parameters in this order:

- name of the FASTA file to extract from (e.g., TLR_CmydPro.fa)
- name of the text file with listed keywords to pull with (e.g., TLR3SF.txt)
- name of the output file where pulled sequences will aggregate (e.g., TLR3SF.fa)

the text file will have one keyword on each line. the FASTA output will include any records that matched at least one of the keywords in the list.

```
TLR3
TLR4
```

```bash
> python pull_id_match.py TLR_CmydPro.fa TLR3SF.txt TLR3SF.fa
```



### Step 6: Use _run_alignments.py_ to align respective subfamilies with ClustalW

> Requires ClustalW2

subfamilies need to be aligned in order to proceed with phylogenetic reconstruction. run_alignments.py can be used to align respective subfamilies.

run_alignments.py takes one parameter:

- name of the text file with listed names of FASTA files (without suffix) to be aligned (e.g., TLR_align.txt)

```
TLR3SF
```

```bash
> python run_alignments.py TLR_align.txt
```

three alignment files are output in three formats: FASTA, NEXUS, PHYLLIP



### Step 7: Determine substitution model

Here we used Prottest to determine the best substitution matrix for tree reconstruction

> Requires ProtTest

maximum likelihood phylogenetic reconstruction requires a specified substitution matrix. Prottest can be used to determine the ideal substitution matrix for that data. requires phyllip alignment input.

```bash
> prottest3 -i TLR3SF_alignment.phy -o TLR3SF_models
```



### Step 8: Generate gene family trees 

Here, we used the CIPRES Science Gateway https://www.phylo.org/portal2 (running RaxxML and MrBayes) to generate ML and Bayesian subfamily trees.

