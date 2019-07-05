# Automating gene family evolution analyses

**_This tutorial takes a list of homologous genes of interest, BLASTs to obtain these sequences from NCBI, combines the outputs, and aligns them with ClustalW._** _Here we use this protocol to get 10 Toll-Like Receptor genes from tetrapods with available genomes on NCBI (N=22). We divide the sequences into TLR gene subfamilies (known a priori), and automatically aligin the subfamily sequences into alignments  to use for phylogenetic reconstruction_



### Step 1: Clone GitHub repo

_The scripts will run straight out of the box (note—these were only tested on unix and linux platforms)_

```bash
> git clone https://github.com/mmoral31/GopAga2.0-Genome.git
```



### Step 2: Make a BLAST database for your species of interest (i.e. target species)

> Requires BLAST+

if you want to query against whole proteomes, search for the target species on NCBI and pull protein FASTAs (e.g., *Chelonia mydas*)

```bash
> makeblastdb -in GCF_000344595.1_CheMyd_1.0_protein.faa -dbtype prot -parse_seqids -out CmydPro
```

**_[repeat for each species of interest]_**



### Step 3: Set up query sequences and input file

**Need to:** 

**1) combine query sequences into a single FASTA**

​	_in this case we’re using TLR protein sequences from human. You want to pick a query taxon that is well curated (i.e. will have the most complete suite of genes with the best orthology assignments)_

**2) set up input text file with gene names**

determine the query sequences to be used and place into a single FASTA file labeled as {gene family}_query.fasta where {gene family} is the name given to the group of genes used as a query (e.g., TLR_query.fasta)

e.g., **_query genes:_**

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

e.g., **_gene family text file:_**

```
TLR3
TLR4
TLR7
```



### Step 4: Use _blastp.py_ to extract query genes from the target species

> Requires BLAST+ and Biopython

_blastp.py_ takes three parameters in this order:

- name of the created BLAST database _without suffix_ (**e.g., CmydPro**)
- name of the gene family for query (e.g., tlr). This is important for distinguishing the reference FASTA (**e.g., TLR_query.fasta**)
- name of the gene family text file (**e.g., TLR.txt**)

```bash
> python blastp.py CmydPro tlr TLR.txt
```

this will output one file per gene family member to record standard output BLAST results and one aggregate file with all pulled sequences from the target species.

- {gene family member}_{db name} for all of the individual records
- {gene family}_{db name}.fa for the aggregate FASTA file (**e.g., TLR_CmydPro.fa**)



### Step 5: Use _pull_id_match.py_ to separate output into subfamilies for analysis

> Requires Biopython

if the gene family of interest has subfamilies you want to align and analyze separately (as in TLRs), use pull_seqids.py to separate and extract a subfamily of interest _[if you don’t want to subset your data you can skip this step]_

_pull_id_match.py_ takes three parameters in this order:

- name of the FASTA file to extract from (**e.g., TLR_CmydPro.fa**)
- name of the text file with listed keywords to pull with (**e.g., TLR3SF.txt**)
- name of the output file where pulled sequences will aggregate (**e.g., TLR3SF.fa**)

The text file will have one keyword on each line. The FASTA output will include any records that matched at least one of the keywords in the list.

e.g., **text file:**

```
TLR3
TLR4
```

```bash
> python pull_id_match.py TLR_CmydPro.fa TLR3SF.txt TLR3SF.fa
```



### Step 6: Use _run_alignments.py_ to align respective subfamilies with ClustalW

> Requires ClustalW2

to align subfamilies for further analysis (e.g., tree reconstructions) use **_run_alignments.py_** 

_run_alignments.py_ takes one parameter:

- name of the text file with listed names of FASTA files (without suffix) to be aligned (e.g., TLR_align.txt) _Note—these files will be aligned separately. This is useful if you want to make many separate alignments_

e.g., **text file:**

```
TLR3SF
```

```bash
> python run_alignments.py TLR_align.txt
```

three alignment files are output in three formats: FASTA, NEXUS, PHYLIP

_Note—ClustalW needs to be in your Bash profile path to be recognized_



### Step 7: Determine substitution model

Here we used Prottest to determine the best substitution matrix for tree reconstruction

> Requires ProtTest

Prottest can be used to determine the best substitution matrix for your data. Requires PHYLIP alignment as input.

```bash
> prottest3 -i TLR3SF_alignment.phy -o TLR3SF_models
```



### Step 8: Generate gene family trees 

Here, we used the CIPRES Science Gateway https://www.phylo.org/portal2 (running RaxxML and MrBayes) to generate ML and Bayesian subfamily trees.

