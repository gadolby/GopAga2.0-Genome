# Automating gene family evolution analyses

### Step 1: Clone in GitHub repository at https://github.com/mmoral31/GopAga2.0-Genome.git

```bash
> git clone https://github.com/mmoral31/GopAga2.0-Genome.git
```



### Step 2: Gather protein FASTA to be searched against and create a BLAST database

> Requires BLAST+

if you want to query against whole proteomes, search for the target species on NCBI and pull protein FASTAs (reference: *Chelonia mydas*)

```bash
> makeblastdb -in GCF_000344595.1_CheMyd_1.0_protein.faa -dbtype prot -parse_seqids -out CmydPro
```



### Step 3: Gather query reference sequences in a single FASTA and set up input text file with gene names

determine the query sequences to be used and place into a single FASTA file labelled as {gene family}_ref.fasta where {gene family} is the name given to the group of genes used as a query (reference: TLR_ref.fasta)

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

create a text file with each line containing the name of one gene family member in the order of the reference FASTA (reference: TLR.txt)

```
TLR3
TLR4
TLR7
```



### Step 4: Use blastp.py to extract gene family members from the target species

> Requires BLAST+ and Biopython

blastp.py takes three parameters in this order:

- name of the created BLAST database without suffix (reference: CmydPro)
- name of the gene family for query (reference: TLR). this is important for distinguishing the reference FASTA (reference: TLR_ref.fasta)
- name of the gene family text file (reference: TLR.txt)

```bash
> python blastp.py CmydPro TLR TLR.txt
```

this will output one file per gene family member to record standard output BLAST results and one aggregate file with all pulled sequences from the target species.

- {gene family member}_{db name} for all of the individual records

- {gene family}_{db name}.fa for the aggregate FASTA file

  

### Step 5: Use pull_id_match.py to separate gene family into subfamilies for alignment

> Requires Biopython

if the gene family of interest has separate subfamilies (as in TLRs), use pull_seqids.py to separate and extract a subfamily of interest

pull_id_match.py takes three parameters in this order:

- name of the FASTA file to extract from (reference: TLR_CmydPro.fa)
- name of the text file with listed keywords to pull with (reference: TLR3SF.txt)
- name of the output file where pulled sequences will aggregate (reference: TLR3SF.fa)

the text file will have one keyword on each line. the FASTA output will include any records that matched at least one of the keywords in the list.

```
TLR3
TLR4
```

```bash
> python pull_id_match.py TLR_CmydPro.fa TLR3SF.txt TLR3SF.fa
```



### Step 6: Use run_alignments.py to align respective subfamilies with ClustalW

> Requires ClustalW2

subfamilies need to be aligned in order to proceed with phylogenetic reconstruction. run_alignments.py can be used to align respective subfamilies.

run_alignments.py takes one parameter:

- name of the text file with listed names of FASTA files (without suffix) to be aligned (reference: TLR_align.txt)

```
TLR3SF
```

```bash
> python run_alignments.py TLR_align.txt
```

three alignment files are output in three formats: FASTA, NEXUS, PHYLLIP



### Step 7: Use Prottest to determine the best substitution matrix for maximum likelihood phylogenetic reconstruction

> Requires ProtTest

maximum likelihood phylogenetic reconstruction requires a specified substitution matrix. Prottest can be used to determine the ideal substitution matrix for that data. requires phyllip alignment input.

```bash
> prottest3 -i TLR3SF_alignment.phy -o TLR3SF_models
```



### Step 8: Use the CIPRES Science Gateway with RaxxML and MrBayes to phylogenetically reconstruct respective subfamilies with maximum likelihood and Bayesian methods

