## Data-types

### FASTQ-files

Sequencing facility might provide data from different technologies in a variety of formats, but most probably one of the formats the data will be delivered in will be the FASTQ format. We will start this tutorial with data in this format. It is a wide text-based format for sequences, they are related to the FASTA format which will see too.

A FASTQ-file can contain many many sequences, each sequence is represented by two data-types, the actual sequence and the quality scores:

```
@an_example_sequence_in_a_fastq_file
ACACATATATACACACATACACACACATATACACACACATATAAACACATATATACATTTATATGCATATATTAATACATATATATTTAAGTTGATGGAGAGTATAACAGAGTTAGGCTGCTTATT
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFF/FFFBFFBFFFFFFFFFFFFFFFFFFFFFFBFFFF//FFFFFFFF<FFFF/FFFFFB/FFF<BF//<</BFFFFF
```

Normally each entry in a FASTQ-file is composed of exactly four lines.

The quality scores are related the quality of sequencing. They are `Q = -10 log<sub>10</sub> p` where `p` is the probability of a wrong base-call, and each letter corresponds to an integer value of Q. These letter can differ a bit depending on the sequencing platform, for recent illumina we have this correspondence table (courtesy of wikipedia):

```
  BCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi
   |     |                               |
   3.....9...............................41
 ```

 ### FASTA-files

 FASTA-files are the mother of all sequence data-types. It is a very simply text base file type compose of one or multiple entries. Each entry is composed of two type parts: A description part, which is a single line starting with `>` containing text information about the following sequence; and a data-part containing sequence data as plain text, this can be spread over multiple lines. FASTQ-files can be easily transformed to FASTA-files by removing the qualities and the `+`-line and replacing `@` by `>`, e.g:

 ```
 >an_example_sequence_in_a_fasta_file
 ACACATATATACACACATACACACACATATACACACACATATAAACACATATATACATTT
 ATATGCATATATTAATACATATATATTTAAGTTGATGGAGAGTATAACAGAGTTAGGCTG
 CTTATT
 ```

 ### SAM/BAM-files

 SAM is a text-based file type used to represent alignments, usually used for NGS-mappers, when one aligns many short reads to a reference.

 It is normally composed of a header indicated by lines starting with the `@` character (which also indicated comment lines, e.g. information that does not correspond to the SAM-format), followed by often many many many lines, one per sequence to be aligned. Each line is split into columns with tabulations and each valuer corresponds to some info. From wikipedia:

 | Col | Field | Type | Brief Description |
 | --- | ----- | ---- | ----------------- |
 |1|QNAME|String|Query template NAME|
 |2|FLAG|Int|bitwise FLAG |
 |3|RNAME|String|References sequence NAME|
 |4|POS|Int|1- based leftmost mapping POSition|
 |5|MAPQ|Int|MAPping Quality|
 |6|CIGAR|String|CIGAR String|
 |7|RNEXT|String|Ref. name of the mate/next read|
 |8|PNEXT |Int|Position of the mate/next read|
 |9|TLEN|Int|observed Template LENgth|
 |10|SEQ|String|segment SEQuence|
 |11|QUAL|String|ASCII of Phred-scaled base QUALity+33|

 or for meore detail : [SAM/BAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

 These files can be massive so they are often compressed using `samtools` (available trough `module` on UPPMAX), and are then often called BAM-files. An example :

 ```
 2167:2:1101:11949:2366  16      contig_99_631734_length_1867_multi_8_in_2_out_1 1493    42      180M    *       0       0       ACCTTACTAGTCGCTGTTGCGAAAGATAATTTTTAGTTGCTAATAATCGGAGCCCTGCATGGCCATCACTTTATACCCTGTTAAACCAGATTTTATGGCTGAAATCGGTGATGTCGACTTGAGTATGCCGCTCTCTGACGAAGACCGCGCCGCCATAAAAGCCGCTTTCTTTAAATACAG    bbbeeeeegggggiigfhiiifhhiiihiiiiiiiiiiiiiiiiiiefhiefhihiiggggfeeeeebbbdddddccbcccccbccccaaccbb`bccccccab`bb_^^^bc_^`^bba`bbbaaaa`aa]]`bc_^ZZb`^]edfffiiiiiiiiiiiiiiiiiigggggeeeecbbb    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:180        YT:Z:UU
 2167:2:1101:11897:2252  4       *       0       0       *       *       0       0       TTTAAATTTCTACATACATACTTTTCTAAATAAACCTGGACGCAATTTAAAAATTTCAAATTAATTTTATTTAATAAAAGCTTTTATTAATTTACGAAAGCGTCTAAAATACTAAAAAAAATACAAAAATATTTAAACTGTGTCTTCATTTAAAAAATCCAAATCAGTAAAATGTGAAATTGTATTTTTTTT        ^__ec`ccgggcefg_cdhhaegggdghhhfcb^gae^eghdf#af`cfd^a_eeeehbbbeghae`fccdgfbgdgedb_bbcbddb_cacd`bac_[ca`aadcdcddb`bcaabdeccccb`ZZ^gg`bbbgfgggebefggffffhgffbfccdhgebe_ihfhhhfhfgahhgf]gcgceceec__^        YT:Z:UU
 2167:2:1101:10183:2426  0       contig_99_3998773_length_1078_multi_3_in_1_out_1        849     42      212M    *       0       0       AACTAACCTTAACAAGCCTTAGCAAGATGTATTTTTGCTTTGTAATAAGCTAGAGTGTTTCTTCCTGGCAAACTTTTTCTAACTAGCTATGACTTCGAAGCTTGTTTGGTTCTAAAAACAATGTAAGCATCGTAGCATTTAAAGAACCAATTGTCAATTTTACTACCTCTTAAAAAACAATTAAAGCATCGTAGCATTTAAAGAACCAATTG    bbbeeeeegggfgiiiihhiiiifhiiiihiiiiiihiiihfgihiiihihf^cgghiiiiiiiiiihigghfhhhiggggggeeeeeeddcdddccbb`acccccccb``ccddccccddcccdcbcccceeeeeggggeggcc\hgfhhihhiihiiiihgfahfageeiiihgeihiiiiiiiifgiiiiiiiiihgggggeeeeebbb    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:212        YT:Z:UU
 ```

 ### GFF-files

 GFF-files are text-files that contain annotation information for sequence data. Basically these are tab-delimited tables, where each line corresponds to an annotation. You can get more information [here](https://www.ensembl.org/info/website/upload/gff.html).
