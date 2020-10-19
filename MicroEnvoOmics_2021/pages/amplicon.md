```
CHECKPOINT
```

## AMPLICON-ANALYSES

### Tools

General amplicon processing and more:
* QIIME
* MOTHUR
* usearch (used in this course, mostly because of convinience)
* dada2 (bioconductor)

Quantification:
* phyloseq (bioconductor)
* metagenomeseq (bioconductor)
* DEseq2 (bioconductor)

Other:
* PiCRUST
* Phylotyping

### Processing the amplicon

Most analyses of 16s rRNA gene amplicons rely on an OTU-table, e.g. a table of counts of microbial-types in your samples. For this we need to cluster the reads into 'OTUs' (Operational Taxonomical Units), these are then counted.

We will use the `usearch` software to do the necessary computations (https://www.drive5.com/usearch/manual/). We are going to be fancy and we are going to use the latest version of `usearch` which is not on UPPMAX. We have it however locally in our project folder in `/proj/g2017026/bin/` (as well as some other software we will use later). For ease of use we will put this folder in our `PATH` environment variable (this contains all the folders which have executables readily available).

```bash
PATH=/proj/g2017026/bin/:$PATH
```

> Optional : which other folders are in the `PATH`

From now one for a few exercise you will pic one of the amplicon libraries to process. We will call it MY_LIBRARY for the next few questions!

#### Merging reads

As we saw, each library has pairs of reads. In the case of a well designed amplicon the PCR product should be shorter than the lengths of the two reads added, in which case they will overlap, and the two pairs can be merged into one longer pseudo-read which will improve the quality of the analyses. This is the case for this dataset.

> Using the `usearch10` program with the `-fastq_mergepairs` option
>
> Merge the pairs of MY_LIBRARY
>
> Optional : play with the parameters to increase/decrease merging rate

<!-- One sample is nice and fine, but for a good analysis of the data we will need to analyze all the reads as one. If we just concatenate the files we will lose the information of which samples the reads came for. Luckily `usearch` can make our life a bit easier!

> Using the `--relabel` and `-fastq_mergepairs` options of `usearch10`
>
> Merge all the libraries.
> How many reads got lost in the merging
>
> Optional : any relationship between the merging rate and some of the sampling-parameters -->

#### QC-ing the reads

The reads have now been merged into pseudo-reads, but we have not done any quality-check yes. For merged paired reads it is convenient to do it after merging, as the merging can increase the quality of the overlapping bases! In all sequencing runs there will be some bad reads, this can be due to many reasons, but usually mostly to the quality of the DNA, e.g. extraction and PCRs, or interactions between spots on the sequencing array. We will use the quality scores from the fastq files to remove reads below as certain quality.

> Using the `-fastq_filter` option of `usearch10`
>
> Remove from MY_LIBRARY the reads with more than one expected error and output it as a FASTA-file
>
> Optional: How does the loss of reads vary with the expected error.


#### Dereplicating

All clustering algorithms can be sensitive to the number of reads, however amplicon data has huge amounts of duplicates due to the PCRs, so one way to help the clustering is to remove these duplicates. Also, as there it is very reasonable to expect duplicates, any reads that appear only once is probably due to either PCR- or sequencing errors. In `usearch` both of these steps are done together.

> Using the `-fastx_uniques` option of `usearch10`
>
> Dereplicate MY_LIBRARY and remove singletons

> Optional: How does the count distribution of the reads look like?

#### Clustering

We are finally ready to do the clustering. The clustering heuristic of `usearch` is pretty fast and incorporates chimera detection as well. Please note that one of the keys of the heuristic is that the reads are aligned, as they are in an illumina amplicon, but might not be in other cases where you'd like to run this tool for!

> Using the `-cluster_otus` option of `usearch10`

> Compute 97% similarity sequence-clusters in MY_LIBRARY

> Optional: vary the identity cutoff, also what happens if you run the clustering without removing the singletons in the previous step

#### Mapping back

Now the whole data of MY_LIBRARY has been reduced to a small number of sequence-clusters (OTUs). Each of them represented by the most abundant read of that cluster. However the information of the total amount of reads has been lost along the way. For this reason we need to map the reads back to the OTUs, which means we need to find for each read the OTU that matches it best and that is more similar then the identity threshold. This will also allow us to recruit back reads lost during QC!

> Using the `-usearch_global` option of `usearch`
> Get some help from : https://www.drive5.com/usearch/manual/cmd_otutab.html
> Map the merged reads to the representative sequences of the clusters, and make an otu-table

> Optional : what are the reads that did not map?

```
CHECKPOINT
```


#### Putting it all together

Well, this is all nice and fine,  but this is only one sample, so not really a table, just a vector. We could run the same thing for every sample separately but we would just get many separate vectors, we could not connect the OTUs easily to each other. What we need to do is to run the clustering step with all the reads of all the sample at the same time!

> Using all the options of `usearch10` that you have used until now,
> Have a look at the `-relabel` option of the `-fastq_mergepairs` function

> Make an OTU-table with all the samples, e.g. run the merging, QC for every sample separately, concatenate the resulting FASTA-files, and follow then with the rest of the script.

> Optional : quantify the read-losses. Are there any experimental factors influencing the loss of reads?



#### New Age Clustering

In recent years there have been many discussions about OTUs and there meaning and the issues linked to them. Hence new modern methods of clustering amplicon reads have been developed. Much talked about recently has been the `dada2`-tool running with `R`, this tool does not use sequence identity to cluster reads, but build a noise model to clean reads from variation that is not biologically relevant. `usearch` also released recently a similar heuristic, the `unoise` algorithm.

> Using the `-unoise3` option of `usearch10`

> redo the OTU-table

> Optional : compare the mapping rates between the two methods


#### Taxonomical Annotation

We have not a table with OTUs and their abundances in each samples, as well as representative reads for all the OTUs. To learn more about our data we want also to annotate the reads. In the case of amplicon, we want to annotate for taxonomy, as we normally know exactly what genes we already have. For 16s rRNA amplicon there are many tools and databases. For sake of simplicity we will use a tool also provided by the `usearch` software, the `sintax`-tool.

> Let's split the class in two using the two OTU-tables

> Using the `-sintax` option of `usearch10`
> and the taxonomy/16s database found in /proj/g2017026/data/rdp_16s_v16_sp.udb

> Give all the OTUs taxonomical classifications.

> Optional : play with the cutoff, and try to install/use an other database.

### Let's make some plots : or where the real game starts!

We won't really have time to dwelve into diverse quantification methods and statistaical analyses, that would need many more prerequisites and much more time.

However, often in bioinformatics you get to use other people scripts to do a variety of tasks, they might come in a variety of languages. In this case I provide you with an R-script to make a little plot with this! It is called `ampliconplot.R` and you can find it in the `scripts`-folder.
`R` is already install on rackham so to run an `R` script, you normally just need to do `Rscript name_of_script.R`.

> Have a look at the `R`-script
> Try to run it!

You might notice that it complains about an abscent package. This is common when running scripts of other people. However, installing scripts in R is pretty easy. The easiest is to open `R` in the terminal, e.g. type `R`. You will see a prompt, and now run the command:

```R
install.packages(“reshape2”) # installs the package
quit() # leaves R
```

### Solution

```bash
MY_PATH=/proj/g2017026/2017_MG_course/home/$USER/amplicon

cd $MY_PATH

for sample in `find /proj/g2017026/2017_MG_course/raw_data/amplicon/ -name "*_R1*.fastq"`;
do
  id=`basename ${sample%%_R1.fastq}`
  out_file=$MY_PATH/$id.fastq
  usearch10 -fastq_mergepairs $sample -fastqout $out_file --relabel "$id:"
  cat $out_file >> $MY_PATH/all_reads.fastq
done

### Question can you say something about the merging rate?! improve the rate maybe?!

usearch10 -fastq_filter all_reads.fastq -fastq_maxee 1 -fastaout all_reads.clean.fasta
usearch10 -fastx_uniques all_reads.clean.fasta -relabel Uniques_ -sizeout --minuniquesize 2 -fastaout all_reads.uniques.fasta
usearch10 -cluster_otus all_reads.uniques.fasta  -relabel OTUs_ -otus  all_reads.OTUs.fasta
usearch10 -usearch_global all_reads.fastq -db all_reads.OTUs.fasta -strand plus -id 0.97 -otutabout OTU_table.txt

### How do hit rates look for specific libs?!

DB=/proj/g2017026/2017_MG_course/data/usearch/rdp_16s_v16_sp.udb

usearch10 -sintax all_reads.OTUs.fasta -db $DB -tabbedout taxonomy.tax -strand both -sintax_cutoff 0.8
usearch10 -unoise3 all_reads.uniques.fasta -zotus all_reads.zOTUs.fasta
usearch10 -usearch_global all_reads.fastq -db all_reads.zOTUs.fasta -strand plus -id 0.97 -otutabout zOTU_table.txt
usearch10 -sintax all_reads.zOTUs.fasta -db ~/Data/tax_dbs/rdp_16s_v16_sp.udb -tabbedout ztaxonomy.tax -strand both -sintax_cutoff 0.8

Rscript ~/MG_course/2017_MG_course/scripts/ampliconplot.R

```

Answer all the questions with yes.

> Try running the script again
> Add any other package possibly needed.
> Make sure the filenames match!
> Edit the script if needed, but make a personal copy first!
>
> Optional : Why did I make separate plots for nasal and fecal data?
>            Make barstacks with other taxonomical levels
