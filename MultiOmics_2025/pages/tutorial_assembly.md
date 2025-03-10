# MultiOmics Course : Assembly tutorial

## Learning goals of this tutorial:

We have  number of goals for this course:
* Understand and recognise basic datatypes involved in metagenomics (fasta, fastq, sam/bam...)
* Understand the differences between the types of molecular data (reads, paired-reads, contigs, genomes...)
* Understand and replicate the presented workflows in metagenomics
* Familiarise yourself with manuals and online ressources

## SETUP

We will use `ssh` (Secure SHell) to connect to this cluster, use you ssh-client of choice. To start, we have to log into the login node. Hopefully you set it all up with the keys.

```bash
ssh  <username>@dardel.pdc.kth.se
```

The first time you do this (which should not be the case this time...) there might be a few questions, just answer the default (yes usually).

If you couldn't set up keys (for example because you don't have a SUPR account for reasons), use kerberos

```bash
kinit -f <username>@NADA.KTH.SE
ssh  <username>@dardel.pdc.kth.se
```

kinit will ask for your password, enter it. When you type your password nothing will be shown on the screen. This a security feature, not a bug.

**Let the etherpad know that you got here. **

Let's navigate to everyones working folder, we do not have enough space in our homes to do the work so we will work on the so-called `scratch`-partition, a shared 'hard-drive' that everyone can use, but is occasionally emptied, so be sure to copy anything you wanna keep back to a safer place:

```bash
cd $SNIC_TMP
pwd 
```

`pwd` will show you the current folder, it should look something like `/cfs/klemming/scratch/m/<username>`. The `$SNIC_TMP`-folder should be accessible by you from every computer in the cluster, and has near infinite space (well, 1.1 Petabytes as of right now, you can check it with `du -h .`).

We could use the module system to run all the programs, but we are gonna be a bit lazy and use something called conda. Let's download it and run the installer with :

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

agree to the license, we can use the default location for the install, but bw aware that the conda-folder can baloon in size when used heavily, and the amount of space in your home is limited.

It will run for a little while and copy files and then ask something ("do you wish to update your shell profile").

Now log out of dardel and in again.

Your terminal should now have a `(base)`-prefix, if it does **update your progress on the etherpad**.

We need to set up a few things

```bash
# install mamba to make things faster
conda install mamba

# making sure we don't get sued for using stuff without the right license
conda config --remove channels defaults
rm ~/miniconda3/.condarc

# install all the things
mamba  env create -f /cfs/klemming/scratch/m/morbu/MultiOmics_2025/other_files/env.yaml

# starting the env
conda activate MetaOmics

# testing if it works
megahit -v

```

You should get `MEGAHIT v1.2.9` or something like that if you get there **check point on the etherpad again**.

Now lets connect to one of our reserved nodes:


```bash

# book a part of a node from the reservation
salloc -t 30:00:00 -c 16 -A edu25.slu -p shared --res edu25-slu-2025-03-11

# check the name of the node
squeue -u <username> | grep shared

# connect to the node and start the environment again
ssh nid<node number>
conda activate MetaOmics
```

So now we are ready to start


```bash
COURSE_ROOT=/cfs/klemming/scratch/m/morbu/MultiOmics_2025
OWN_ROOT=$COURSE_ROOT/work_folder/$USER
DB_PATH=$COURSE_ROOT/dbs/
mkdir $OWN_ROOT
cd $OWN_ROOT

```

Now you have a nice little folder in our course folder, just for you, and only you can write in it.

We have prepared a dataset for this course, the libraries can be found at `$COURSE_ROOT/mock_libraries`


```bash

ls $COURSE_ROOT/mock_libraries 

```

This will show you all the `FASTQ`-files we have, they are paired, so there are two per sample. I have prepared a little google doc with a line per sample, you can find it [here](https://docs.google.com/spreadsheets/d/1J0VG6eFK8hgaSmQ96EjCKVV2Ya_OFDqENx_ZHu7O1Lo/edit?usp=sharing), it has some metadata and the paths to the files again. 

I would like each of yous to pick one of the samples on the file and put in your name in the appropriate column of the file, and then you can make a copy of the file to your folder.

```bash

MY_SAMPLE=<sample_name>
cp $COURSE_ROOT/mock_libraries/${MY_SAMPLE}_R1.fastq.gz .
cp $COURSE_ROOT/mock_libraries/${MY_SAMPLE}_R2.fastq.gz .

```
## Read based analysis

### Some QC

Before we do any analysis we  will quality check and filter the reads some, so we gona use `fastQC for an overview of read quality and fastp` to remove bad reads chop off adapters and some more clean up.

```bash
fastqc ${MY_SAMPLE}_R1.fastq.gz ${MY_SAMPLE}_R2.fastq.gz
fastp -w 16 --in1 ${MY_SAMPLE}_R1.fastq.gz --in2 ${MY_SAMPLE}_R2.fastq.gz --out1 ${MY_SAMPLE}_clean_R1.fastq.gz --out2 ${MY_SAMPLE}_clean_R2.fastq.gz --json ${MY_SAMPLE}_fastp.json --html ${MY_SAMPLE}_fastp.html
fastqc ${MY_SAMPLE}_clean_R1.fastq.gz ${MY_SAMPLE}_clean_R2.fastq.gz
```

Get the output files and have a look at them and discuss in your group. [insert the scp part here] **When you bored of them check-in on the etherpad**


### First Taxonomical profiling

First we can have a little look at what is in the samples directly, we will used kraken2 for this:

```bash

# if it's the first time you you run kraken, you will want to run this
sh /cfs/klemming/home/m/morbu/miniconda3/envs/MetaOmics/opt/krona/updateTaxonomy.sh

# classify reads
k2 classify --db $DB_PATH/k2_db/  --output ${MY_SAMPLE}.kraken --use-names --threads 16 --report ${MY_SAMPLE}.kreport  ${MY_SAMPLE}_R1.fastq.gz ${MY_SAMPLE}_R2.fastq.gz

# make a pretty figure
ktImportTaxonomy -t 5 -m 3 -o ${MY_SAMPLE}.html ${MY_SAMPLE}.kreport

# do some transformation to make a useful file
bracken -d $DB_PATH/k2_db -i ${MY_SAMPLE}.kreport -o ${MY_SAMPLE}.bracken -r 150  -l G t 20

```

3 files are worth looking at here. The `.kraken` has a read-wise annotation, were each line is a read with it's classification, which can be useful get get reads matching specific taxa of interest. You can look at it like this:

```bash
less -S ${MY_SAMPLE}.kraken

# similarly for the braken-file

less -S ${MY_SAMPLE}.bracken
```

(leave less by pressing `q`). Braken takes the kraken-report and summarises it at a specific taxonomic level, which is useful to make abundance tables, we summed the data to Genus-level. Put the path to this file into the google-doc, so we can later maybe make one of these.
The last  one is an `html`-file that we can look at with the browser on you computer, for thaat we will have to download it. My way of choice is using scp, your milage might vary.

```bash

#from you own computer

scp <dardel-username>@dardel.pdc.kth.se:/cfs/klemming/scratch/m/morbu/MultiOmics_2025/work_folder/<dardel-username>/<your-sample>.html <wherever-you-want-to-put-your-file>

```

Make sure the values between `<>` are replaced by the appropriate values (copy and paste is your best friend, just saying...). You can then open the file on your computer. From the outputs find the percentage of classified reads and
add it to the google-sheet. ** Once your wrote the stuff into the google-sheet, check in on the etherpad **

### Comparing samples

We can compare these taxonomic profiles to each other [optional task for someone who has too much time/energy, get the braken results from all the students and turn it into an abundance table, maybe even make some plots]. An other option is to use * min-hashes *, it's a fancy-ish computer science thing developed to compare texts and used very effitiently do detect plagiarism. Turns out, it works very well to compare genomic data too, and get similarity scores.

Compute the min-hash of your sample.

```bash
sourmash sketch dna -p k=31,scaled=1000,abund --merge $MY_SAMPLE ${MY_SAMPLE}_R1.fastq.gz ${MY_SAMPLE}_R2.fastq.gz -o ${MY_SAMPLE}.sig
```

Then ask your neighbour where their signature file is, and use sourmash to compare the two:

```bash

sourmash compare <your_signature_file> <your_friends_file>

```

Put the path to your signature file into the google-doc, so if there an overzealous student, they can do a similarity matrix. Anyhow **check into the etherpad** if you here


## Genome-resolved Metagenomics

Let's move on.

### Assembly

Let's assemble the reads into longer contiguous pieces of DNA (contigs). Assembling (meta)genomes can take numerous steps and numerous tries. The usual way is to perform a variety of assemblies using different tools and combinations of samples, and to evaluate them to find the ones that seam more useful/good.

No one solution fits all, but we will try to look at some of the way to do and to judge assemblies.

For the purpose of this tutorial we will only use the `megahit`-assembler. This is mostly due to practical considerations, it is both fast and memory-efficient, additionally it actually is a very good assembler. However it does not necessarily produce the best assembly, so for your personal project consider using other assemblers on as well as `megahit` (for example `metaSPAdes`).

We will start with simple assemblies of our libraries, meaning that we consider our libraries all independent and assemble them as such.

```bash
megahit -1 ${MY_SAMPLE}_clean_R1.fastq.gz -2 ${MY_SAMPLE}_clean_R2.fastq.gz -t 20 -o ${MY_SAMPLE}_assembly
```

Now this is done, you should have an assembly `FASTA`-file in the `${MY_SAMPLE}_assembly`-folder, called `final.contigs.fa`.

We should get some data on it. otherwise it's just a bunch of `ATGC`, well at least before we do some annotation but that is tomorrow.

So let's some stats on the assembly with `quast`:

```bash

quast -o ${MY_SAMPLE}_quast ${MY_SAMPLE}_assembly/final.contigs.fa -t 20

```

There are some fancy `html`-output, but for us now the `report.tsv` should be enough to get an idea and put some data into the sheet.

Now we got an assembly, we want to cluster the contigs in there into Metagenome Assembled Genomes (MAGs), we'll use metabat, because. But first we need to map our sample back to the assembly, and then use a script to jumble it around.

```bash

bowtie2-build ${MY_SAMPLE}_assembly/final.contigs.fa ${MY_SAMPLE}_idx
bowtie2 -x ${MY_SAMPLE}_idx -1 ${MY_SAMPLE}_clean_R1.fastq.gz -2 ${MY_SAMPLE}_clean_R2.fastq.gz -S ${MY_SAMPLE}_bowtie2.sam --threads 20
samtools view  -b -S -@20  ${MY_SAMPLE}_bowtie2.sam | samtools sort -@ 24 -o ${MY_SAMPLE}_bowtie2.bam -

jgi_summarize_bam_contig_depths --outputDepth ${MY_SAMPLE}_depth.tsv  ${MY_SAMPLE}_bowtie2.bam 
metabat2 --inFile ${MY_SAMPLE}_assembly/final.contigs.fa --outFile ${MY_SAMPLE}_bowtie_bins/${MY_SAMPLE}_bin -a ${MY_SAMPLE}_depth.tsv  -t 20

```

Hopefully you got at leaast a bin out (I cheated, made sure all have a bin at least)

Now we'd still like to know what those bins are and how good they are.

```bash

export GTDBTK_DATA_PATH=$DB_PATH/gtdbtk/release220
gtdbtk classify_wf  --genome_dir ${samp}_bowtie_bins -x .fa --out_dir ${samp}_classification --cpus 20 --skip_ani_screen
checkm lineage_wf -x .fa  ${samp}_bowtie_bins/  ${samp}_checkm/ > ${samp}_checkm.txtb


```





#### QC the assembly

A variety of tools is available to QC genomic assemblies, however many of them are not appropriate for metagenomic data, as they presuppose a single organism. However for most metagenomic applications, a few factors will already be indicative of the quality of the assembly.

One of the first factors is the amount of reads that align to the assembly. The proportion of reads mapping our assembly gives us a good idea of how representative of the sampled environment the assembly is.

To get this information we will need a mapping tool. We will use `BBmap`. Mappers, like `BBmap`, don't only return the proportion of mapping reads, but also the specific alignments, these are output as SAM-files, a particularly horrible and wasteful file type.

> Using `bbmap.sh` loaded with `module`, use the `bamscript=`-option
> Map your library to your assembly. What proportion of your reads map to your assembly? Any comments about the size of the SAM-file?
> Optional : map subsets of the other samples to your assembly, how does it look like? (checkout the `reads` option, and also the  `showprogress` one , mostly because I like it)

SAM-files are HUUUUUGE, we do not like that, BAM-files, the compressed version of SAM is a bit better.

> Use the script written by the `bamscript=` option of `bbmap.sh`.
> What are the generated files?
> Keep the BAM-file and other file, we will need them later.

However this does not say anything about the quality of the `contigs` directly, only about how much of the data has been assembled, but it could all just have assembled into very small contigs!

We need to look at how the contigs "look-like". A first simple way of looking at this is to simply do some stats on contig lengths and coverage. There are many ways probably to do this, but we will use a script that has been written by yours truly. We will use this to learn a bit about scripting and stuff.

You can find the script in [here] (../scripts/seq_stats.R), this is a pretty simple `R`-script, make a copy in your working directory.

To run an `R`-script, you can just do:

 ```bash
Rscript my_script.R
 ```

 Let's try to run the script you just obtained from a friendly collaborator.

 > run the script

Hmmm, you probably got errors. Many scripting languages have systems to load additional packages of functions. For R you have to use the `install.packages` function to get packages from CRAN, also there is something called `bioconductor`.

> Open R
> install the missing packages, try to run the script

Probably you still got errors ... Well it's probably time to open the script and have a look.

> Fix the script, the error messages obtained should give you an idea of the problem
> Now it should run, you should get some basic info about your assembly, as well as a pretty plot
> Optional : Personalize the script a bit. And run the script with all the assemblies.

We have now some basic data about our few simple assemblies.

> Any particular feelings about the obtained assemblies? Also, what about the plot?!

One other way to check out an assembly is to inspect the k-mer graph. Bandage is a cool tool to do that!

> Install Bandage and use [this] (https://github.com/voutcn/megahit/wiki/Visualizing-MEGAHIT's-contig-graph) to inspect your assembly!

#### Alternate assembly.

So we are "lucky" our assemblies look kinda nice ... we can now make some other assemblies. We could normalize the reads to remove micro-diversity and errors. We could split the reads into different subsets (e.g. eukaryotic and prokarotic for example). Or coassemble multiple libraries together.

> Discuss some strategies
> OPTIONAL : make an other assembly.
