# MultiOmics Course, March 2025

## Learning goals of this tutorial:

We have  number of goals for this course:
* Understand and recognise basic datatypes involved in metagenomics (fasta, fastq, sam/bam...)
* Understand the differences between the types of molecular data (reads, paired-reads, contigs, genomes...)
* Understand and replicate the presented workflows in metagenomics
* Familiarise yourself with manuals and online ressources


## UPPMAX :

All the exercices here will be done using dardel (https://www.pdc.kth.se/hpc-services/computing-systems/dardel-hpc-system/dardel-1.1043529).

The rackham-cluster of uppmax is a distributed HPC cluster, meaning it is a lot of powerful computers networked together. It has 4 login-nodes and many more computing nodes.

### SSH

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

mkdir $OWN_ROOT
cd $OWN_ROOT

```

Now you have a nice little folder in our course folder, just for you, and only you can write in it.

We have prepared a dataset for this course, the libraries can be found at `$COURSE_ROOT/mock_libraries`


```bash

ls $COURSE_ROOT/mock_libraries 

```

This will show you all the `FASTQ`-files we have, they are paired, so there are two per sample. I have prepared a little google doc with a line per sample, you can find it [here](https://docs.google.com/spreadsheets/d/1J0VG6eFK8hgaSmQ96EjCKVV2Ya_OFDqENx_ZHu7O1Lo/edit?usp=sharing). 


We will not spend much time on the login node. we just want to know which computer is booked for each of us.

> Optional question: what does the `-X` in the command do.

Rackham, as most UPPMAX computer systems, runs a queuing tool called SLURM (Simple Linux Utility for Resource Management, (https://slurm.schedmd.com/). This tool is used to distribute computing resources to all users. A computer has been booked for each participant, SLURM knows 'its name'.

On Rackham, to know which resources you have requested (or someone else requested for you) the `jobinfo` command is used.

> Use the jobinfo command to find the name of your dedicated computer

Once you have the name of the computer (something looking like r123), you can log directly onto that computer. That name is the name of the computer in the local network of rackham, so you can simply connect with :

```bash
ssh -X r123
```

You are now on the computer that you can call `$HOME` for the next 3 days (actually only todfay, you will have to do that every day, we have a different booking for each day).

The data we are gonna use for this course is in our common project folder, e.g. `/proj/g2019027/`, more specifically in  `/proj/g2019027/2019_MG_course/raw_data`.

> Using the `ls`-  or even better `tree`-function
> How many fastq-files can you find in this folder, and list all of them?!

As we see, there is plenty of libraries, each of us should pick (at least) one:

> Go to (https://etherpad.wikimedia.org/p/SLU_metagenomics_course_2019) [our virtualwhiteboardthing]
> put your name next to a sample of your choice! This is gonna be your sample!
> What data, and depth was your sample taken at?


### The module system and other software prerequisites

Rackham uses a module system for programs (http://www.uppmax.uu.se/resources/software/installed-software/). So many bioinformatic programs are available without installing anything. We will use this as much as possible during this course. You can list all available module with `module list`. Let's load the `bioinfo-tools`-module, which is a prerequisite for most of the tools we will use in the next couple of days:

```bash
module load bioinfo-tools
```

Modules are loaded on Rackham until you disconnect from your terminal, so one you loaded them you do not need to rerun these commands. But if you disconnect, they will not be loaded anymore, so you will need to reload on you next connect!
## Processing your reads

### QC-ing

 As mentioned, the raw-data we will use is available in the `/proj/g20190272019_MG_course/raw_data`, make sure you name you sample here [https://etherpad.wikimedia.org/p/SLU_metagenomics_course_2019] (our virtualwhiteboardthing). This is unpublished data from a lake in Switzerland, multiple time points as well as depth points. We have already preprocessed and subset it some to guarantee some quality of asssembly, but let's check it out!

 > Use the linux command  `ls`, `wc`, and `du` :
  > What is the content of this folder? How much data do we have?!
  > Optional : Can you guess what the file names mean?!

The first step is to quality check these libraries. From this mornings tutorial with Domenico, you should have learned some tools, it's time to use them!

> Use the tools from the previous tutorial to QC your reads
> What can you learn from the returned report?!
> Can you say something about the number of reads?
> Is your library good?
> Very optional : Anything to say about the diversity of the reads in the libraries?

### Read trimming

Once we have some idea of how problematic the libraries are, and potential issues, we can use some tools to clean the libraries a bit reduce noise and increase quality

A very commonly used tool for this is `trimmomatic`, it can do a large variety of operations on reads.

> Load `trimmomatic` with UPMMAX's `module`-system

> Clean the reads!

> Optional : re-run `FastQC` on the cleaned libraries and discuss

## Genome-resolved Metagenomics

### Assembly

Let's assemble the reads into longer contiguous pieces of DNA (contigs). Assembling (meta)genomes can take numerous steps and numerous tries. The usual way is to perform a variety of assemblies using different tools and combinations of samples, and to evaluate them to find the ones that seam more useful/good.

No one solution fits all, but we will try to look at some of the way to do and to judge assemblies.

For the purpose of this tutorial we will only use the `megahit`-assembler. This is mostly due to practical considerations, it is both fast and memory-efficient, additionally it actually is a very good assembler. However it does not necessarily produce the best assembly, so for your personal project consider using other assemblers on as well as `megahit` (for example `metaSPAdes`).

We will start with simple assemblies of our libraries, meaning that we consider our libraries all independent and assemble them as such.

> Using the `megahit`-assembler
> Assemble your the reads of your sample!
> Optional : How would you run a different assembly?


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
