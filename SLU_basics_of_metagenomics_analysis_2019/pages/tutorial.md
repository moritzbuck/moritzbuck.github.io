# Metagenomics Course, Fall 2019

## Learning goals of this tutorial:

We have  number of goals for this course:
* Understand and recognise basic datatypes involved in metagenomics (fasta, fastq, sam/bam...)
* Understand the differences between the types of molecular data (reaads, paired-reads, contigs, genomes...)
* Understand and replicate the presented workflows in metagenomics
* Familiarise yourself with manuals and online ressources


## UPPMAX :

All the exercices here will be done using UPPMAX (https://www.uppmax.uu.se/) more specifically the Rackham cluster.

The rackham-cluster of uppmax is a distributed HPC cluster, meaning it is a lot of powerful computers networked together. It has 4 login-nodes and many more computing nodes.

### SSH

We will use `ssh` (Secure SHell) to connect to this cluster, use you ssh-client of choice. To start, we have to log into the login node.

```bash
ssh -X MY_USERNAME@rackham.uppmax.uu.se
```

The first time you do this (which should not be the case this time...) there might be a few questions, but eventually you will be asked for your password, enter it. When you type your password nothing will be shown on the screen. This a security feature, not a bug.

We will not spend much time on the login node. we just want to know which computer is booked for each of us.

> Optional question: what does the `-X` in the command do.

Rackham, as most UPPMAX computer systems, runs a queuing tool called SLURM (Simple Linux Utility for Resource Management, https://slurm.schedmd.com/). This tool is used to distribute computing resources to all users. A computer has been booked for each participant, SLURM knows 'its name'.

On Rackham, to know which resources you have requested (or someone else requested for you) the `jobinfo` command is used.

> Use the jobinfo command to find the name of your dedicated computer

Once you have the name of the computer (something looking like r123), you can log directly onto that computer. That name is the name of the computer in the local network of rackham, so you can simply connect with :

```bash
ssh -X r123
```

You are now on the computer that you can call `$HOME` for the next 3 days (actually only todfay, you will have to do that every day, we have a different booking for each day).

The data we are gonna use for this course is in our common project folder, e.g. `/proj/g2019027/`), more specifically in  `/proj/g2019027/2019_MG_course/raw_data`.

> Using the `ls`-  or even better `tree`-function
> How many fastq-files can you find in this folder, and list all of them?!

As we see, there is plenty of libraries, each of us should pick (at least) one:

> Go to (https://etherpad.wikimedia.org/p/SLU_metagenomics_course_2019)[our virtualwhiteboardthing]
> put your name next to a sample of your choice! This is gonna be your sample!
>
> What data, and depth was your sample taken at?


### The module system and other software prerequisites

Rackham uses a module system for programs (http://www.uppmax.uu.se/resources/software/installed-software/). So many bioinformatic programs are available without installing anything. We will use this as much as possible during this course. You can list all available module with `module list`. Let's load the `bioinfo-tools`-module, which is a prerequisite for most of the tools we will use in the next couple of days:

```bash
module load bioinfo-tools
```

Modules are loaded on Rackham until you disconnect from your terminal, so one you loaded them you do not need to rerun these commands. But if you disconnect, they will not be loaded anymore, so you will need to reload on you next connect!

## Processing your reads

### QC-ing

 As mentioned, the raw-data we will use is available in the `/proj/g2019027/raw_data`. This is unpublished data from a lake in Switzerland, multiple time points as well as depth points. We have already preprocessed and subset it some to guarantee some quality of asssembly, but let's check it out!

 > Use the linux command  `ls`, `wc`, and `du` :
 >
 > What is the content of this folder? How much data do we have?!
 >
 > Optional : Can you guess what the file names mean?!

The first step is to quality check these libraries. From this mornings tutorial with Domenico, you should have learned some tools, it's time to use them!

> Use the tools from the previous tutorial to QC your reads
>
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

#### Assemble

We will start with simple assemblies of our libraries, meaning that we consider our libraries all independent and assemble them as such.

> Using the `megahit`-assembler

> Assemble your the reads of your sample!

> Optional : How would you run a different assembly?


#### QC the assembly

A variety of tools is available to QC genomic assemblies, however many of them are not appropriate for metagenomic data, as they presuppose a single organism. However for most metagenomic applications, a few factors will already be indicative of the quality of the assembly.

One of the first factors is the amount of reads that align to the assembly. The proportion of reads mapping our assembly gives us a good idea of how representative of the sampled environment the assembly is.

To get this information we will need a mapping tool. We will use `BBmap`. Mappers, like `BBmap`, don't only return the proportion of mapping reads, but also the specific alignments, these are output as SAM-files, a particularly horrible and wasteful file type.

> Using `bbmap.sh` loaded with `module`, use the `bamscript=`-option
>
> Map your library to your assembly. What proportion of your reads map to your assembly? Any comments about the size of the SAM-file?
>
> Optional : map subsets of the other samples to your assembly, how does it look like? (checkout the `reads` option, and also the  `showprogress` one , mostly because I like it)

SAM-files are HUUUUUGE, we do not like that, BAM-files, the compressed version of SAM is a bit better.

> Use the script written by the `bamscript=` option of `bbmap.sh`.
>
> What are the generated files?
>
> Keep the BAM-file and other file, we will need them later.

However this does not say anything about the quality of the `contigs` directly, only about how much of the data has been assembled, but it could all just have assembled into very small contigs!

We need to look at how the contigs "look-like". A first simple way of looking at this is to simply do some stats on contig lengths and coverage. There are many ways probably to do this, but we will use a script that has been written by yours truly. We will use this to learn a bit about scripting and stuff.

You can find the script in [script-folder], this is a very simple python script, make a copy in your working directory.

To run a python script, normally you'd just have to do:

 ```bash
 python my_script.py
 ```

 Let's try to run the script you just obtained from a friendly collaborator.

 > run the script

Hmmm, you probably got an error. Many scripting languages have systems to load additional packages of functions. The system used in python is called `pip`. So you can install packages for `python` doing:

```bash
pip install my_package
```

However this normally needs superuser permission, so you will have to run this with the `--user` option, which will install the package into your `home`

> Using `pip`

> install the missing packages and run the script, use it to get some basic info about your assembly

> Optional : Personalize the script a bit. And run the script with all the assemblies.

We have now some basic data about our few simple assemblies.

> Any particular feelings about the obtained assemblies?

Some of the assemblies as you can see have decent mapping rate, however, the contigs composing them are not very good. Often, mostly if sequencing depth, things assemble well, but the lack of coverage breaks the assembly in many places. Also certain assemblers output also very small contigs which are often not usable for downstream analyses, but give an impression of high-mapping rates. Hence, assemblies are often length filtered. I am sure there are plenty of tools to do that, we however provided a home made script again, it also cleans up the names of the contigs a bit.

> Using the `fasta_filter.py` script

> Filter out small contigs from your assembly

> Optional : ammend the script to give personalised names to your contigs

Now, we just have to see how the characteristics of this 'new'-assembly are.

> Re-compute mapping rate and metagenome-characteristics


#### Alternate assembly.

So some assemblies are OK, some are not ... we can now make some other assemblies. We could normalize the reads to remove micro-diversity and errors. We could split the reads into different subsets (e.g. eukaryotic and prokarotic for example). Or coassemble multiple libraries together.

> Discuss some strategies, and make an other assembly.
> Run the mapping and the `python`-script once it is done


### Binning the contigs

We are more interested in extracting genomes in our analysis, so we will use a binning approach with our assembly. To extract genomic bins, or Metagenome Assembled Genomes (MAGs), or Genomes From Metagenomes (GFMs), we will use the `MetaBAT` tool.

This tool will need the sorted BAM-files we made to run.

> For you genomic-assembly of choice. Collect sorted indexed BAM-files corresponding to that assembly in a folder. (Remember the `bamscript=` option of `bbmap`)

> Using `runMetaBat.sh`, loaded with `module load MetaBaT`. Have a look at the FASTA-files generated.

> Optional : run the binning on an other assembly!

##### Evaluating bins

Now we each hopefully got a number of bins. Meaning that the assembly has been split into bags of contigs that somehow look/vary similarly. This process is however far from perfect, some bins will be-merged MAGs, others will contain viral data, or Eukaryotic chromosome, or simply trash.

We need to evaluate the bins for that. A very simple first filter is to remove bins that are too small or too large.

> compute rough size estimates of your bins (hint, the file size correlate well with the genome size), and pic a bin that has a size more or less corresponding to a microbial genome.

> optional : compute other stats on these bins, let your imagination run wild (or maybe use the script I provided earlier)

Now we should all have a bin we hope is a genome! To check this we will use a tool called `checkm`, fundamentally it uses single copy marker genes to check if a bin/MAG is complete, and if it is contaminated.

> Load `checkm`. You will need once again UPPMAXs `module` system.

Once it seams like it is running we want to run the `lineage_wf`, the documentation of `checkm` is a bit more confusing than other tools... So I advice you to look at [this specific page](https://github.com/Ecogenomics/CheckM/wiki/Quick-Start) and don't get lost in the others!

> Using the `lineage_wf` of `checkm`

> Compute completeness/contamination for your MAG of choice. Did you get any bonus information?

> Optional : compute this for more MAGs, and try the [taxon_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows)

> Collect all 'good' MAGs into a shared folder

##### Bin annotation

Now we should have a collection of MAGs that we can further analyze. The first step is to predict genes again, as right now we only have raw genomic sequences. We will use a different tool this time, one of my all-time-favorites : `prokka`.

This tool does gene prediction as well as some quiet good annotations, and is actually quiet easy to run!

> Use `prokka` loaded with `module`

> Predict genes and annotate your MAG!

> Optional : use `prokka` on some of the bins that did not pass the previous quality checks!

`prokka` produces a number of output files that all kind of represent similar things. Mostly variants of FASTA-files, one with the genome again, one with the predicted proteins, one with the genes of the predicted proteins. Also it renames all the sequence with nicer IDs! Additionally a very useful file generated is a GFF-file, which gives more information about the annotations then just the names you can see in the FASTA-files.


You will explore this tomorrow some more with Maliheh.

We now know more about the genes your MAG contains, however we do not really know who we have?! `checkm` might have given us an indication but it is only approximative.

Taxonomical classification for full genomes is not always easy for MAGs, often the 16S gene is missing as it assembles badly, and which other genes to use to for taxonomy is not always evident. One option is to take many genes and to make a tree.

One tool implementing such a thing is `phylophlan`, it has a reference tree based on a concatenation of reference genes, tries to find homologues of these genes in your bin, and then to fit  them into the reference tree.

> Install `phylophlan`

> Use it to identify your MAG(s)

##### A bit of phylogenomics

If we have time!
