# Metagenomics Course, Fall 2019


### Binning the contigs

Nice, we have genomic-data with some patterns in it. Now we want to extract genomes from it, so we will use a binning approach with our assembly. To extract genomic bins, or Metagenome Assembled Genomes (MAGs), or Genomes From Metagenomes (GFMs), we will use the `MetaBAT` tool.

This tool will need the sorted BAM-files we made to run.

> For you genomic-assembly of choice. Cfind the sorted indexed BAM-file you made assembly in a folder. (Remember the `bamscript=` option of `bbmap`), now is the time!
> Using `runMetaBat.sh`, loaded with `module load MetaBaT`. Have a look at the FASTA-files generated.
> Optional : run the binning on an other assembly, or with more mappings!

#### Evaluating bins

Now we each should have all gotten a number of bins. Meaning that the assembly has been split into bags of contigs that somehow look/vary similarly. This process is however far from perfect, some bins will be-merged MAGs, others will contain viral data, or Eukaryotic chromosome, or simply trash.

We need to evaluate the bins for that. A very simple first filter is to remove bins that are too small or too large.

> compute rough size estimates of your bins (hint, the file size correlate well with the genome size)

Also, the `R`-script previously used colours by bin (yey!), so let's make our plot prettier:

> Change the script so it finds your bins, and run it again, it should make the figure prettier!
> optional : compute other stats on these bins, let your imagination run wild (or maybe use the script I provided earlier)

Now we should all have a bin we hope is a genome! To check this we will use a tool called `checkm`, fundamentally it uses single copy marker genes to check if a bin/MAG is complete, and if it is contaminated.

> Load `checkm`. You will need once again UPPMAXs `module` system.

Once it seams like it is running we want to run the `lineage_wf`, the documentation of `checkm` is a bit more confusing than other tools... So I advice you to look at [this specific page](https://github.com/Ecogenomics/CheckM/wiki/Quick-Start) and don't get lost in the others!

> Using the `lineage_wf` of `checkm`
> Compute completeness/contamination for your bins. Did you get any bonus information?
> Optional : [taxon_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows)
> Collect all your 'good' MAGs into the shared folder `/proj/g2019027/2019_MG_course/MAG_collection/`, make sure there is no other MAG there with that name.

#### Functional annotation

Now we should have a collection of MAGs that we can further analyze. The first step is to predict genes as right now we only have raw genomic sequences. We will use one of my all-time-favorites : `prokka`.

This tool does gene prediction as well as some decent and usefull annotations, and is actually quiet easy to run!

> Use `prokka` loaded with `module`
> Predict genes and annotate your MAGs!
> Optional : use `prokka` on some of the bins that did not pass the previous quality checks!

`prokka` produces a number of output files that all kind of represent similar things. Mostly variants of FASTA-files, one with the genome again, one with the predicted proteins, one with the genes of the predicted proteins. Also it renames all the sequence with nicer IDs! Additionally a very useful file generated is a GFF-file, which gives more information about the annotations then just the names you can see in the FASTA-files.

The annotations of `prokka` are good but not very complete for environmental bacteria. Let's run an other tool I like a lot, eggNOGmapper. This is a bit heavier in computation and it is not on  `UPPMAX`, so you will have to [install](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2) it.

> Install eggNOG-mapper
> run it on at least one MAG (don't be greedy, it is not fast!)
> OPTIONAL : undestand the output of it ...

You will explore this tomorrow some more with Maliheh.

#### Taxonomic annotation

We now know more about the genes your MAG contains, however we do not really know who we have?! `checkm` might have given us an indication but it is only approximative.

Taxonomic classification for full genomes is not always easy for MAGs, often the 16S gene is missing as it assembles badly, and which other genes to use to for taxonomy is not always evident. Typically marker genes, min-hashes or k-mer databases are used as reference. It is often problematic for environmental data as the databases are not biased into our direction! We will use a min-hash database I compiled specifically for this (based on other tools and the full-datasets) using a tool called [sourmash](https://sourmash.readthedocs.io/en/latest/).

> Install sourmash!
> Use the lca function of sourmash with the database found here /proj/g2019027/2019_MG_course/dbs/our_taxonomical_database.json

This database was made by MAGs annotated with an other [tool](https://github.com/Ecogenomics/GTDBTk), `gtdbtk`, it uses marker genes and loads of data. It is a bit heavy, and tricky to run/install but much more sensitive.

> Optional: install and run gtdbtk, you can use the database here `/proj/g2019027/2019_MG_course/dbs/gtdbtk`
