





for f in `find . -name "*_R1*"`;
do
  id=`dirname $f | cut -f2-3 -d"/" | tr '/' '_'`
  out_file=`dirname $f | cut -f2-3 -d"/"`/$id
  if [[ -e $out_file.fastq ]] ; then
    i=0
    while [[ -e $out_file-$i.fastq ]] ; do
        let i++
    done
    out_file=$out_file-$i
  fi
  std_out=$out_file.log
  out_file=$out_file.fastq
  usearch10 -fastq_mergepairs $f -fastqout $out_file --relabel "$id:" 2> $std_out > /dev/null
  echo $out_file
  cat $out_file >> all_reads.fastq
done

### Question can you say something about the merging rate?! improve the rate maybe?!


usearch10 -fastq_filter all_reads.fastq -fastq_maxee 1 -fastaout all_reads.clean.fasta
usearch10 -fastx_uniques all_reads.clean.fasta -relabel Uniques_ -sizeout --minuniquesize 2 -fastaout all_reads.uniques.fasta
usearch10 -cluster_otus all_reads.uniques.fasta  -relabel OTUs_ -otus  all_reads.OTUs.fasta
usearch10 -usearch_global all_reads.fastq -db all_reads.OTUs.fasta -strand plus -id 0.97 -otutabout OTU_table.txt

### How do hit rates look for specific libs?!

usearch10 -sintax all_reads.OTUs.fasta -db ~/Data/tax_dbs/rdp_16s_v16_sp.udb -tabbedout taxonomy.tax -strand both -sintax_cutoff 0.8

usearch10 -unoise3 all_reads.uniques.fasta -zotus all_reads.zOTUs.fasta
usearch10 -usearch_global all_reads.fastq -db all_reads.zOTUs.fasta -strand plus -id 0.97 -otutabout zOTU_table.txt
usearch10 -sintax all_reads.zOTUs.fasta -db ~/Data/tax_dbs/rdp_16s_v16_sp.udb -tabbedout ztaxonomy.tax -strand both -sintax_cutoff 0.8
