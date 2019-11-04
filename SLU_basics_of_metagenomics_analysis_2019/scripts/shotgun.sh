9-915, meetandgreet
1.15hrs of Chris greening on applied metagenomics



module load bioinfo-tools
module load FastQC
module load trimmomatic
module load megahit
module load seqtk

course_folder=~/people/0023_anoxicencyclo/course/
### prepping data

nb_reads=3000000
raws=0000_raws/0100_reads/genomic2/Loc

cd ~/people/0023_anoxicencyclo/

mkdir -p $course_folder/raws

mkdir -p $course_folder/used_genomes
mkdir -p $course_folder/all_bins


for f in `cat $course_folder/keepers4courseWithMonos`; do cp 4500_assembly_analysis/genomics/genomes/${f}.fna $course_folder/used_genomes/; done

cat $course_folder/used_genomes/*.fna $course_folder/all_genomes.fna

bbmap.sh ref=$course_folder/all_genomes.fna path=$course_folder/mock/
for sample in `ls $raws* | cut -f4 -d/  | cut -f1 -d_ | uniq`;
do
  echo $sample;
  mkdir -p $course_folder/mock/$sample
  if [ ! -f $course_folder/mock/$sample/${sample}.sam ];
  then
    echo "mapping them"
    bbmap.sh mappedonly=t showprogress=100000 ref=$course_folder/all_genomes.fna path=$course_folder/mock/ in=${raws%%Loc}/${sample}_L001_R1.fastq.gz in2=${raws%%Loc}/${sample}_L001_R2.fastq.gz out=$course_folder/mock/$sample/${sample}.sam bamscript=$course_folder/maps/$sample/${sample}.sh 2> $course_folder/mock/$sample/${sample}.log
  fi
  if [ ! -f $course_folder/mock/$sample/${sample}.sub.fwd.fastq ];
  then
    echo "extracting them"
    cat $course_folder/mock/$sample/${sample}.sam | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$11"\n+\n"$12}' > $course_folder/mock/$sample/${sample}.fwd.fastq
    cat $course_folder/mock/$sample/${sample}.sam | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$11"\n+\n"$12}' > $course_folder/mock/$sample/${sample}.rev.fastq
    echo "subsetting them"
    seqtk sample -s100 $course_folder/mock/$sample/${sample}.fwd.fastq $nb_reads > $course_folder/mock/$sample/${sample}.sub.fwd.fastq
    seqtk sample -s100 $course_folder/mock/$sample/${sample}.rev.fastq $nb_reads > $course_folder/mock/$sample/${sample}.sub.rev.fastq
  fi
done


mkdir -p $course_folder/cleans/

for sample in `ls $raws* | cut -f4 -d/  | cut -f1 -d_ | uniq`;
do
  echo $sample;
  mkdir -p $course_folder/mock/$sample
  if [ ! -f $course_folder/cleans/${sample}_1U ];
  then
    echo "QC them"
    java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -phred33 -threads 20  $course_folder/mock/$sample/${sample}.sub.fwd.fastq $course_folder/mock/$sample/${sample}.sub.rev.fastq  -baseout $course_folder/cleans/${sample} ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  fi
done


mkdir -p $course_folder/assemblies
for sample in `ls $raws* | cut -f4 -d/  | cut -f1 -d_ | uniq`;
do
  if [ ! -f $course_folder/assemblies/$sample/final.contigs.fa ]
  then
    echo assembling $sample;
    megahit --continue -t 20 --min-contig-len 1500 -1 $course_folder/cleans/${sample}_1P  -2 $course_folder/cleans/${sample}_2P -o $course_folder/assemblies/$sample
  fi
done


for sample in `ls $raws* | cut -f4 -d/  | cut -f1 -d_ | uniq`;
do
  if [[ -f $course_folder/assemblies/$sample/final.contigs.fa && ! -f $course_folder/maps/$sample/${sample}_sorted.bam ]]
  then
    echo maapping back $sample;
    mkdir -p $course_folder/maps/$sample
    bbmap.sh ref=$course_folder/assemblies/${sample}/final.contigs.fa path=$course_folder/maps/$sample 2> /dev/null
    bbmap.sh ref=$course_folder/assemblies/${sample}/final.contigs.fa path=$course_folder/maps/$sample in=$course_folder/cleans/${sample}_1P in2=$course_folder/cleans/${sample}_2P out=$course_folder/maps/$sample/${sample}.sam bamscript=$course_folder/maps/$sample/${sample}.sh 2> $course_folder/maps/$sample/${sample}.log
    $course_folder/maps/$sample/${sample}.sh
  fi
done

for sample in `ls $raws* | cut -f4 -d/  | cut -f1 -d_ | uniq`;
do
  if [[ -f $course_folder/maps/$sample/${sample}_sorted.bam && ! -f $course_folder/bins/$sample/final.contigs.fa.depth.txt ]]
  then
    echo binning $sample;
    mkdir -p $course_folder/bins/$sample
    runMetaBat.sh $course_folder/assemblies/${sample}/final.contigs.fa $course_folder/maps/$sample/${sample}_sorted.bam 2> /dev/null > /dev/null
    mv final.contigs* $course_folder/bins/$sample
  fi
done

for f in `ls $course_folder/bins/Loc0*/final.contigs.fa.metabat-bins/*.fa`;
do
  cp $f $course_folder/all_bins/`echo $f | rev | cut -f3 -d/ | rev`_`basename $f`;
done

cd $course_folder/all_bins/
for f in `ls *.fa`
do
  if [[ ! -f ${f%%.fa}/${f%%.fa}.gff ]]
  then
    echo prokka $f
    prokka --prefix ${f%%.fa} --cpus 20  --outdir ${f%%.fa} $f 2> /dev/null ;
  fi
done

for f in `ls *.fa`; 
do
if [[ ! -f ${f%%.fa}/${f%%.fa}.emapper.annotations ]]
  then
    echo emap  $f;
    emapper.py --cpu 20 -i ${f%%.fa}/${f%%.fa}.faa -o ${f%%.fa}/${f%%.fa} -m diamond;
  fi
done


sourmash compute -k 31 --scaled 10000 -p 20 *.fa
sourmash lca classify --db /home/moritz/people/0023_anoxicencyclo/1500_coasses/Loclat/assemblies/megahit/binning/metabat/sourmash/local.lca.json --query  *.sig > class.tax
checkm lineage_wf -t 20 -x fa . checkm > checkm.txt
