#!/bin/bash
####Modified by marcguery on Fri 16 Apr 17:33:43 CEST 2021########


################################PARAMETERS################################

bam=$1
kmer=$2

n=$LSB_JOBINDEX
if [ -z "$n" ] ; then
  n=5;
fi

################################################################


################################CONFIG################################

##Data files
MASTERLIST=MasterList.txt #NA
LITTLEPF3D7=/nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Feb2015/little.Pf3D7.Feb2015.txt #NA
LITTLEPF3D7V4=/nfs/pathdata2/Plasmodium/falciparum/3D7/Reference/Little.Pf3D7V4.txt #NA
BAMFOLDER=/nfs/team112_data00/pf.bwa/bam/
VARCONFIG=/lustre/scratch118/infgen/team133/tdo/Pfalciparum/VAR/Assembly.Version3/config.txt #NA
CAFOLDER=CA/
GAPCLOSEFOLDER=10-gapclose/
SCFFASTA=genome.scf.fasta

bam=$( head -n $n $MASTERLIST | tail -n 1 | awk '{ print $1}');
name=$( head -n $n $MASTERLIST | tail -n 1 |  awk '{ print $2}');


##Tools
#sh
ANNOTVARGENES=~/Bin/vargene.AnnotationVARgenes.assembly.V0.1.sh
MASPOSTIMP=~/Bin/vargene.AssemblyPipeline.masurca.PostImprovement.NoReadextraction.sh
ASSSTATS=~tdo/Bin/assstats.sh #NA
#Perl
MERGEOSAM=~/Bin/vargene.mergeOverlappingSamtoolsRegions.pl
SAMTOFASTQ=~/Bin/sam2fastq.2files.pl #NA
JOINSOLEX=~tdo/Bin/join.solexaBack.actual.sga.pl #NA

#Common
SAMTOOLS1_2=samtools1.2
SAMTOOLS=samtools
FASTAQ=fastaq
SGA=sga
MASURCA=/nfs/pathogen003/tdo/bin/MaSuRCA-2.3.2/bin/masurca
VELVETH=~/bin/velveth

#Others
BSUB=bsub #NA
VELVETHTDO=~tdo/bin/velveth.tdo #NA
VELVETG=~/bin/velvetg #NA
VELVETGTDO=~tdo/bin/velvetg.tdo #NA

################################################################


################################SCRIPT################################

####Version 4: read correction####
#set -e

if [ -z "$bam" ] ; then
	echo "Please don't use that script!"
	exit 1;
fi

if [ -z "$kmer" ] ; then
	kmer=71
fi

if [[ $bam == *cram* ]] ; then
   $SAMTOOLS1_2 view -F 2304 -b $bam > $name.bam 
else
	ln -s $bam $name.bam
fi

$SAMTOOLS index $name.bam
#ln -s $bam.bai $name.bam.bai

echo "$n -- $bam\n";
#name=$(echo $bam | sed 's/.bam//g')

mkdir -p Assembly.$name
cd Assembly.$name

finalDIR=$(pwd)
deleteIT=0
lfs setstripe -c -1 .

###
###test 3D7/
#$SAMTOOLS view -F 4 ../$bam  MAL4:171451-180412  > MappedVAR.$bam.sam

bam=$name.bam

# General test
insertSize=350


# ADAPT ofr _v3 self mappeind
region1=$(egrep "stevor|PIR|PfEMP1" $LITTLEPF3D7 | \
			perl -nle '@ar=split(/\t/); @bb=split(/\.\./,$ar[4]);($a) = $bb[0] =~ /\(*(\d+)/; ($b)= $bb[(scalar(@bb)-1)] =~ /(\d+)\)*/; $chr=substr($ar[0],3,5); print "Pf3$chr:".($a-500)."-".($b+500)'  | \
			awk '{ s=s" "$1} END { print s }' | \
			perl $MERGEOSAM )

region=$(egrep "stevor|PIR|PfEMP1" $LITTLEPF3D7 | \
			perl -nle '@ar=split(/\t/); @bb=split(/\.\./,$ar[4]);($a) = $bb[0] =~ /\(*(\d+)/; ($b)= $bb[(scalar(@bb)-1)] =~ /(\d+)\)*/; $chr=substr($ar[0],3,5); print "Pf3$chr\_v3:".($a-500)."-".($b+500)'  | \
			awk '{ s=s" "$1} END { print s }' | \
			perl $MERGEOSAM )

region2=$(egrep "stevor|PIR|PfEMP1|rif" $LITTLEPF3D7V4 | \
			perl -nle '@ar=split(/\t/); @bb=split(/\.\./,$ar[4]);($a) = $bb[0] =~ /\(*(\d+)/; ($b)= $bb[(scalar(@bb)-1)] =~ /(\d+)\)*/; $chr=substr($ar[0],3,5); print "$ar[5]:".($a-500)."-".($b+500)'  | \
			awk '{ s=s" "$1} END { print s }' | \
			perl $MERGEOSAM i)

#region=$(egrep "stevor|PIR|PfEMP1" $LITTLEPF3D7 | \
			# perl -nle '@ar=split(/\t/); @bb=split(/\.\./,$ar[4]);($a) = $bb[0] =~ /\(*(\d+)/; ($b)= $bb[(scalar(@bb)-1)] =~ /(\d+)\)*/; $chr=substr($ar[0],3,5); print "Pf3$chr\_v3:".($a-500)."-".($b+500)'  | \
			# awk '{ s=s" "$1} END { print s }' | \
			# perl $MERGEOSAM )

#echo $region

if [ "0" == "1" ] ; then
	#rm NonMapped.$bam.sam MappedVAR.$bam.sam
	for x in `ls $BAMFOLDER/$name*bam`; do 
		echo " getting reads from $x";
		$SAMTOOLS view -F 4 $x $region $region1 $region2 >> MappedVAR.$bam.sam
		$SAMTOOLS view -f 4 $x >> NonMapped.$bam.sam

		insertSize=$($SAMTOOLS view -f 2 $x  MAL14:200000-1000000 | \
						cut -f 9 | \
						sort -n | \
						awk 'BEGIN{c=0;} {if($1>0) {a[c++]=$1; }} END{ median=a[int(c/2)]; print median}')

		#insertSize=$($SAMTOOLS view -f 2 $BAMFOLDER/$x  MAL14:200000-1000000 | \
						# cut -f 9 | \
						# sort -n | \
						# awk 'BEGIN{c=0;} {if($1>0) {a[c++]=$1; }} END{ median=a[int(c/2)]; print median}')
	done

	RL=$($SAMTOOLS view  $x  MAL14 | head -n 1 | cut -f 10 | perl -nle 'print length($_)');

else
	# touch NonMapped.$bam.sam
	if [ "1" == "1" ] ; then
		echo "using local bam file "
		$SAMTOOLS view  ../$bam $region $region1 $region2 > MappedVAR.$bam.sam
		### get all the relevant positions

		### changed June 2019
		# insertSize set to 500 
		# RL adapter
		$SAMTOOLS view -f 4 ../$bam > NonMapped.$bam.sam
		insertSize=$($SAMTOOLS view -f 2 ../$bam  Pf3D7_04_v3:200000-1000000 | \
						cut -f 9 | \
						sort -n | \
						awk 'BEGIN{c=0;} {if($1>0) {a[c++]=$1; }} END{ median=a[int(c/2)]; print median}')
		
		RL=$($SAMTOOLS view  ../$bam  | \
				head -n 1 | \
				cut -f 10 | \
				perl -nle 'print length($_)');
	
	else 
		#longer smal test  MAL4:169451-182412 MAL4:550144-624339  MAL4:936588-994282 MAL4:1156045-1185848
		$SAMTOOLS view -F 4 ../$bam  MAL4:31153-71936  > MappedVAR.$bam.sam
		insertSize=$($SAMTOOLS view -f 2 ../$bam  MAL14:200000-1000000 | \
						cut -f 9 | \
						sort -n | \
						awk 'BEGIN{c=0;} {if($1>0) {a[c++]=$1; }} END{ median=a[int(c/2)]; print median}')
		
		RL=$($SAMTOOLS view  ../$bam  MAL14 | \
		head -n 1 | \
		cut -f 10 | \
		perl -nle 'print length($_)');
	fi
fi

if [ -z $RL ] ; then 
	RL=100
fi 

if [ -z $insertSize ] ; then 
	insertSize=500
fi 

if [ "$RL" -gt 80 ] ; then
	kmer=71
else
	kmer=61
fi

echo "insert size is $insertSize. The readlength is $RL. "
echo "kmer is $kmer"

cat MappedVAR.$bam.sam NonMapped.$bam.sam  | perl $SAMTOFASTQ rawread_tmp

$FASTAQ enumerate_names rawread_tmp_1.fastq rawreads_1.fastq
$FASTAQ enumerate_names rawread_tmp_2.fastq rawreads_2.fastq

mv rawread_tmpSE.fastq rawreadsSE.fastq
rm rawread_tmp*fastq

#cat MappedVAR.$bam.sam   | perl $SAMTOFASTQ $bam
#rm MappedVAR.$bam.sam NonMapped.$bam.sam
fastq=$bam

### read correction
#$SGA preprocess -m 51 --permute-ambiguous -f 3 -q 3 -p 1 rawreads_1.fastq rawreads_2.fastq > tmp.fastq

#$SGA index -t 2 --disk=2000000  tmp.fastq
#$SGA correct -k 51 -x 3 -t 2  -o read.corrected.fastq tmp.fastq; 
#perl $JOINSOLEX read.corrected.fastq $fastq

### kick off the masurca
cp $VARCONFIG config.txt
$MASURCA config.txt
./assemble.sh

cp $CAFOLDER/$GAPCLOSEFOLDER/$SCFFASTA SC.$name.fasta

#mv $GAPCLOSEFOLDER .
# uncomment the next 1,2,4 lines
rm -rf $CAFOLDER
rm assemble.sh config.txt environment.sh gapClose.err \
	genome.uid global_arrival_rate.txt guillaumeKUnitigsAtLeast32bases_all.fasta \
	guillaumeKUnitigsAtLeast32bases_all.jump.fasta k_u_hash_0 *.sam meanAndStdevByPrefix.pe.txt \
	pA.linking.frg pA.renamed.fastq pe.cor.fa pe.cor.log pe_data.tmp pe.linking.fa quorum.err \
	quorum_mer_db.jf rawreadsSE.fastq runCA1.out runCA2.out runCA3.out super1.err superReadSequences_shr.frg \
	tigStore.err unitig_cov.txt unitig_layout.txt

rm -rf work1
#cp $GAPCLOSEFOLDER/$SCFFASTA SC.$name.fasta
#rm -rf $GAPCLOSEFOLDER
#$ANNOTVARGENES SC.$name.fasta empty.bam $name

cd ..
$BSUB -R "select[type==X86_64 && mem > 9000] rusage[mem=9000,tmp=4000]" -M9000 \
	-o Output/out.Improvement.$name.o -e Output/out.Improvement.$name.e -n 4 \
	-R "span[hosts=1]" \
	bash $MASPOSTIMP $name

###now improve the assembly by default
### do de novo assembly using predifined k-mer of 71
#for i in  61 71 81  ; do 
#	NAME=k.denovo.$fastq.$i
#	$VELVETH $NAME $i -shortPaired  -fastq -separate  "$fastq"_1.fastq "$fastq"_2.fastq 
#	$VELVETHTDO $NAME $i -shortPaired  -fastq -separate  rawreads_1.fastq rawreads_2.fastq 
#	echo "$VELVETG $NAME -exp_cov auto -cov_cutoff 5 -very_clean yes -ins_length $insertSize  -min_contig_lgth 500 "
#	$VELVETGTDO $NAME -exp_cov auto -cov_cutoff 5 -very_clean yes -ins_length $insertSize  -min_contig_lgth 300  
#done

#if [ ! -f "k.denovo.$fastq.$kmer/contigs.fa" ] ; then
#	echo "Error in velvet, maybe read not long enough"
#	exit 125;
#fi
#$ASSSTATS k.denovo.$fastq.$kmer/contigs.fa  > stats.VelvetAssembly.txt

################################################################
