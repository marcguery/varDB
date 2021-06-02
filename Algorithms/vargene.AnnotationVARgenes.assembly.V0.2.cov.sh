#!/bin/bash
####Modified by marcguery on Fri 16 Apr 17:33:43 CEST 2021########


################################PARAMETERS################################

genome=$1;
readroot=$2
name=$3

AUGUSTUS_CONFIG_PATH=/nfs/pathogen003/tdo/Tools/Augustus/config; export AUGUSTUS_CONFIG_PATH
SMALT_PARAMETER=" -n 1 -y 0.7 " ; export SMALT_PARAMETER

################################################################


################################CONFIG################################

##Data files
tmp=$$;
olddir=$PWD
dir="/tmp/tdo/dir.tdo.$tmp/"
#dir="local"
REF3D7REFAA=/nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Oct2011/Pf3D7_v3.aa.fasta
UPSSEQ=/nfs/pathdata2/Plasmodium/falciparum/3D7/Reference/ups.sequences.txt
VARSEQ=/nfs/pathdata2/Plasmodium/falciparum/3D7/Reference/PfV2.internalVAR.fasta

##Tools

#Commons
AUGUSTUS=augustus
SAMTOOLS=samtools
FORMATDB=formatdb
BLASTALL=blastall
JAVA=/software/java/bin/java

#sh
ANNOTATEAUGUSTUS=~/Bin/annotation.wrapperAnnotateAugustus.sh #NA
LITTLEEMBL=~tdo/Bin/little.allEMBL2xx.sh #NA
LITTLESMALT=~/Bin/little.smalt.bam.sh #NA
#Python

#Perl
GETANNOTFASTA=getAnnoFasta.pl #NA
SEPSEQ=SeperateSequences.pl #NA
MAINRATT=~/Bin/ratt/main.ratt.pl #NA
PSUUNION=psu_union.pl #NA
EMBLMOREDESC=/nfs/users/nfs_t/tdo/Bin/EMBL2aa_moreDescription.pl #NA
EMBLCDSREICH=/nfs/users/nfs_t/tdo/Bin/EMBL2CDS.Reichenowi.previousID.pl #NA
FASTAEXTRACTOR=fasta_extractor.pl #NA
FASTA2SINGLELINE=fasta2singleLine.pl #NA
DSIDCYS=~/Bin/DSID_cys.pl #NA
COUNTCONTIGSIZE=/nfs/users/nfs_t/tdo/Bin//Countcountigsize.pl #NA
RAWREADNORM=~tdo/Bin/vargene.rawReadCountNormalization.pl
FASTA2GTF=~/Bin/fasta2gtf.pl #NA

#Others
MARKDUPL=$ICORN2_HOME/MarkDuplicates.jar

################################################################


################################SCRIPT################################

echo "Old dir $olddir   New Dir $dir" 

mkdir -p $dir
cp $genome $dir
cd $dir
cp -s $olddir/$readroot*q .
### Getting some annotation. So far based on augustus
$AUGUSTUS --species=3D7_genefamilies $genome  > out.augustus.txt
### rename the output to the correct ID
cat out.augustus.txt | \
    perl -e 'my $n=shift; while(<STDIN>){s/g(\d+)/$n.g$1/g; print }' $name > out2.augustus.txt

#$GETANNOTFASTA out2.augustus.txt; mv out.augustus.aa Genes.$name.aa.fasta
### to embl
$ANNOTATEAUGUSTUS $genome $REF3D7REFAA out2.augustus.txt

mkdir Seq
cd Seq
$SEPSEQ ../$genome  &> /dev/null
cd ..
mkdir EMBL
for x in ` grep '>' $genome | sed 's/>//g' ` ; do 
    perl $MAINRATT  doEMBL EMBL/$x Augustos/$x.embl Seq/$x
done &> Annotation.txt
#$PSUUNION EMBL/*embl > $name.union.embl
rm -rf Seq/

### get a union file and all the genoems with VAR annotation
list=$(ls -S EMBL/* | awk '{s=s" "$1} END{print s}')
$PSUUNION $list > Union.$name.embl
rm -rf EMBL/;
rm -rf Augustos*/
#perl $EMBLMOREDESC Union.$name.embl > $name.all.aa.fasta

touch All.$name.aa.fasta All.$name.na.fasta

$LITTLEEMBL $name
mv All.$name.aa.fasta $name.all.aa.fasta
mv All.$name.na.fasta $name.all.nt.fasta
#perl $EMBLCDSREICH Union.$name.embl $name.all.nt.fasta 0
grep -i PfEMP1 $name.all.aa.fasta | awk '{ print  $1 }' | sed 's/>//g' >  id.withPfEMP1.annotation.txt

### get just the genes annotated with Pfemp1
$FASTAEXTRACTOR -i $name.all.aa.fasta -f id.withPfEMP1.annotation.txt -s $name.VAR.aa.fasta &> /dev/null
$FASTAEXTRACTOR -i $name.all.nt.fasta -f id.withPfEMP1.annotation.txt -s $name.VAR.nt.fasta &> /dev/null


### get sis data
$FASTA2SINGLELINE $name.VAR.aa.fasta | \
    perl -nle 'if (/>/){print } elsif(/(DI\wD[I|V]\wR.{10,160}DY\wPQ\w{2}R)/) { print $1} else { print "NA"};' | \
    grep -B 1 DI | grep -v "^\-\-" > dbl.$name.fasta 
perl $DSIDCYS dbl.$name.fasta cys.$name.txt  &> /dev/null
cut -f 3 cys.$name.txt | \
    sort | \
    grep -v cys | \
    uniq -c  | \
    perl -nle 's/\s+/\t/g; print '> cys.$name.summary.txt

### differential expression
# it will be jsut done on the ..nt...
# map reads with smalt

genome1k=$name.VAR.nt.0.5kb.fasta
genome=$name.VAR.nt.fasta
$COUNTCONTIGSIZE $genome 500 $genome1k

$LITTLESMALT $genome1k 20 8 $readroot\_1.fastq $readroot\_2.fastq Res.Mapped 1000
$JAVA -XX:MaxPermSize=512m -Xmx2000m -XX:ParallelGCThreads=1 -XX:-UseParallelGC -XX:+UseSerialGC -jar \
    $MARKDUPL VALIDATION_STRINGENCY=SILENT M=/dev/null ASSUME_SORTED=TRUE \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=Res.Mapped.bam O=Res.Mapped.markdup.bam

## We will count the reads per var gene longer than 1kb
$SAMTOOLS view -F 1028 Res.Mapped.bam | \
    cut -f 3 | \
    uniq -c | \
    perl -nle '/^\s*(\d+)\s+(\S+)/; print "$2\t$1" ' > ReadCount.0.5k.txt
amount=$(awk '{ s+=$2} END { print s}' ReadCount.0.5k.txt )
$SAMTOOLS faidx $genome1k
cat ReadCount.0.5k.txt | \
perl $RAWREADNORM $genome1k.fai  $amount > ReadCount.$name.txt
#$FASTA2GTF $genome > $genome.gtf
#awk '$2=="AUGUSTUS"' out2.augustus.txt | \
    # egrep "transcript|exon" | \
    # sed 's/CDS/exon/g' | \
    # perl -e 'while(<STDIN>){ if (/exon/){@ar=split(/\t/); $CDS{$ar[0]}++; chomp($_); print "$_ exon_number \"$CDS{$ar[0]}\"\n"  }}' > ForCufflinks.gtf

#bam=Res.Mapped.bam
#cufflinks  --no-update-check -b $genome -q  -G $genome.gtf -o cuff.coverage $bam


### blast of the UPS sequences / interal regions.
$FORMATDB -p F -i $genome
$BLASTALL -m 8 -F F -p blastn -e 1e-2 -d $genome -i $UPSSEQ -o comp.ups.$name.blast
$BLASTALL -m 8 -p blastn -e 1e-80 -d $genome -i $VARSEQ -o comp.internal.$name.blast


cp Read* cys*  *VAR*a *.all.*.fasta  Union*embl  $olddir

rm out*txt
rm Union*embl.??.fasta *fai
rm $genome* formatdb.log  id.withPfEMP1.annotation.txt
rm *fastq 
cp Res.Mapped.bam $olddir
zip -r $olddir/Annotation.zip *

echo "Old dir $olddir   New Dir $dir" 

cd $olddir
rm -rf $dir/
echo "done."
################################################################
