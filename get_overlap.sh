#! /bin/bash

gff2bed_path="/mnt/bay3/sswang/software/NGS/basic_processing/bedops/"
gff2bed=$gff2bed_path/gff2bed
export PATH=$gff2bed_path:$PATH
export PATH="/home/sswang/huche/help/Jan/SL/scripts":$PATH

blast2bed="blast2bed.pl"

e_value=10000
#############################################################################
while [ $# -gt 0 ]; do
	case $1 in
		--outdir)
			outdir=$2;	shift;
			;;
		--upstream)
			upstream=$2;	shift;
			;;
		--downstream)
			downstream=$2;	shift;
			;;
		--feature)
			feature=$2;	shift;
			;;
		--gff)
			gff1=$2;	shift;
			;;
		--blast_result)
			blast_result=$2;	shift;
			;;
		--e_value)
			e_value=$2;	shift;
			;;
	esac
	shift
done

if [ -z $outdir ]; then
	echo "outdir should be specified by '--outdir'";
	exit;
fi
[ ! -e $outdir ] && mkdir $outdir

#############################################################################
bed1="$outdir/from_gff.bed"
bed2="$outdir/from_blast.bed"
gff1_basename=`basename $gff1`
new_gff1="$outdir/$gff1_basename"

create_new_gff.rb --gff $gff1 --upstream $upstream --feature $feature > $new_gff1
$gff2bed < $new_gff1 > $bed1

blast2bed.pl --blastm8 $blast_result --out_prefix from_blast --e_value $e_value
mv from_blast.bed $outdir

bedtools intersect -a $bed1 -b $bed2 -s

