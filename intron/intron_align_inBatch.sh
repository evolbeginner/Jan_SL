#! /bin/bash


########################################################################
intron_align=$SL/scripts/intron/intron_align.rb
intron_sliding_cutoff=5


########################################################################
while [ $# -gt 0 ]; do
	case $1 in
		--indir)
			indir=$2
			shift
			;;
		--outdir)
			outdir=$2
			shift
			;;
		--gff)
			gff=$2
			shift
			;;
		--intron_sliding|--intron_sliding_cutoff)
			intron_sliding_cutoff=$2
			shift
			;;
		*)
			echo "Error param!"
			echo "$1"
			echo "Exiting ......"
			exit
			;;
	esac
	shift
done


mkdir -p $outdir


########################################################################
for i in $indir/*pep.aln; do
	basename=`basename $i`;
	a=`grep -o -P '(.+)(?=.pep.aln)' <<< $basename`;
	b=`sed 's/[-]/|/' <<< $a`;
	grep -P $b $gff > $outdir/test.gff;
	ruby2.1 $intron_align -i $i --gff $outdir/test.gff --intron_sliding $intron_sliding_cutoff;
	#ruby2.1 $intron_align -i test_SL/symbB.v1.2.036605-symbB.v1.2.021957.pep.aln --gff 1.gff
done > $outdir/result


