#! /bin/bash


####################################################################
dirname=`dirname $0`
cd $dirname;
wd=`pwd`
cd -

SL=~/huche/help/Jan/SL
blast_exe=~/tools/self_bao_cun/blast_from_fasta/blast_exe.sh
find_relict_SLs=$SL/scripts/find_relict_SLs.rb 
getUpstreamGff=$wd/getUpstreamGff.rb
blast2gff=$wd/blast2gff.rb
extract_slList=$wd/extract_slList.rb
filter_BLAST_sl=$wd/filter_BLAST_sl.rb
filter_blast_result=$wd/../filter_blast_result.rb
cpu=2
evalue=1e-10
#query_start_max=10


####################################################################
while [ $# -gt 0 ]; do
	case $1 in
		-i|--in)
			infile=$2
			shift
			;;
		--db_fasta)
			db_fastas=(${db_fastas[@]} $2)
			shift
			;;
		--sl)
			sl_file=$2
			shift
			;;
		--outdir)
			outdir=$2
			shift
			;;
		--evalue)
			evalue=$2
			shift
			;;
		--evalue_sl|--evalue_SL)
			evalue_sl=$2
			shift
			;;
		--query_start_max)
			query_start_max=$2
			shift
			;;
		--cpu|--CPU)
			cpu=$2
			shift
			;;
		--force)
			is_force=true
			;;
		*)
			echo "Unknown arguments! Exiting ......"
			exit
	esac
	shift
done


[ -z $infile ] && echo "infile not give! Exiting ......" && exit
[ -z $db_fastas ] && echo "db_fasta not give! Exiting ......" && exit
[ -z $sl_file ] && echo "sl_file not give! Exiting ......" && exit
[ -z $outdir ] && echo "outdir not give! Exiting ......" && exit
[ -z $evalue ] && echo "evalue not give! Exiting ......" && exit
[ -z $evalue_sl ] && echo "evalue_sl not give! Exiting ......" && exit
[ -z $query_start_max ] && echo "query_start_max not give! Exiting ......" && exit


if [ -d $outdir ]; then
	if [ ! -z $is_force ]; then
		rm -rf $outdir
	fi
fi
mkdir -p $outdir

mkdir -p $outdir/fasta
db_fasta=$outdir/fasta/input.fasta
cat ${db_fastas[@]} > $outdir/fasta/input.fasta


####################################################################
echo "Blasting against EST ......"
$blast_exe --input $infile --db $db_fasta --evalue $evalue --type nucl --outdir $outdir/BLAST_cds --outfmt 6 --CPU $cpu

echo "Finding SL relict ......"
ruby2.1 $find_relict_SLs -i $db_fasta --outdir $outdir/BLAST_sl/ori --SL_seq $sl_file --quiet --remove_null --blast_outfmt 1 --bl2seq --strand both


ruby2.1 $filter_blast_result --indir $outdir/BLAST_sl/ori --outdir $outdir/BLAST_sl/selected/e10 --action cp --blast_outfmt 1 --aligned_length_range 10-50 --e_value 0-$evalue_sl --seq_files $db_fasta --force

cat $outdir/BLAST_sl/selected/e10/* > $outdir/BLAST_sl/selected/e10.blast
awk 'BEGIN{OFS="\t"}!/^#/{if($11<=10){print}}' $outdir/BLAST_sl/selected/e10.blast | sponge $outdir/BLAST_sl/selected/e10.blast
awk 'BEGIN{OFS="\t"}!/^#/{if($11<=1){print}}' $outdir/BLAST_sl/selected/e10.blast | sponge $outdir/BLAST_sl/selected/e1.blast
awk 'BEGIN{OFS="\t"}!/^#/{if($11<=0.1){print}}' $outdir/BLAST_sl/selected/e10.blast | sponge $outdir/BLAST_sl/selected/e0.1.blast
awk 'BEGIN{OFS="\t"}!/^#/{if($11<=1e-3){print}}' $outdir/BLAST_sl/selected/e10.blast | sponge $outdir/BLAST_sl/selected/e0.001.blast
awk 'BEGIN{OFS="\t"}!/^#/{if($11<=1e-5){print}}' $outdir/BLAST_sl/selected/e10.blast | sponge $outdir/BLAST_sl/selected/e0.00001.blast


####################################################################
ruby $getUpstreamGff -i $outdir/BLAST_cds/blast_result --seq $infile --subject_file $db_fasta --query_start_max $query_start_max > $outdir/subject.gff

cut -f 9 $outdir/subject.gff | sort | uniq | sed 's/ID=//' > $outdir/subject.list
sort $outdir/subject.gff | sponge $outdir/subject.gff
sort $outdir/subject.list | sponge $outdir/subject.list


for i in 10 1 0.1 0.001 0.00001; do
	new_outdir=$outdir/e$i
	mkdir $new_outdir
	ruby $filter_BLAST_sl --gff $outdir/subject.gff -b $outdir/BLAST_sl/selected/e$i.blast --sl_seq $sl_file > $outdir/BLAST_sl/selected/filtered.e$i.blast
	ruby $blast2gff -i $outdir/BLAST_sl/selected/filtered.e$i.blast > $new_outdir/sl.gff
	bedtools intersect -a $outdir/subject.gff -b $new_outdir/sl.gff -s -F 1 > $new_outdir/bedtoolsIntersect
#cut -f 9 $outdir/sl.gff | sed 's/ID=//'| sort | uniq | join - $outdir/subject.gff | sort -k9 | cut -d " " -f 9| sed 's/ID=//' > $outdir/sl.list
	ruby $extract_slList --sl_gff $new_outdir/sl.gff --subject_gff $outdir/subject.gff > $new_outdir/sl.list
done


####################################################################
rm -rf $db_fasta


