#! /bin/bash


########################################################
calculate_overlaps=~/tools/self_bao_cun/others/calculate_overlaps.rb
filter_blast_result=~/huche/help/Jan/SL/scripts/filter_blast_result.rb


########################################################
while [ $# -gt 0 ]; do
	case $1 in
		--prefix)
			prefix=$2
			shift
			;;
		--no_exon)
			is_no_exon=true
			;;			
	esac
	shift
done

[ -z $prefix ] && echo "prefix has to be given! Exiting ......" && exit ;


########################################################
ruby2.1 $filter_blast_result --indir BLAST_SL/ori-2/ --outdir BLAST_SL/selected-2/CDS_upstream300_bl2seq_D1_E0-0.00001 --action cp --blast_outfmt 1 --aligned_length_range 10-50 --e_value 0-0.00001 --seq_files /home/sswang/huche/help/Jan/SL/data/genomic_data/Bschlosseri/result_sequences/Bschlosseri.CDS_upstream300.fasta --force

ruby2.1 $filter_blast_result --indir BLAST_SL/ori-2/ --outdir BLAST_SL/selected-2/CDS_upstream300_bl2seq_D1_E0-0.001 --action cp --blast_outfmt 1 --aligned_length_range 10-50 --e_value 0-0.001 --seq_files /home/sswang/huche/help/Jan/SL/data/genomic_data/Bschlosseri/result_sequences/Bschlosseri.CDS_upstream300.fasta --force

cat BLAST_SL/selected-2/CDS_upstream300_bl2seq_D1_E0-0.001/* > BLAST_SL/selected-2/CDS_upstream300_bl2seq_D1_E0-0.001.blast8


########################################################
cd gene_lists
[ ! -d bl2seq ] && mkdir bl2seq;
mv * bl2seq;
mkdir -p blastn/1e-5
ls ../BLAST_SL/selected-2/CDS_upstream300_bl2seq_D1_E0-0.00001/ > blastn/1e-5/genes_with_SLs.all.e0-0.00001.AlnLength10-50.list
[ -d bl2seq/blastn ] && rm -rf bl2seq/blastn
cd -


cd SL-parent_pairs/
[ ! -d bl2seq ] && mkdir bl2seq;
mv * bl2seq;
mkdir -p blastn/1e-5
cd blastn/1e-5
ruby $calculate_overlaps --i1 ../../bl2seq/final_pairs --i2 ../../../gene_lists/blastn/1e-5/genes_with_SLs.all.e0-0.00001.AlnLength10-50.list --show 12 --content 1 --no_report > final_pairs
ruby $calculate_overlaps --i1 ../../bl2seq/*yn00 --i2 final_pairs --f2 3 --show 12 --content 1 --no_report > ${prefix}_pair_yn00
cd ../../
[ -d bl2seq/blastn ] && rm -rf bl2seq/blastn
cd ../


if [ -z $is_no_exon ]; then
	cd exon_counts
	[ ! -d bl2seq ] && mkdir bl2seq;
	mv * bl2seq
	mv bl2seq/${prefix}.exon_counts ./
	mkdir -p blastn/1e-5
	ruby $calculate_overlaps --i1 ./*.exon_counts --i2 ../gene_lists/blastn/1e-5/genes_with_SLs.all.e0-0.00001.AlnLength10-50.list --show 12 --content 1 --no_report > blastn/1e-5/${prefix}_SL.all.exon_counts
	[ -d bl2seq/blastn ] && rm -rf bl2seq/blastn
	cd -
fi


