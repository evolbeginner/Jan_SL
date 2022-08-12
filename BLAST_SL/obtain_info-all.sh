for i in blast_res/*/*; do

	if grep 'kmer0' <<< $i >/dev/null; then
		bash /mnt/bay3/sswang/lian_xi/SL/scripts/BLAST_SL/get_distr.sh $i/blast_1st | sed 's/\t/:/' | ruby ~/tools/self_bao_cun/others/transpose.rb -i - > $i/info-all;
		continue
	fi
	
	cat ${i/kmer*/kmer0}/info-all > $i/info-all;

	bash /mnt/bay3/sswang/lian_xi/SL/scripts/BLAST_SL/get_distr.sh $i/blast_1st | sed 's/\t/:/' | ruby ~/tools/self_bao_cun/others/transpose.rb -i - >> $i/info-all

	awk 'BEGIN{OFS="\t"}{if(NR==1){a="R"}; if(NR==2){a="E"}; print a, $0}' $i/info-all | sponge $i/info-all

	no_genes=`wc -l $i/blast_1st | awk '{print $1}'`

	echo $no_genes >> $i/info-all

	tac $i/info-all | sponge $i/info-all
	cat $i/info-all

done
