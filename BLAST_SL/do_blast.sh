for h in shuffled/*; do B=`basename $h`; for j in $h/*; do b=`basename $j`; echo "$B $b"; ruby /mnt/bay3/sswang/lian_xi/SL/results/Skawagutii/BLAST_SL/new/run_blast_in_batch_SL.rb --indir $j --outdir blast/$B/$b --force --cpu 16; ruby /mnt/bay3/sswang/lian_xi/SL/results/Skawagutii/BLAST_SL/new/run_blast_in_batch_SL.rb --indir $j --outdir blast/$B-gapped/$b --force --cpu 16 --with_gap; done ; done
