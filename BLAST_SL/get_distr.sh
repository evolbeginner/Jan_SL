###################################################
infile=$1


###################################################
echo -ne "5-\t"
awk 'function abs(v) {return v < 0 ? -v : v} {gap_ext=abs(abs($8-$7)-abs($10-$9))-$6; a=$4-$5*3-$6*5-gap_ext*2; if(a<=5){c++}}END{print c/NR}' $infile


###################################################
function output_interval(){
	for a in 5-8 9-10 11-12 13-14; do
		i=`echo $a | awk -F"-" '{print $1}'`
		j=`echo $a | awk -F"-" '{print $2}'`
		export i; export j
		echo -ne "$a\t"
		awk 'function abs(v) {return v < 0 ? -v : v} {gap_ext=abs(abs($8-$7)-abs($10-$9))-$6; a=$4-$5*3-$6*5-gap_ext*2; if(a>=ENVIRON["i"] && a<=ENVIRON["j"]){c++}}END{print c/NR}' $infile
	done
}


###################################################
for i in `seq 6 14`; do
	export i=$i
	echo -ne "$i\t"
	#awk 'function abs(v) {return v < 0 ? -v : v} {gap_ext=abs($8-$7-($10-$9)-$6); a=$4-$5*3-$6*5-gap_ext*2; if(a>=ENVIRON["i"] && a<=ENVIRON["i"]){c++}}END{print c/NR}' $infile
	awk 'function abs(v) {return v < 0 ? -v : v} {gap_ext=abs(abs($8-$7)-abs($10-$9))-$6; a=$4-$5*3-$6*5-gap_ext*2; if(a>=ENVIRON["i"] && a<=ENVIRON["i"]){c++}}END{print c/NR}' $infile
done


# 11-12, 13-14
# output_interval

# 15+
export i=15
echo -ne "$i+\t"
awk 'function abs(v) {return v < 0 ? -v : v} {gap_ext=abs(abs($8-$7)-abs($10-$9))-$6; a=$4-$5*3-$6*5-gap_ext*2; if(a>=ENVIRON["i"]){c++}}END{print c/NR}' $infile


