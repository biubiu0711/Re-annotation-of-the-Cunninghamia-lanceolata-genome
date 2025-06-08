#!/bin/bash

input=/home/qq/lwb/data/final_latest_gtf/rmats/group
output=/home/qq/lwb/data/final_latest_gtf/rmats/result
rmat=/home/qq/lwb/anaconda3/envs/rmats-lwb/rMATS
rnainput=/home/qq/work/lwb/data/chinafir_data/rna_result/samtools_output
gtf=/home/qq/lwb/data/final_latest_gtf/2nd_gtf/cf.2nd.gtf

source /home/qq/lwb/anaconda3/bin/activate rmats-lwb

#for i in fir_ck.txt fir_2_c.txt fir_2_o.txt fir_26_c.txt fir_26_o.txt fir_74_c.txt fir_74_o.txt;do
for i in CK.txt 2OW.txt 26OW.txt 74OW.txt;do
        #echo $i
	for j in 2CW.txt 2OW.txt 26CW.txt 26OW.txt 74CW.txt 74OW.txt;do
		#echo $i 'vs'$j
		if [ "$i" == "$j" ];then
			echo 'same group'
			continue
		elif [ $i == '2OW.txt' -a $j != '2CW.txt' ];then
			echo 'not 2c vs 2o' 
			continue
		elif [ $i == '26OW.txt' -a $j != '26CW.txt' ];then
			echo 'not 26c vs 26o'
			continue
		elif [ $i == '74OW.txt' -a $j != '74CW.txt' ];then
			echo 'not 74c vs 74o'
			continue
		#elif [ $i == 'fir_74_o.txt' -o $i == 'fir_2_o.txt' -o $i == 'fir_26_o.txt' ];then
		#	continue
		fi
		i_1=${i%%.*}
		j_1=${j%%.*}
		dir=$j_1'vs'$i_1
		echo $dir
        	mkdir $output/$dir
       		mkdir $output/tmp
        	python $rmat/rmats.py --b1 $input/$j --b2 $input/$i --gtf $gtf -t paired --readLength 150 --nthread 20 --od $output/$dir --tmp $output/tmp > $output/$dir/log 2>&1
		for k in A3SS A5SS RI SE MXE;do
			file_name_1=$k'.MATS.JC.txt'
			file_name_2=$k'.MATS.JCEC.txt'
			replace_1=$k'.MATS.JC.awk.txt'
			replace_2=$k'.MATS.JCEC.awk.txt'
			awk -F '\t' '{if (NR != 1) $3=$2}1' $output/$dir/$file_name_1 > $output/$dir/$replace_1
			awk -F '\t' '{if (NR != 1) $3=$2}1' $output/$dir/$file_name_2 > $output/$dir/$replace_2
			sed -i 's/ /\t/g' $output/$dir/$replace_1
			sed -i 's/ /\t/g' $output/$dir/$replace_2
			sed -i 's/\"//g' $output/$dir/$replace_1
			sed -i 's/\"//g' $output/$dir/$replace_2
		done
		rm -r $output/tmp
        	echo $dir ' finish'
	done
done

