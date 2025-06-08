genome=/home/qq/lwb/genome/Lachesis_assembly_changed.fa
input_fa=/home/qq/lwb/data/chinafir_data/rnaseq
input_fa_3rd=/home/qq/lwb/data/chinafir_data/longreads/fasta/hq
#gtf=/home/qq/lwb/data/chinafir_data/new_gtf/lr2rmats/lr_merge_orf.gtf
gtf=/home/qq/lwb/data/final_latest_gtf/2nd_gtf/cf.2nd.gtf
qualitative_output=/home/qq/lwb/data/final_latest_gtf/suppa/qualitative
quantitative_output=/home/qq/lwb/data/final_latest_gtf/suppa/quantitative
species=cff

source /home/qq/lwb/anaconda3/bin/activate suppa

if [ ! -d $qualitative_output ];then
	mkdir $qualitative_output
	suppa.py generateEvents -i $gtf -o $qualitative_output/$species.events -e SE SS MX RI FL -f ioe
	awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' $qualitative_output/*.ioe > $qualitative_output/$species.all.events.ioe
	wc -l $qualitative_output/$species.all.events.ioe
fi

if [ ! -d $quantitative_output ];then
	mkdir $quantitative_output
fi


if [ ! -d $quantitative_output/salmon_index ];then
	gffread $gtf -g $genome -w $quantitative_output/$species.trans.fa
	trans_fa=$quantitative_output/$species.trans.fa
	mkdir $quantitative_output/salmon_index
	salmon index -t $trans_fa -i $quantitative_output/salmon_index/$species.transcripts.salmon.index
fi

#suppa.py generateEvents -i $gtf -o $qualitative_output/$species.events -e SE SS MX RI FL -f ioe
#awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' $qualitative_output/*.ioe > $qualitative_output/$species.all.events.ioe
#wc -l $qualitative_output/$species.all.events.ioe


#for group in 2nd 3rd;do
for group in 3rd;do
	echo "$group start"
	mkdir $quantitative_output/$group
	if [ ! -d $quantitative_output/$group/salmon_output ];then
		mkdir $quantitative_output/$group/salmon_output
		mkdir $quantitative_output/$group/salmon_output/log
		index=$quantitative_output/salmon_index/$species.transcripts.salmon.index
		if [ $group == '2nd' ];then
			sample_name_list=(Cla-CK-Rep1_L3_175A75 Cla-CK-Rep2_L3_318X18 Cla-2h-CW-Rep1_L3_319X19 Cla-2h-CW-Rep2_L3_320X20 Cla-2h-OW-Rep1_L3_321X21 Cla-2h-OW-Rep2_L3_322X22 Cla-26h-CW-Rep1_L3_323X23 Cla-26h-CW-Rep2_L3_324X24 Cla-26h-OW-Rep1_L3_325X25 Cla-26h-OW-Rep2_L3_326X26 Cla-74h-CW-Rep1_L3_328X28 Cla-74h-CW-Rep2_L3_329X29 Cla-74h-OW-Rep1_L4_330X30 Cla-74h-OW-Rep2_L4_332X32)
		elif [ $group == '3rd' ];then
			sample_name_list=(fir_ck_rep1 fir_ck_rep2 fir_2h_cw_rep1 fir_2h_cw_rep2 fir_2h_ow_rep1 fir_2h_ow_rep2 fir_26h_cw_rep1 fir_26h_cw_rep2 fir_26h_ow_rep1 fir_26h_ow_rep2 fir_74h_cw_rep1 fir_74h_cw_rep2 fir_74h_ow_rep1 fir_74h_ow_rep2)
		fi
 		#for i in Cla-CK-Rep1_L3_175A75 Cla-CK-Rep2_L3_318X18 Cla-2h-CW-Rep1_L3_319X19 Cla-2h-CW-Rep2_L3_320X20 Cla-2h-OW-Rep1_L3_321X21 Cla-2h-OW-Rep2_L3_322X22 \
 	       #Cla-26h-CW-Rep1_L3_323X23 Cla-26h-CW-Rep2_L3_324X24 Cla-26h-OW-Rep1_L3_325X25 Cla-26h-OW-Rep2_L3_326X26 Cla-74h-CW-Rep1_L3_328X28 \
 	       #Cla-74h-CW-Rep2_L3_329X29 Cla-74h-OW-Rep1_L4_330X30 Cla-74h-OW-Rep2_L4_332X32;do
	       for i in "${sample_name_list[@]}";do
		if [ $group == '2nd' ];then
			sample=${i%%_*}
	 		salmon quant -i $index -l A --validateMappings --gcBias -1 $input_fa/${i}.R1.fastq.gz -2 $input_fa/${i}.R2.fastq.gz -p 10 -o $quantitative_output/$group/salmon_output/${sample} > $quantitative_output/$group/salmon_output/log/${sample}_salmon.log 2>&1
		elif [ $group == '3rd' ];then
			salmon quant -i $index -l A --validateMappings --gcBias -r $input_fa_3rd/${i}_hq.fasta  -p 10 -o $quantitative_output/$group/salmon_output/${i} > $quantitative_output/$group/salmon_output/log/${i}_salmon.log 2>&1
		fi
	 	done
		cd $quantitative_output/$group/salmon_output
		if [ $group == '2nd' ];then
			multipleFieldSelection.py -i Cla*/quant.sf -k 1 -f 4 -o $quantitative_output/$group/salmon_output/iso_tpm.txt
		elif [ $group == '3rd' ];then
			multipleFieldSelection.py -i fir*/quant.sf -k 1 -f 4 -o $quantitative_output/$group/salmon_output/iso_tpm.txt
		fi
		perl -alne '{/(\|.*\|)\t/; ;s/$1//g;s/\|//g;print}' $quantitative_output/$group/salmon_output/iso_tpm.txt > $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt
	fi

	suppa.py psiPerEvent -i $qualitative_output/$species.all.events.ioe -e $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt -o $quantitative_output/$group/project_events 1>psiPerEvent_log.txt

	mkdir $quantitative_output/$group/group
	cut -f 1-3 $quantitative_output/$group/project_events.psi > $quantitative_output/$group/group/ck.psi
	cut -f 1-3 $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt > $quantitative_output/$group/group/ck.tpm

	cut -f 1,4-5 $quantitative_output/$group/project_events.psi > $quantitative_output/$group/group/2cw.psi
	cut -f 1,4-5 $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt > $quantitative_output/$group/group/2cw.tpm

	cut -f 1,6-7 $quantitative_output/$group/project_events.psi > $quantitative_output/$group/group/2ow.psi
	cut -f 1,6-7 $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt > $quantitative_output/$group/group/2ow.tpm
	
	cut -f 1,8-9 $quantitative_output/$group/project_events.psi > $quantitative_output/$group/group/26cw.psi
	cut -f 1,8-9 $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt > $quantitative_output/$group/group/26cw.tpm

	cut -f 1,10-11 $quantitative_output/$group/project_events.psi > $quantitative_output/$group/group/26ow.psi
	cut -f 1,10-11 $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt > $quantitative_output/$group/group/26ow.tpm

	cut -f 1,12-13 $quantitative_output/$group/project_events.psi > $quantitative_output/$group/group/74cw.psi
	cut -f 1,12-13 $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt > $quantitative_output/$group/group/74cw.tpm

	cut -f 1,14-15 $quantitative_output/$group/project_events.psi > $quantitative_output/$group/group/74ow.psi
	cut -f 1,14-15 $quantitative_output/$group/salmon_output/iso_tpm_formatted.txt > $quantitative_output/$group/group/74ow.tpm

	mkdir $quantitative_output/$group/result
	mkdir $quantitative_output/$group/result/diff_csv
	
	# for i in ck 2cw 26cw 74cw;do
	for i in ck 2ow 26ow 74ow;do
		# for j in 2cw 2ow 26cw 26ow 74cw 74ow;do
		for j in 2cw 2ow 26cw 26ow 74cw 74ow;do
			if [ "$i" == "$j" ];then
				echo 'same group'
				continue
			elif [ $i == '2ow' -a $j != '2cw' ];then
				echo 'not 2c vs 2o'
       				continue
			elif [ $i == '26ow' -a $j != '26cw' ];then
	      			echo 'not 26c vs 26o'
       	                 	continue
                	elif [ $i == '74ow' -a $j != '74cw' ];then
                        	echo 'not 74c vs 74o'
                        	continue
                	fi
			echo $i' vs '$j >> $quantitative_output/$group/result/statistic_0.05
			echo $i' vs '$j >> $quantitative_output/$group/result/statistic_0.01
			suppa.py diffSplice -m empirical -gc -i $qualitative_output/$species.all.events.ioe --save_tpm_events -p $quantitative_output/$group/group/$j.psi $quantitative_output/$group/group/$i.psi -e $quantitative_output/$group/group/$j.tpm $quantitative_output/$group/group/$i.tpm -o $quantitative_output/$group/result/$j'vs'$i'_diffSplice'
			cat $quantitative_output/$group/result/$j'vs'$i'_diffSplice.dpsi' | perl -alne '{print if $F[2] <0.05}' | wc >> $quantitative_output/$group/result/statistic_0.05
			cat $quantitative_output/$group/result/$j'vs'$i'_diffSplice.dpsi' | perl -alne '{print if $F[2] <0.01}' | wc >> $quantitative_output/$group/result/statistic_0.01
		done
	done
	echo "$group  finish"
done
echo 'all finish'
