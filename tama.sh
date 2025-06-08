tama=/home/lfgu/lwb/software/tama-master
bam_1=/home/lfgu/lwb/chinafir/tama/fir_output/new_pipeline/bam/cf.allsample.sorted.bam
genome_1=/home/lfgu/lwb/chinafir/genome/Lachesis_assembly_changed.fa
prefix=cf_newpi
output=/home/lfgu/lwb/chinafir/tama/fir_output/new_pipeline
samtools=/home/lfgu/lwb/software/samtools-0.1.13/samtools
protein_db=/home/lfgu/lwb/chinafir/UniProt/uniprot_sprot.fasta
db_name=uniprot
blast_db_dir=/home/lfgu/lwb/chinafir/UniProt


source /home/lfgu/lwb/software/conda/bin/activate tama_lwb

#	$samtools view -Sb $sam_1/$prefix'_'$i'.sam' -o $output/$i/$i'.bam'
#        $samtools sort $output/$i/$i'.bam' $output/$i/$i'_sort'
#        $samtools index $output/$i/$i'_sort.bam'
#	python2 $tama/tama_collapse.py -b BAM -s $output/$i/$i'_sort.bam' -f $genome_1 -p $prefix
if [ ! -d $output/01-collapse ];then
	mkdir $output/01-collapse
	python2 $tama/tama_collapse.py -b BAM -rm low_mem -s $bam_1 -f $genome_1 -p $prefix -d merge_dup -x no_cap -m 0 -a 0 -z 0 -sj sj_priority -lde 30 -sjt 30
	echo 'collapse finish'
	echo '1' > $output/01-collapse/filelist.txt
	find $output/01-collapse -name *trans_read.bed | xargs >>  $output/01-collapse/filelist.txt
	echo 'trans_read' >> $output/01-collapse/filelist.txt
	cat $output/01-collapse/filelist.txt | xargs > filelist
	sed -i 's/ /\t/g' filelist
fi

##python2 $tama/tama_merge.py -f file_list -p chuizhishan -m 0 -a 0 -z 0 -d merge_dup
filelist=$output/01-collapse/filelist

if [ ! -d $output/02-readsupport ];then
	mkdir $output/02-readsupport
	cd $output/02-readsupport
        python2 $tama/tama_go/read_support/tama_read_support_levels.py -f $filelist -m $output/01-collapse/chinafir_merge.txt -o $prefix
        echo 'read-support finish'
fi

if [ ! -d $output/03-readsfilter ];then
	mkdir $output/03-readsfilter
	cd $output/03-readsfilter
	python2 $tama/tama_go/filter_transcript_models/tama_remove_single_read_models_levels.py -b $output/01-collapse/'chinafir.bed' -r $output/02-readsupport/$prefix'_read_support.txt' -o $prefix'.filtered'
	echo 'filter finish'
fi

if [ ! -d $output/04-orfpredict ];then
	mkdir $output/04-orfpredict
	cd $output/04-orfpredict
	#bedtools getfasta -name -split -s -fi $genome_1 -bed $output/03-readsfilter/$prefix'.filtered.bed' -fo $output/04-orfpredict/$prefix'.tama.fa'
	bedtools getfasta -name -split -s -fi ${genome_1} -bed $output/03-readsfilter/$prefix'.filtered.bed' -fo $output/04-orfpredict/$prefix'.tama.fa'
	orf=${prefix}'.orf.fasta'
	python2 ${tama}/tama_go/orf_nmd_predictions/tama_orf_seeker.py -f $output/04-orfpredict/$prefix'.tama.fa' -o $output/04-orfpredict/${orf}

	conda deactivate
	#makeblastdb -in ${protein_db} -dbtype prot -out ${blast_db_dir}/${db_name}
	blast_out=${prefix}'.blast'
	blastp -evalue 1e-10 -num_threads 40 -db ${blast_db_dir}/${db_name} -query $output/04-orfpredict/${orf} -out $output/04-orfpredict/${blast_out}

	source /home/lfgu/lwb/software/conda/bin/activate tama_lwb
	orf_blast_out=${prefix}'.orf.blast.out'
	python2 ${tama}/tama_go/orf_nmd_predictions/tama_orf_blastp_parser.py -b ${blast_out} -o ${orf_blast_out}
	orf_bed=${prefix}'.orf.bed'
	python2 ${tama}/tama_go/orf_nmd_predictions/tama_cds_regions_bed_add.py -p ${orf} -a $output/03-readsfilter/$prefix'.filtered.bed' -f $output/04-orfpredict/$prefix'.tama.fa' -o $output/04-orfpredict/${orf_bed}
	echo 'orf predict finish'
fi

if [ ! -d $output/05-keep_1_representative_model ];then
        mkdir $output/05-keep_1_representative_model
        cd $output/05-keep_1_representative_model
        python2 $tama/tama_go/filter_transcript_models/tama_filter_primary_transcripts_orf.py -b $output/04-orfpredict/$prefix'.orf.bed' -o $output/05-keep_1_representative_model/$prefix'.represent.bed'
	python2 $tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_orf_nmd.py $output/05-keep_1_representative_model/$prefix'.represent.bed' $output/05-keep_1_representative_model/$prefix'.final.gtf'
        echo 'finish'
fi

exit


if [ ! -d $output/06-remove_fragment ];then
	mkdir $output/06-remove_fragment 
	cd $output/06-remove_fragment
	python2 $tama/tama_go/filter_transcript_models/tama_remove_fragment_models.py -f $output/05-keep_1_representative_model/$prefix'.represent.bed' -o $output/06-remove_fragment/$prefix'.remove.fragment'
	echo 'remove fragment finish'
fi



if [ ! -d $output/07-reformat ];then
	mkdir $output/07-reformat
	cd $output/07-reformat
	python2 $tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_orf_nmd.py $output/05-keep_1_representative_model/$prefix'.represent.bed' $output/07-reformat/$prefix'.final.gtf'
	echo 'reformat finish'
fi
