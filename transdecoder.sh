td=/home/qq/lwb/software/TransDecoder-v5.7.1
gtf=/home/qq/work/lwb/final_gtf/tama.gtf
#trans_fa=/home/qq/lwb/data/chinafir_data/annotation/transcript/chinafir.transcripts.3rd.fa
#input=
output=/home/qq/work/lwb/final_gtf/test
#db=/home/qq/work/lwb/UniProtKB
db=/home/qq/work/lwb/UniRef
td_util=/home/qq/lwb/software/TransDecoder-v5.7.1/util
ge=/home/qq/lwb/genome/Lachesis_assembly_changed.fa

cd $output

$td_util/gtf_to_alignment_gff3.pl $gtf > $output/cf.trans.gff3

$td_util/gtf_genome_to_cdna_fasta.pl $gtf $ge > $output/chinafir.trans.fasta

trans_fa=$output/chinafir.trans.fasta

$td/TransDecoder.LongOrfs -t $trans_fa
fa=${trans_fa##*/}
echo $fa
#blastp -query $output/$fa'.transdecoder_dir/longest_orfs.pep' -db $db/unipro -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > $output/blastp.outfmt6
blastp -query $output/$fa'.transdecoder_dir/longest_orfs.pep' -db $db/uniref50 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 40 > $output/blastp.outfmt6


$td/TransDecoder.Predict -t $trans_fa --retain_blastp_hits $output/blastp.outfmt6


