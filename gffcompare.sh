ref_gtf=/home/qq/lwb/genome/chinafir.gtf
tama_gtf=/home/qq/work/lwb/final_gtf/cf_final.final.gtf
str_gtf=/home/qq/lwb/data/chinafir_data/new_gtf/align/stringtie_not_ref_gtf/stringtied_merged.gtf
#lr_gtf=/home/qq/lwb/data/chinafir_data/new_gtf/lr2rmats/np_gtf/cf_np.lr.gtf

gffcompare -r $tama_gtf -o cf $str_gtf
