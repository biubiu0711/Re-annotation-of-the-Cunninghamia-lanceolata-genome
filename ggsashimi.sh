gg=/home/qq/lwb/software/ggsahimi/ggsashimi.py
tsv=/home/qq/lwb/data/final_latest_gtf/ggsashimi/2nd/26owcw.csv
tsv=/home/qq/lwb/data/final_latest_gtf/ggsashimi/2nd/2owcw.csv
#region='Cula01:475285053-475356208'
#475287527-475325651:475325728-475325882
#475287527-475325651 475325728-475325882 475287527-475325882
gtf=/home/qq/lwb/data/final_latest_gtf/2nd_gtf/cf.2nd.gtf
#gene='MSTRG.3216'
#region='Cula01:475287441-475326011' #M3216
gene='MSTRG.13057'
region='Cula04:47816290-47816740'
region='Cula04:47812869-47816740'
gene='MSTRG.10301'
region='Cula03:289733357-289734494'



$gg -g $gtf -b $tsv -c ${region} -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 --fix-y-scale -o $gene -F pdf -s ANTISENSE
