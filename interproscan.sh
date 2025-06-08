pep=/home/qq/work/lwb/final_gtf/transdecoder/chinafir.trans.fasta.transdecoder.2.pep
interpro=/home/qq/lwb/software/interproscan/interproscan-5.67-99.0/interproscan.sh
output=/home/qq/work/lwb/final_gtf/interpro/result




bash $interpro -iprlookup -goterms -dp -appl Pfam -appl PRINTS -appl PANTHER -appl ProSiteProfiles -appl SMART -f TSV -i $pep -o $output/cf.final.fasta.ipscan

echo 'finish'

