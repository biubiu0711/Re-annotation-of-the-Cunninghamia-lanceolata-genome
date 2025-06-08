# Re-annotation-of-the-Cunninghamia-lanceolata-genome
Here are some scripts

![技术路线v3](https://github.com/user-attachments/assets/8fbb0418-4707-4edd-bfb5-193d479edf00)

## annotation
tama.sh + gffcompare.sh + transdecoder.sh + interproscan.sh

## differential expression
hisat2.sh + DE.R

## WGCNA
WGCNA.R

## AS
suppa.sh + rmats.sh + ggsashimi.sh  
suppa.sh use for Pacbio data
rmats.sh usr for Rna-seq data
ggsashimi.sh can plot target gene in rmats result

## APA
Chinafir.ipynb(part APA  

## gene feature and gene expression
intron_function.py include some tools use for analyze of the gene feature and gene expression,  
you can use use_intron_func.py to use it, input file is gtf file

## final
Chinafir.ipynb include all analysis and plotting scripts, in jupyter notebook format.
