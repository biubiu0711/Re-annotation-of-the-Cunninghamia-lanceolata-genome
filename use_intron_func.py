import intron_function as infunc
import os

input_dir = r'J:\data\annotation_final_edition\日本晴基因组'
gtf = r'NIP-T2T.gtf'
gtf_name_prefix = 'rice'
output_1 = r'J:\data\annotation_final_edition\日本晴基因组'

reverse = False
gap = 1000

for extract_type in ['exon', 'intron']:
    output_dir = output_1 + '\\' + extract_type
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    for site in ['first', 'last']:
        input_gtf = input_dir + '\\' + gtf

        infunc.prepare(input_gtf, output_gtf=output_dir + '\\' + gtf_name_prefix + '.' + extract_type + '.gtf')

        infunc.extract_longest_trans_represent_gene(input_gtf=output_dir + '\\' + gtf_name_prefix + '.' + extract_type + '.gtf',
                                                    output_gtf=output_dir + '\\' + gtf_name_prefix + '.represent_gene.gtf')
        if extract_type == 'intron':
            infunc.extract_intron(input_gtf=output_dir + '\\' + gtf_name_prefix + '.represent_gene.gtf',
                                  output_gtf=output_dir + '\\' + gtf_name_prefix + '.' + extract_type + '.gtf',
                                  reverse=reverse)

        infunc.extract_proximal_distal_intron(input_gtf=output_dir + '\\' + gtf_name_prefix + '.' + extract_type + '.gtf',
                                              output_gtf=output_dir + '\\' + gtf_name_prefix + '.proximal_distal.' + extract_type + '.gtf', gap=gap,
                                              extract_type=extract_type)

        infunc.extract_proximal_distal_distance_intron(input_gtf=output_dir + '\\' + gtf_name_prefix  + '.' + extract_type + '.gtf',
                                                       output_gtf=output_dir + '\\' + gtf_name_prefix + '.proximal_distal_distance.' + extract_type + '.gtf',
                                                       gap=gap, extract_type=extract_type)

        infunc.extract_FirstOrLast_IntronOrExon_lenth(input_gtf=output_dir + '\\' + gtf_name_prefix + '.' + extract_type + '.gtf',
                                                      output_gtf=output_dir + '\\' + gtf_name_prefix + '.' + site + '_' + extract_type + '_length.gtf',
                                                      extract_type=extract_type, site=site)

#input_gtf = r'J:\data\bioinformatic\chinafir\tama\new_pipeline\czs_uncorrect_no_filter\czs_nopolish_nofilter.represent.final.gtf'

#infunc.count_overlap_gene(input_gtf)

