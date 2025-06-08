def prepare(input_gtf, output_gtf):
    r = open(output_gtf, 'w')

    with open(input_gtf) as f:
        for line in f:
            line_rep = line.strip().split('\t')
            g_type = line_rep[2]
            info = line_rep[8].split('; ')
            info_new = info[0] + '; ' + info[1]
            line_rep[8] = info_new
            if g_type == 'transcript' or g_type == 'exon':
                r.write('\t'.join(line_rep) + '\n')

    r.close()

def extract_intron(input_gtf, output_gtf, reverse=True):
    # 输出为内含子的gtf
    r = open(output_gtf, 'w')

    intron_end = None

    with open(input_gtf) as f:
        for line in f:
            line_rep = line.strip().split('\t')
            strand = line_rep[6]
            if line_rep[2] == 'transcript':
                exon_c = 0
                r.write(line)
            elif line_rep[2] == 'exon':
                exon_c += 1
            else:
                continue
            if reverse:
                if strand == '+':
                    if exon_c == 1:
                        intron_star = int(line_rep[4]) + 1
                    elif exon_c > 1:
                        intron_end = int(line_rep[3]) - 1
                        line_rep[2] = 'intron'
                        line_rep[3], intron_star = str(intron_star), int(line_rep[4]) + 1
                        line_rep[4] = str(intron_end)
                        r.write('\t'.join(line_rep) + '\n')
                elif strand == '-':
                        if exon_c == 1:
                            intron_end = int(line_rep[3]) - 1
                        elif exon_c > 1:
                            intron_star = int(line_rep[4]) + 1
                            line_rep[2] = 'intron'
                            line_rep[4], intron_end = str(intron_end), int(line_rep[3]) - 1
                            line_rep[3] = str(intron_star)
                            r.write('\t'.join(line_rep) + '\n')
            elif not reverse:
                if exon_c == 1:
                    intron_star = int(line_rep[4]) + 1
                elif exon_c > 1:
                    intron_end = int(line_rep[3]) - 1
                    line_rep[2] = 'intron'
                    line_rep[3], intron_star = str(intron_star), int(line_rep[4]) + 1
                    line_rep[4] = str(intron_end)
                    r.write('\t'.join(line_rep) + '\n')


    r.close()

def extract_longest_trans_represent_gene(input_gtf, output_gtf):
    # 选出最长的转录本作为基因的代表
    trans_dit = {}
    gene_trans = {}

    with open(input_gtf) as f:
        for line in f:
            line_rep = line.strip().split('\t')
            if line_rep[2] == 'transcript':
                trans_id = line_rep[8].split('; ')[0].split('"')[1]
                gene_id = line_rep[8].split('; ')[1].split('"')[1]
                lenth = int(line_rep[4]) - int(line_rep[3])
                if gene_id in trans_dit:
                    if trans_dit[gene_id] < lenth:
                        trans_dit[gene_id] = lenth
                        gene_trans[gene_id] = trans_id
                else:
                    trans_dit[gene_id] = lenth
                    gene_trans[gene_id] = trans_id
            else:
                continue

    r = open(output_gtf, 'w')

    with open(input_gtf) as f:
        for line in f:
            line_rep = line.strip().split('\t')
            trans_id = line_rep[8].split('; ')[0].split('"')[1]
            gene_id = line_rep[8].split('; ')[1].split('"')[1]
            line_type = line_rep[2]
            if trans_id in gene_trans[gene_id]:
                r.write(line.strip() + '\n')
            # if line_type == 'gene' or line_type == 'transcript':
            #     r.write(line.strip() + '\n')
            #     gene_id = line_rep[8].split('; ')[1].split('"')[1]
            #     trans_longest = gene_trans[gene_id]
            #     continue
            # if trans_id == trans_longest:
            #     r.write(line.strip() + '\n')
    r.close()


def extract_proximal_distal_intron(input_gtf, output_gtf, gap=1000, extract_type='intron'):
    # 输入文件为extract_intron的输出文件
    # 提取启动子-近端、终止子-近端内含子，输出为GTF + ‘\t' + 此内含子为近端还是远端
    r = open(output_gtf, 'w')

    with open(input_gtf) as f:
        for line in f:
            line_rep = line.strip().split('\t')
            strand = line_rep[6]
            gene_type = line_rep[2]
            if gene_type == 'transcript':
                r.write(line.strip() + '\t' + 'site' + '\n')
                if strand == '+':
                    gene_start = int(line_rep[3])
                    gene_end = int(line_rep[4])
                elif strand == '-':
                    gene_start = int(line_rep[4])
                    gene_end = int(line_rep[3])
            elif gene_type == extract_type:
                if strand == '+':
                    proximal_site = int(line_rep[3])
                    distal_site = int(line_rep[4])
                    if proximal_site - gap <= gene_start:
                        r.write(line.strip() + '\t' + 'proximal' + '\n')
                    elif distal_site + gap >= gene_end:
                        r.write(line.strip() + '\t' + 'distal' + '\n')
                    else:
                        r.write(line.strip() + '\t' + 'none' + '\n')
                elif strand == '-':
                    proximal_site = int(line_rep[4])
                    distal_site = int(line_rep[3])
                    if proximal_site + gap >= gene_start:
                        r.write(line.strip() + '\t' + 'proximal' + '\n')
                    elif distal_site - gap <= gene_end:
                        r.write(line.strip() + '\t' + 'distal' + '\n')
                    else:
                        r.write(line.strip() + '\t' + 'none' + '\n')
    r.close()


def extract_proximal_distal_distance_intron(input_gtf, output_gtf, gap=1000, extract_type='intron'):
    # 输入文件为extract_intron的输出文件
    # 提取启动子-近端、终止子-近端内含子的距离，输出为GTF + ‘\t' + 启动子近端内含子距离 + '\t' + 终止子近端内含子距离
    r = open(output_gtf, 'w')
    gene_dict = {}

    with open(input_gtf) as f:
        for line in f:
            line_rep = line.strip().split('\t')
            strand = line_rep[6]
            gene_type = line_rep[2]
            if gene_type == 'transcript':
                gene_id = line_rep[8].split('; ')[0].split('"')[1]
                res = line.strip()
                gene_dict[gene_id] = {}
                gene_dict[gene_id]['transcript'] = res
                gene_dict[gene_id]['proximal_dis'] = None
                gene_dict[gene_id]['distal_dis'] = None
                if strand == '+':
                    gene_start = int(line_rep[3])
                    gene_end = int(line_rep[4])
                elif strand == '-':
                    gene_start = int(line_rep[4])
                    gene_end = int(line_rep[3])
            elif gene_type == extract_type:
                if strand == '+':
                    proximal_site = int(line_rep[3])
                    distal_site = int(line_rep[4])
                    if proximal_site - gap <= gene_start:
                        proximal_distance = proximal_site - gene_start
                        if not gene_dict[gene_id]['proximal_dis'] or gene_dict[gene_id][
                            'proximal_dis'] >= proximal_distance:
                            gene_dict[gene_id]['proximal_dis'] = proximal_distance
                    elif distal_site + gap >= gene_end:
                        distal_distance = gene_end - distal_site
                        if not gene_dict[gene_id]['distal_dis'] or gene_dict[gene_id]['distal_dis'] >= distal_distance:
                            gene_dict[gene_id]['distal_dis'] = distal_distance
                elif strand == '-':
                    proximal_site = int(line_rep[4])
                    distal_site = int(line_rep[3])
                    if proximal_site + gap >= gene_start:
                        proximal_distance = gene_start - proximal_site
                        if not gene_dict[gene_id]['proximal_dis'] or \
                                gene_dict[gene_id]['proximal_dis'] >= proximal_distance:
                            gene_dict[gene_id]['proximal_dis'] = proximal_distance
                    elif distal_site - gap <= gene_end:
                        distal_distance = distal_site - gene_end
                        if not gene_dict[gene_id]['distal_dis'] or \
                                gene_dict[gene_id]['distal_dis'] >= distal_distance:
                            gene_dict[gene_id]['distal_dis'] = distal_distance

    for key in gene_dict.keys():
        string = gene_dict[key]['transcript']
        if not gene_dict[key]['proximal_dis']:
            string += '\t' + 'None'
        else:
            string += '\t' + str(gene_dict[key]['proximal_dis'])
        if not gene_dict[key]['distal_dis']:
            string += '\t' + 'None'
        else:
            string += '\t' + str(gene_dict[key]['distal_dis'])
        r.write(string + '\n')

    r.close()


def extract_FirstOrLast_IntronOrExon_lenth(input_gtf, output_gtf, extract_type='intron', site='first'):
    # 提取转录本第一个或最后一个内含子或外显子的长度
    r = open(output_gtf, 'w')
    gene_dict = {}

    with open(input_gtf) as f:
        for line in f:
            line_rep = line.strip().split('\t')
            strand = line_rep[6]
            gene_type = line_rep[2]
            gene_id = line_rep[8].split('; ')[0].split('"')[1]
            if gene_type == 'transcript':
                res = line.strip()
                gene_dict[gene_id] = {}
                gene_dict[gene_id]['transcript'] = res
                gene_dict[gene_id]['intron_start'] = None
                gene_dict[gene_id]['length'] = None
            elif gene_type == extract_type:
                intron_length = int(line_rep[4]) - int(line_rep[3]) + 1
                if strand == '+':
                    intron_start = int(line_rep[3])
                    if gene_id in gene_dict:
                        if site == 'last':
                            if not gene_dict[gene_id]['intron_start'] or gene_dict[gene_id]['intron_start'] < intron_start:
                                gene_dict[gene_id]['intron_start'] = intron_start
                                gene_dict[gene_id]['length'] = intron_length
                        elif site == 'first':
                            if not gene_dict[gene_id]['intron_start'] or gene_dict[gene_id]['intron_start'] > intron_start:
                                gene_dict[gene_id]['intron_start'] = intron_start
                                gene_dict[gene_id]['length'] = intron_length
                    else:
                        print('ERROR')
                        break
                elif strand == '-':
                    intron_start = int(line_rep[4])
                    if gene_id in gene_dict:
                        if site == 'last':
                            if not gene_dict[gene_id]['intron_start'] or gene_dict[gene_id]['intron_start'] > intron_start:
                                gene_dict[gene_id]['intron_start'] = intron_start
                                gene_dict[gene_id]['length'] = intron_length
                        elif site == 'first':
                            if not gene_dict[gene_id]['intron_start'] or gene_dict[gene_id]['intron_start'] < intron_start:
                                gene_dict[gene_id]['intron_start'] = intron_start
                                gene_dict[gene_id]['length'] = intron_length
                    else:
                        print('ERROR')
                        break
    for key in gene_dict.keys():
        gene = gene_dict[key]['transcript'] + '\t' + str(gene_dict[key]['length'])
        r.write(gene + '\n')

    r.close()




def count_overlap_gene(input_gtf):
    # 检测基因注释中有无重复的基因
    count = 0
    gene = {}
    overlap_id = []

    with open(input_gtf) as f:
        for line in f:
            line_rep = line.strip().split('\t')
            line_type = line_rep[2]
            gene_id = line_rep[8].split('; ')[0].split('"')[1]
            chr = line_rep[0]
            gene_start = int(line_rep[3])
            gene_end = int(line_rep[4])
            strand = line_rep[6]
            if line_type != 'gene':
                continue
            if len(gene.keys()) == 0:
                gene[gene_id] = {'chr':chr, 'gene_start':gene_start, 'gene_end':gene_end, 'strand':strand}
            else:
                overlap = False
                for key in gene.keys():
                    if gene[key]['chr'] == chr and gene[key]['strand'] == strand:
                        if gene[key]['gene_start'] <= gene_start and gene[key]['gene_end'] >= gene_end:
                            count += 1
                            overlap_id.append(gene_id)
                            overlap = True
                            #print(gene_id + '  ' + key + ' overlap')
                        elif gene[key]['gene_start'] >= gene_start and gene[key]['gene_end'] <= gene_end:
                            count += 1
                            overlap_id.append(key)
                            overlap = True
                            #print(gene_id + '  ' + key + ' overlap')
                            #gene.pop(key)
                            gene[key] = {'chr':chr, 'gene_start':gene_start, 'gene_end':gene_end, 'strand':strand}
                if not overlap:
                    gene[gene_id] = {'chr': chr, 'gene_start': gene_start, 'gene_end': gene_end, 'strand': strand}
    print(len(gene.keys()))
    print(count)
    #print(overlap_id)


def stringtie_rename(input_gtf, output_gtf):
    # 将stringtie组装后基因名（MSTRG.*）改为gtf中的基因名
    r = open(output_gtf, 'w')

    with open(input_gtf) as f:
        for line in f:
            if line[0] == '#':
                continue
            line_rep = line.strip().split('\t')
            tmp_list = line_rep[8].split('; ')
            trans_id = tmp_list[1].split('"')[1]
            gene_id = trans_id.split('.')[0]
            line_rep[8] = 'gene_id "' + gene_id + '"; ' + '; '.join(tmp_list[1:])
            r.write('\t'.join(line_rep) + '\n')
    r.close()
