import os
from Bio import SeqIO
from Bio import Entrez
import primer3
from pandas import DataFrame
from pandas import read_csv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML


def readFASTA(FASTA):
    Index = []
    Seq = []
    n = 0
    for i in FASTA:
        if '>' in i and n == 0:
            txt = ''
            Index.append(i.replace('>', '').replace('\n', ''))  # replace('\n','')
            n += 1
        elif '>' in i and n > 0:
            Seq.append(txt)
            txt = ''
            Index.append(i.replace('>', '').replace('\n', ''))
        else:
            txt += i.strip()
    Seq.append(txt)
    return Seq, Index


def extract_pairs(out_address):
    primer3_result_df = read_csv(out_address, header=0, index_col=0)
    left_primer = list(primer3_result_df.loc['PRIMER_LEFT_SEQUENCE'].values)
    right_primer = list(primer3_result_df.loc['PRIMER_RIGHT_SEQUENCE'].values)
    primer_list = []
    while left_primer:
        f = left_primer.pop(0)
        r = right_primer.pop(0)
        primer_list.append([f, r])
    return primer_list


def get_target(email: str, gene_id: dict, gene_name=''):
    Entrez.email = email
    if not gene_name:
        filename = gene_id['id']
    else:
        filename = gene_name
    if not os.path.isfile(filename):
        net_handle = Entrez.efetch(db=gene_id['db'], id=gene_id['id'], rettype=gene_id['rettype'],
                                   retmode=gene_id['retmode'])
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print(filename + ' Data Saved.')
    else:
        print('Data found in given local address.')
    print('Parsing...')
    record = SeqIO.read(filename, 'gb')
    seq = str(record.seq)
    return record, seq


def design_primer(seq_args, global_args, out_address):
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
    print('Succeeded. A total of ', primer3_result['PRIMER_PAIR_NUM_RETURNED'], 'primer pair(s) are designed.',
          '\nFor forward primer: ', primer3_result['PRIMER_LEFT_EXPLAIN'], '\nFor reverse primer: ',
          primer3_result['PRIMER_RIGHT_EXPLAIN'], '\nFor primer pairs: ', primer3_result['PRIMER_PAIR_EXPLAIN'])
    primer3_result_table_dict = {}
    for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        primer_id = str(i)
        for key in primer3_result:
            if '_' + primer_id + '_' in key:
                info_tag = key.replace("_" + primer_id, "")
                try:
                    primer3_result_table_dict[info_tag]
                except:
                    primer3_result_table_dict[info_tag] = []
                finally:
                    primer3_result_table_dict[info_tag].append(primer3_result[key])
    index = []
    for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        index.append("PRIMER_PAIR_" + str(i + 1))
    primer3_result_df = DataFrame(primer3_result_table_dict, index=index)
    primer3_result_df = primer3_result_df.T
    primer3_result_df.to_csv(out_address)
    print('Primer pair(s) Designed: \n', primer3_result_df)
    return primer3_result_df


def extract_pairs(out_address):
    primer3_result_df = read_csv(out_address, header=0, index_col=0)
    left_primer = list(primer3_result_df.loc['PRIMER_LEFT_SEQUENCE'].values)
    right_primer = list(primer3_result_df.loc['PRIMER_RIGHT_SEQUENCE'].values)
    primer_list = []
    while left_primer:
        f = left_primer.pop(0)
        r = right_primer.pop(0)
        primer_list.append([f, r])
    return primer_list


def blastn(query_address: str, db_address: str, out_address1: str = '', evalue=0.001, identity=18, task='blastn',
           dust='yes'):
    blastn_cline = NcbiblastnCommandline(query=query_address, db=db_address, evalue=evalue, outfmt=5, out=out_address1,
                                         task=task, dust=dust)
    stout, stderr = blastn_cline()
    result_handle = open(out_address1)
    blast_record = NCBIXML.read(result_handle)
    result_handle.close() # 关闭文件流
    e_value_thresh = evalue  # set E_value or other parameter and judge if exist
    identities = identity  # set identity for alignments,for primer design:length of primer-2 is recommended
    count = 0  # count number of blast hits
    # 修改列表，用来存储所有符合条件的结合位点详细信息
    hit_locations = [] 
    
    for alignment in blast_record.alignments:
        chrom_name = alignment.title.split()[-1] if ' ' in alignment.title else alignment.title
        
        for hsp in alignment.hsps:
            if hsp.expect <= e_value_thresh and hsp.identities >= identities:
                
                # 准确记录比对在基因组上的绝对起点和终点（不管正负链，确保 start < end）
                if hsp.sbjct_start <= hsp.sbjct_end:
                    strand = '+'
                    start_pos = hsp.sbjct_start
                    end_pos = hsp.sbjct_end
                else:
                    strand = '-'
                    start_pos = hsp.sbjct_end
                    end_pos = hsp.sbjct_start
                
                hit_locations.append({
                    'chr': chrom_name,
                    'start': start_pos,
                    'end': end_pos,       
                    'strand': strand
                })
                
    return hit_locations


def check_pcr_products(f_hits, r_hits, max_product_size=5000):
    """
    根据F和R在基因组上的所有结合位点，交叉计算能扩增出多少个潜在条带
    """
    products = []
    for f in f_hits:
        for r in r_hits:
            # 条件 1：必须在同一条染色体上
            if f['chr'] != r['chr']:
                continue
                
            p_size = 0
            
            # 条件 2：方向必须面对面，且计算包含引物全长的真实 PCR 产物大小
            # 情况 A：F 在正链(+)，R 在负链(-) -> F 在左边，R 在右边
            if f['strand'] == '+' and r['strand'] == '-':
                if f['start'] < r['end']:
                    # 真实长度 = 右侧引物的终止端 - 左侧引物的起始端 + 1
                    p_size = r['end'] - f['start'] + 1
                    
            # 情况 B：F 在负链(-)，R 在正链(+) -> R 在左边，F 在右边
            elif f['strand'] == '-' and r['strand'] == '+':
                if r['start'] < f['end']:
                    # 真实长度 = 右侧引物的终止端 - 左侧引物的起始端 + 1
                    p_size = f['end'] - r['start'] + 1
            
            # 条件 3：产物长度在常规PCR可扩增范围内
            if 0 < p_size <= max_product_size:
                products.append({
                    'chr': f['chr'],
                    'size': p_size
                })
                
    return products