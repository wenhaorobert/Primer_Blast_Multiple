import os
import sys
import pandas as pd
import primer3
from PrimerBlast_Tools import Primer_Blast_Fuction as pb

db_address = '/home/wr/document/genome/NCBI-zm-5/zm5'
# db_address = '/home/wr/document/genome/NCBI-zm-5/zm5rna'
# db_address = '/home/wr/document/genome/MaizeGDB-W22-v2/w22dna'
print("现在在用的基因组数据库是", db_address)
out_dir = '/home/wr/github_software/Primer_Blast_Multiple/Output/'
    
def evaluate_primer_quality(seq, name="引物"):
    """
    利用 primer3 的底层热力学接口，计算引物物化参数：
    长度、GC含量、Tm值、发夹结构(Hairpin)和二聚体(Dimer)
    """
    length = len(seq)
    
    # 1. 计算 GC 含量
    gc_count = seq.count('G') + seq.count('C')
    gc_content = (gc_count / length) * 100 if length > 0 else 0
    
    # 2. 计算 Tm 值
    tm = primer3.calc_tm(seq)
    
    # 3. 计算发夹结构 (Hairpin) 的最高解链温度
    hairpin_res = primer3.calc_hairpin(seq)
    hairpin_tm = hairpin_res.tm
    
    # 4. 计算自身二聚体 (Self-Dimer) 的最高解链温度
    dimer_res = primer3.calc_homodimer(seq)
    dimer_tm = dimer_res.tm

    # 设立推荐指标预警
    warnings = []
    if not (18 <= length <= 27):
        warnings.append(f"长度异常({length}bp)")
    if not (40 <= gc_content <= 60):
        warnings.append(f"GC偏离({gc_content:.1f}%)")
    if not (55 <= tm <= 65):
        warnings.append(f"Tm偏离({tm:.1f}℃)")
    if hairpin_tm > 40:
        warnings.append(f"易形成发夹(Hairpin Tm: {hairpin_tm:.1f}℃)")
    if dimer_tm > 40:
        warnings.append(f"易形成二聚体(Dimer Tm: {dimer_tm:.1f}℃)")
        
    quality_status = "✅ 优秀" if len(warnings) == 0 else f"⚠️ 建议注意 ({', '.join(warnings)})"
    
    return {
        'length': length,
        'gc': gc_content,
        'tm': tm,
        'hairpin_tm': hairpin_tm,
        'dimer_tm': dimer_tm,
        'status': quality_status
    }

def check_pcr_products_advanced(f_hits, r_hits, max_product_size=5000):
    """
    升级版产物检查逻辑：
    不仅计算 F + R，还全面计算 F + F 以及 R + R 的自扩增条带。
    """
    products = []

    def find_bands(hits1, hits2, is_same_set=False):
        bands = []
        for i in range(len(hits1)):
            start_idx = i + 1 if is_same_set else 0
            for j in range(start_idx, len(hits2)):
                h1 = hits1[i]
                h2 = hits2[j]
                
                # 1. 必须在同一条染色体上
                if h1['chr'] != h2['chr']:
                    continue
                
                p_size = 0
                # 2. 方向必须面对面
                # 情况 A：h1在正链(+), h2在负链(-) -> h1在左，h2在右
                if h1['strand'] == '+' and h2['strand'] == '-':
                    if h1['start'] < h2['end']:
                        p_size = h2['end'] - h1['start'] + 1
                        
                # 情况 B：h1在负链(-), h2在正链(+) -> h2在左，h1在右
                elif h1['strand'] == '-' and h2['strand'] == '+':
                    if h2['start'] < h1['end']:
                        p_size = h1['end'] - h2['start'] + 1

                # 3. 长度在可扩增范围内
                if 0 < p_size <= max_product_size:
                    if is_same_set:
                        type_label = "引物自扩"
                    else:
                        type_label = "常规F+R"
                        
                    bands.append({
                        'chr': h1['chr'],
                        'size': p_size,
                        'type': type_label
                    })
        return bands

    # 交叉检索三种扩增可能性
    products.extend(find_bands(f_hits, r_hits, is_same_set=False))
    products.extend(find_bands(f_hits, f_hits, is_same_set=True))
    products.extend(find_bands(r_hits, r_hits, is_same_set=True))
                
    return products


def main():
    print("="*60)
    print("        引物特异性(Primer-BLAST) 终端快速验证工具        ")
    print("="*60)

    # 1. 基础交互
    print("\n[引物输入]")
    f_seq = input("请输入 Forward 引物序列 (5'->3'): ").strip().upper()
    r_seq = input("请输入 Reverse 引物序列 (5'->3'): ").strip().upper()
    
    pid = "Manual_Pair_1"
    max_size = 5000
    evalue = 100  # 放大 E-value，保证短序列匹配不漏点

    print("\n" + "-"*50)
    print(f"正在检查引物对: {pid} ...")
    print(f"Forward: {f_seq}")
    print(f"Reverse: {r_seq}")
    print("-"*50)
    
    # 写标准 FASTA 格式临时文件
    tmp_f_path = os.path.join(out_dir, "tmp_verify_f.fasta")
    tmp_r_path = os.path.join(out_dir, "tmp_verify_r.fasta")
    
    with open(tmp_f_path, "w") as f:
        f.write(f">f_{pid}\n{f_seq}\n")
    with open(tmp_r_path, "w") as f:
        f.write(f">r_{pid}\n{r_seq}\n")

    # 运行 BLAST
    # 💡 提示：如果发现染色体名字依然不准确，请进入你的 Primer_Blast_Fuction.py 的 blastn 函数中，
    # 将 chrom_name = alignment.title.split()[-1] 修改为 chrom_name = alignment.title.split()[0]
    f_hits = pb.blastn(tmp_f_path, db_address=db_address,
                       out_address1=os.path.join(out_dir, 'tmp_verify_F.xml'),
                       evalue=evalue, identity=len(f_seq)-2, task='blastn-short', dust='no')
                       
    r_hits = pb.blastn(tmp_r_path, db_address=db_address,
                       out_address1=os.path.join(out_dir, 'tmp_verify_R.xml'),
                       evalue=evalue, identity=len(r_seq)-2, task='blastn-short', dust='no')

    # 清理生成的临时 Fasta 文件，保持目录干净
    if os.path.exists(tmp_f_path): os.remove(tmp_f_path)
    if os.path.exists(tmp_r_path): os.remove(tmp_r_path)

    f_hit_count = len(f_hits)
    r_hit_count = len(r_hits)

    # 交叉计算条带
    predicted_bands = check_pcr_products_advanced(f_hits, r_hits, max_product_size=max_size)
    total_bands = len(predicted_bands)

    # 格式化输出条带详情：[大小, 染色体, 类型]
    band_details_list = [f"{b['size']}bp({b['chr']}, {b['type']})" for b in predicted_bands]
    band_details_str = ", ".join(band_details_list) if band_details_list else "无产物"

    # 2. 终端直接输出详细结果
    print(f"\n[BLAST 比对统计]")
    print(f"  -> Forward 引物在基因组共匹配到: {f_hit_count} 个位点")
    print(f"  -> Reverse 引物在基因组共匹配到: {r_hit_count} 个位点")
    print(f"  -> 预测潜在 PCR 条带总数: {total_bands} 个")

    print(f"\n[特异性评估结果]")
    if total_bands == 0:
        print("  -> ❌ 无法扩增 (No Product)")
    elif total_bands == 1:
        if predicted_bands[0]['type'] == "常规F+R":
            if f_hit_count > 4 or r_hit_count > 4:
                print(f"  -> ⚠️ 背景过高 (High Background) | 预测条带: {band_details_str}")
            else:
                print(f"  -> ✅ 特异性良好 (Specific) | 预测单条带: {band_details_str}")
        else:
            print(f"  -> ❌ 错误非特异扩增 | 预测条带: {band_details_str}")
    else:
        print(f"  -> ❌ 多条带脱靶 (Non-Specific)! 预测产物详情: {band_details_str}")
    print("="*60 + "\n")

    f_qual = evaluate_primer_quality(f_seq, "Forward")
    r_qual = evaluate_primer_quality(r_seq, "Reverse")
    
    print("\n[一、 引物热力学及物化质量评估]")
    print(f"🔹 Forward 引物报告:")
    print(f"  - 长度: {f_qual['length']} bp  |  GC含量: {f_qual['gc']:.1f}%  |  Tm值: {f_qual['tm']:.1f} ℃")
    print(f"  - 发夹折叠 (Hairpin Tm): {f_qual['hairpin_tm']:.1f} ℃  |  自身二聚体 (Dimer Tm): {f_qual['dimer_tm']:.1f} ℃")
    print(f"  - 质量评级: {f_qual['status']}")
    
    print(f"\n🔹 Reverse 引物报告:")
    print(f"  - 长度: {r_qual['length']} bp  |  GC含量: {r_qual['gc']:.1f}%  |  Tm值: {r_qual['tm']:.1f} ℃")
    print(f"  - 发夹折叠 (Hairpin Tm): {r_qual['hairpin_tm']:.1f} ℃  |  自身二聚体 (Dimer Tm): {r_qual['dimer_tm']:.1f} ℃")
    print(f"  - 质量评级: {r_qual['status']}")

if __name__ == '__main__':
    main()