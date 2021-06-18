import os
import argparse
def execute_cmd(cmd):
    print("Info: ", cmd)
    os.system(cmd)
def bgzip_and_index(vcf, bgzip, tabix):
    bgzip_cmd = "{} {}".format(bgzip, vcf)
    tabix_cmd = "{} {}".format(tabix, vcf+".gz")
    execute_cmd(bgzip_cmd)
    execute_cmd(tabix_cmd)
    return vcf+".gz"
def bgunzip(vcf):
    if "gz" in vcf:
        unzip_cmd = "gunzip {}".format(vcf)
        execute_cmd(unzip_cmd)
        return vcf.replace(".gz","")
    else:
        return vcf
def sc_file(d,df):
    # f_name = f.split("/")[-1]
    ds = os.path.join(d, df)
    # mv_cmd = "mv {} {}".format(f,ds)
    # execute_cmd(mv_cmd)
    return ds
def i_phase(spechap,extract,bgzip,tabix,bam, vcf, out_dir, name, ref):
    vcf = bgunzip(vcf)
    lst_out = os.path.join(out_dir, name+".lst")
    lst_sorted_out = os.path.join(out_dir, name+".sorted.lst")
    phased_vcf = os.path.join(out_dir, name+".phased.vcf")
    sort_cmd = "sort -n -k3 {} > {}".format(lst_out,lst_sorted_out)
    e_cmd = "{} --vcf {} --bam {} -o {} --ref {}".format(extract, vcf, bam, lst_out, ref)
    execute_cmd(e_cmd)
    execute_cmd(sort_cmd)
    vcf = bgzip_and_index(vcf, bgzip, tabix)
    s_cmd = "{} -f {} -v {} -o {}".format(spechap, lst_sorted_out, vcf, phased_vcf)
    execute_cmd(s_cmd)
    return phased_vcf
def e_lst(extract, bams, vcf, out_dir, name, ref):
    lst_a = ""
    lst_l = os.path.join(out_dir,name+".all.lst")
    lst_l_sorted = os.path.join(out_dir,name+".all.sorted.lst")
    for i in range(len(bams)):
        b = bams[i]
        lst = os.path.join(out_dir, name+".tmp"+str(i)+".lst")
        e_cmd = "{} --vcf {} --bam {} -o {} --ref {}".format(extract, vcf, b, lst,ref)
        lst_a = lst_a+" "+lst
        execute_cmd(e_cmd)
    cat_cmd = "cat {} > {}".format(lst_a, lst_l)
    sort_cmd = "sort -n -k3 {} > {}".format(lst_l,lst_l_sorted)
    execute_cmd(cat_cmd)
    execute_cmd(sort_cmd)
    return lst_l_sorted
def phase_with_lst(spechap, lst, vcf, out_file, bgzip, tabix):
    vcf = bgzip_and_index(vcf,bgzip,tabix)
    s_cmd = "{} -f {} -v {} -o {} --keep_phasing_info".format(spechap, lst, vcf, out_file)
    execute_cmd(s_cmd)
    bgzip_and_index(out_file,bgzip,tabix)
    return out_file+".gz"
def main():
    parser = argparse.ArgumentParser("trio phase")
    parser.add_argument(
        '--father_v', help='Father VCF file in gz format, indexed, seperated by space', required=False)
    parser.add_argument(
        '--mother_v', help='Mother VCF file in gz format, indexed, seperated by space', required=False)
    parser.add_argument(
        '--child_v', help='Child VCF file in gz format, indexed, seperated by space', required=True)
    parser.add_argument(
        '--father_b', help='Father bam file indexed', required=False)
    parser.add_argument(
        '--mother_b', help='Mother bam file indexed', required=False)
    parser.add_argument(
        '--child_b', help='Child bam file indexed', required=True)
    parser.add_argument(
        '--spechap', help='spechap path', required=True)
    parser.add_argument(
        '--extractHairs', help='extractHairs path', required=True)
    parser.add_argument(
        '--bgzip', help='bgzip path', required=True)
    parser.add_argument(
        '--tabix', help='tabix path', required=True)
    parser.add_argument('-o', '--out_dir', help='Out dir', required=True)
    parser.add_argument('--ref', help='reference', required=True)
    parser.add_argument(
        '--step', help='which step', required=False)
    args = parser.parse_args()

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    if args.step == "2":
        print("skip individual phase")
        m_phased_v1= bgunzip(args.mother_v)
        f_phased_v1= bgunzip(args.father_v)
        c_phased_v1= bgunzip(args.child_v)
        c_phased_v2 = sc_file(args.out_dir,"child.vcf")
        m_phased_v2 = sc_file(args.out_dir,"mother.vcf")
        f_phased_v2 = sc_file(args.out_dir,"father.vcf")
    else:
        print("individual phase")
        m_phased_v1=""
        f_phased_v1=""
        c_phased_v1 = i_phase(args.spechap, args.extractHairs, args.bgzip, args.tabix, args.child_b, args.child_v, args.out_dir, "child", args.ref)
        if args.mother_v and args.mother_b:
            m_phased_v1 = i_phase(args.spechap, args.extractHairs,args.bgzip, args.tabix, args.mother_b, args.mother_v, args.out_dir, "mother", args.ref)
        if args.father_v and args.father_b:
            f_phased_v1 = i_phase(args.spechap, args.extractHairs,args.bgzip, args.tabix, args.father_b, args.father_v, args.out_dir, "father", args.ref)
        c_phased_v2 = c_phased_v1+".trio.vcf"
        m_phased_v2 = m_phased_v1+".trio.vcf"
        f_phased_v2 = f_phased_v1+".trio.vcf"
    # print("Raw phase with only vcf")
    # raw_cmd = ""
    # if not args.mother_v:
    #     raw_cmd = "python merge_family.py -f {} -c {} -o {}".format(f_phased_v1, c_phased_v1, args.out_dir)
    # if not args.father_v:
    #     raw_cmd = "python merge_family.py -m {} -c {} -o {}".format(m_phased_v1, c_phased_v1, args.out_dir)
    # else:
    #     raw_cmd = "python merge_family.py -m {} -f {} -c {} -o {}".format(m_phased_v1, f_phased_v1, c_phased_v1, args.out_dir)
    # if execute_cmd(raw_cmd):
    #     print(raw_cmd,"running error")
    #     exit


    print("phasing child ...")
    bams = []
    if args.mother_b:
        bams.append(args.mother_b)
    if args.father_b:
        bams.append(args.father_b)
    c_lst = e_lst(args.extractHairs, bams, c_phased_v1, args.out_dir, "child", args.ref)
    c_out = phase_with_lst(args.spechap, c_lst, c_phased_v1, c_phased_v2, args.bgzip, args.tabix)

    print("phasing parent...")
    if args.mother_b and args.mother_v:
        bams = [args.child_b]
        m_lst = e_lst(args.extractHairs, bams, m_phased_v1, args.out_dir, "mother", args.ref)
        m_out = phase_with_lst(args.spechap, m_lst, m_phased_v1, m_phased_v2, args.bgzip, args.tabix)
    if args.father_b and args.father_v:
        bams = [args.child_b]
        f_lst = e_lst(args.extractHairs, bams, f_phased_v1, args.out_dir, "father", args.ref)
        f_out = phase_with_lst(args.spechap, f_lst, f_phased_v1, f_phased_v2, args.bgzip, args.tabix)
    
    print("Phase with only vcf")
    raw_cmd = ""
    if not args.mother_v:
        raw_cmd = "python merge_family.py -f {} -c {} -o {}".format(f_phased_v2+".gz", c_phased_v2+".gz", args.out_dir)
    if not args.father_v:
        raw_cmd = "python merge_family.py -m {} -c {} -o {}".format(m_phased_v2+".gz", c_phased_v2+".gz", args.out_dir)
    else:
        raw_cmd = "python merge_family.py -m {} -f {} -c {} -o {}".format(m_phased_v2+".gz", f_phased_v2+".gz", c_phased_v2+".gz", args.out_dir)
    if execute_cmd(raw_cmd):
        print(raw_cmd,"running error")
        exit

if __name__ == "__main__":
    main()
