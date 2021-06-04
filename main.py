import os
import argparse
def bgzip_and_index(vcf, bgzip, tabix):
    bgzip_cmd = "{} {}".format(bgzip, vcf)
    tabix_cmd = "{} {}".format(tabix, vcf+".gz")
    os.system(bgzip_cmd)
    os.system(tabix_cmd)
def bgunzip(vcf):
    unzip_cmd = "gunzip {}".format(vcf)
    os.system(unzip_cmd)
    return vcf.replace(".gz","")
def i_phase(spechap,extract,bgzip,tabix,bam, vcf, out_dir, name):
    vcf = bgunzip(vcf)
    lst_out = os.path.join(out_dir, name+".lst")
    lst_sorted_out = os.path.join(out_dir, name+".sorted.lst")
    phased_vcf = os.path.join(out_dir, name+".phased.vcf")
    sort_cmd = "sort -n -k3 {} > {}".format(lst_out,lst_sorted_out)
    e_cmd = "{} --vcf {} --bam {} -o {}".format(extract, vcf, bam, lst_out)
    s_cmd = "{} -f {} -v {} -o {}".format(spechap, lst_sorted_out, vcf, phased_vcf)
    os.system(e_cmd)
    os.system(sort_cmd)
    bgzip_and_index(vcf)
    os.system(s_cmd)
    return phased_vcf
def e_lst(extract, bams, vcf, out_dir, name):
    lst_a = ""
    lst_l = os.path.join(out_dir,name+".all.lst")
    lst_l_sorted = os.path.join(out_dir,name+".all.sorted.lst")
    for i in range(len(bams)):
        b = bams[i]
        lst = os.path.join(out_dir, name+".tmp"+str(i)+".lst")
        e_cmd = "{} --vcf {} --bam {} -o {}".format(extract, vcf, b, lst)
        lst_a = lst_a+" "+lst
        os.system(e_cmd)
    cat_cmd = "cat {} > {}".format(lst_a, lst_l)
    sort_cmd = "sort -n -k3 {} > {}".format(lst_l,lst_l_sorted)
    print(cat_cmd)
    print(sort_cmd)
    os.system(cat_cmd)
    os.system(sort_cmd)
    return lst_l_sorted
def phase_with_lst(spechap, lst, vcf, out_file, bgzip, tabix):
    bgzip_and_index(vcf,bgzip,tabix)
    s_cmd = "{} -f {} -v {} -o {} --keep_phasing_info".format(spechap, lst, vcf+".gz", out_file)
    os.system(s_cmd)
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
    args = parser.parse_args()

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    print("individual phase")
    m_phased_v1=""
    f_phased_v1=""
    c_phased_v1 = i_phase(args.spechap, args.extractHairs, args.bgzip, args.tabix, args.child_b, args.child_v, args.out_dir, "child")
    if args.mother_v and args.mother_b:
        m_phased_v1 = i_phase(args.spechap, args.extractHairs,args.bgzip, args.tabix, args.mother_b, args.mother_v, args.out_dir, "mother")
    if args.father_v and args.father_b:
        f_phased_v1 = i_phase(args.spechap, args.extractHairs,args.bgzip, args.tabix, args.father_b, args.father_v, args.out_dir, "father")

    # print("Raw phase with only vcf")
    # raw_cmd = ""
    # if not args.mother_v:
    #     raw_cmd = "python merge_family.py -f {} -c {} -o {}".format(f_phased_v1, c_phased_v1, args.out_dir)
    # if not args.father_v:
    #     raw_cmd = "python merge_family.py -m {} -c {} -o {}".format(m_phased_v1, c_phased_v1, args.out_dir)
    # else:
    #     raw_cmd = "python merge_family.py -m {} -f {} -c {} -o {}".format(m_phased_v1, f_phased_v1, c_phased_v1, args.out_dir)
    # if os.system(raw_cmd):
    #     print(raw_cmd,"running error")
    #     exit


    c_phased_v2 = c_phased_v1+".trio.vcf"
    m_phased_v2 = m_phased_v1+".trio.vcf"
    f_phased_v2 = f_phased_v1+".trio.vcf"

    print("phasing child ...")
    bams = []
    if args.mother_b:
        bams.append(args.mother_b)
    if args.father_b:
        bams.append(args.father_b)
    c_lst = e_lst(args.extractHairs, bams, c_phased_v1, args.out_dir, "child")
    phase_with_lst(args.spechap, c_lst, c_phased_v1, c_phased_v2, args.bgzip, args.tabix)

    print("phasing parent...")
    if args.mother_b and args.mother_v:
        bams = [args.child_b]
        m_lst = e_lst(args.extractHairs, bams, m_phased_v1, args.out_dir, "mother")
        phase_with_lst(args.spechap, m_lst, m_phased_v1, m_phased_v2, args.bgzip, args.tabix)
    if args.father_b and args.father_v:
        bams = [args.child_b]
        f_lst = e_lst(args.extractHairs, bams, f_phased_v1, args.out_dir, "father")
        phase_with_lst(args.spechap, f_lst, f_phased_v1, f_phased_v2, args.bgzip, args.tabix)
    
    print("Phase with only vcf")
    raw_cmd = ""
    if not args.mother_v:
        raw_cmd = "python merge_family.py -f {} -c {} -o {}".format(f_phased_v2+".gz", c_phased_v2+".gz", args.out_dir)
    if not args.father_v:
        raw_cmd = "python merge_family.py -m {} -c {} -o {}".format(m_phased_v2+".gz", c_phased_v2+".gz", args.out_dir)
    else:
        raw_cmd = "python merge_family.py -m {} -f {} -c {} -o {}".format(m_phased_v2+".gz", f_phased_v2+".gz", c_phased_v2+".gz", args.out_dir)
    if os.system(raw_cmd):
        print(raw_cmd,"running error")
        exit

if __name__ == "__main__":
    main()
