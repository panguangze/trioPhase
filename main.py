import os
from subprocess import call
import argparse

def i_phase(spechap,extract,bam, vcf, out_dir, name):
    lst_out = os.path.join(out_dir, name+".lst")
    lst_sorted_out = os.path.join(out_dir, name+".sorted.lst")
    sort_cmd = "sort -n -k3 {} > {}".format(lst_out,lst_sorted_out)
    call(sort_cmd, shell=True)
    phased_vcf = os.path.join(out_dir, name+".phased.vcf")
    e_cmd = "{} --vcf {} --bam {} -o {}".format(extract, vcf, bam, lst_sorted_out)
    s_cmd = "{} -f {} -v {} -o {}".format(spechap, lst_out, vcf, phased_vcf)
    call(e_cmd, shell=True)
    call(s_cmd, shell=True)
    return phased_vcf
def e_lst(extract, bams, vcf, out_dir):
    lst_a = ""
    lst_l = os.path.join(out_dir,"all.lst")
    lst_l_sorted = os.path.join(out_dir,"all.sorted.lst")
    for i in range(len(bams)):
        b = bams[i]
        lst = os.path.join(out_dir, "tmp"+str(i)+".lst")
        e_cmd = "{} --vcf {} --bam {} -o {}".format(extract, vcf, b, lst)
        lst_a = lst_a+" "+lst
        call(e_cmd, shell=True)
    cat_cmd = "cat {} > {}".format(lst_a, lst_l)
    sort_cmd = "sort -n -k3 {} > {}".format(lst_l,lst_l_sorted)
    call(cat_cmd, shell=True)
    call(sort_cmd, shell=True)
    return lst_l_sorted
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
        '--spechap', help='spechap paht', required=True)
    parser.add_argument(
        '--extractHairs', help='extractHairs path', required=True)
    parser.add_argument('-o', '--out_dir', help='Out dir', required=True)    
    
    args = parser.parse_args()
    print("individual phase")
    c_phased_v1 = i_phase(args.spechap, args.extractHairs, args.child_b, args.child_v, args.out_dir, "child")
    if args.mother_v and args.mother_b:
        m_phased_v1 = i_phase(args.spechap, args.extractHairs, args.mother_b, args.mother_v, args.out_dir, "mother")
    if args.father_v and args.father_b:
        f_phased_v1 = i_phase(args.spechap, args.extractHairs, args.father_b, args.father_v, args.out_dir, "father")

    print("Raw phase with only vcf")
    raw_cmd = ""
    if not args.mother_v:
        raw_cmd = "python merge_family.py -f {} -c {} -o {}".format(f_phased_v1, c_phased_v1, args.out_dir)
    if not args.father_v:
        raw_cmd = "python merge_family.py -m {} -c {} -o {}".format(m_phased_v1, c_phased_v1, args.out_dir)
    else:
        raw_cmd = "python merge_family.py -m {} -f {} -c {} -o {}".format(m_phased_v1, f_phased_v1, c_phased_v1, args.out_dir)
    if call(raw_cmd, shell=True):
        print(raw_cmd,"running error")
        exit

if __name__ == "__main__":
    main()
