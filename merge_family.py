import argparse
import vcf
from record import Record, PhaseSet, ChromosomoHaplotype
import sys
import os
import copy

def merge_set(c_info, c_phase_set, child=True):
    phase_set_keys = list(c_info.chromo_phase_set.keys())
    for phase_set_key in phase_set_keys:
        phase_set = c_info.chromo_phase_set[phase_set_key]
        if phase_set.origin == 2:
            phase_set.build_origin(child)
        if phase_set.origin != 2:
            if phase_set.origin == 1:
                phase_set.flip()
            c_info.connect_phase_set(c_phase_set, phase_set)


def child_haplotype(c_info, f_info = None,
                    m_info = None):
    c_phase_set, f_phase_set, m_phase_set = merge_unphased_snp(c_info, f_info, m_info)
    merge_set(c_info, c_phase_set)
    c_phase_set.finalize_phaseset_label()
    if f_info != None:
        merge_set(f_info, f_phase_set, False)
        f_phase_set.finalize_phaseset_label()
    if m_info != None:
        merge_set(m_info, m_phase_set, False)
        f_phase_set.finalize_phaseset_label()

def trio_heter_source(c_r, f_r, m_r):
    a = [f_r.hap0, f_r.hap1, m_r.hap0, m_r.hap1]
    if a.count(c_r.hap0) == 2 and a.count(c_r.hap1) == 2:
        return False
    elif a.count(c_r.hap0) == 1:
        if a.index(c_r.hap0) == 0 or a.index(c_r.hap0) == 1:
            c_r.origin = 0
        elif a.index(c_r.hap0) == 2 or a.index(c_r.hap0) == 3:
            c_r.origin = 1
        return True
    elif a.count(c_r.hap1) ==1:
        if a.index(c_r.hap1) == 0 or a.index(c_r.hap1) == 1:
            c_r.origin = 1
        elif a.index(c_r.hap1) == 2 or a.index(c_r.hap1) == 3:
            c_r.origin = 0
        return True
    return False


# parent hap origin
def p_hap_source(c_r, p_r):
    if c_r.hap0 == p_r.hap0:
        p_r.origin = 0
    if c_r.hap0 == p_r.hap1:
        p_r.origin = 1


# child hap origin
def c_hap_source(c_r, p_r, origin):
    if c_r.hap0 == p_r.hap0:
        c_r.origin = origin
    if c_r.hap1 == p_r.hap0:
        c_r.origin = abs(origin - 1)

def hete_to_ref_hom(rec):
    res = copy.deepcopy(rec)
    res.hap0 = res.ref
    res.hap1 = res.ref
    res.origin = -1
    return res

# def phase_with_large_seg(c_rec, p_rec, c_info, p_info:ChromosomoHaplotype):
#     # c_phase_set = c_info.chromo_phase_set[c_rec.ps]
#     p_phase_set = p_info.chromo_phase_set[p_rec.ps]
#     poses = list(p_phase_set.keys())
#     idx = poses.index(p_rec.pos)
#     n = 0;
#     p_hap0 = ''
#     c_hap0 = ''
#     # for i in poses[idx:min(idx+10, len(poses))]:
#     #     if i in c_info.chromo_record.keys() and c_info.chromo_record[i].phased():
#     #         p_tmp = c_info.chromo_record[i]
#     #         p_hap0= p_hap0+p_tmp.hap0

#     for i in poses[max(idx - 10, 0), min(idx+10, len(poses))]:
        

        
def merge_unphased_snp(c_info, f_info, m_info):
    # empty set , will filled by unphased snp
    c_phase_set = PhaseSet(1, 0)
    f_phase_set = PhaseSet(1, 0)
    m_phase_set = PhaseSet(1, 0)
    for pos in c_info.chromo_record.keys():
        c_rec = c_info.chromo_record[pos]
        f_rec = None
        m_rec = None
        if f_info != None and pos in f_info.chromo_record.keys():
            f_rec = f_info.chromo_record[pos]
        else:
            f_rec = hete_to_ref_hom(c_rec)
        if m_info != None and pos in m_info.chromo_record.keys():
            m_rec = m_info.chromo_record[pos]
        else:
            m_rec = hete_to_ref_hom(c_rec)
        if pos =='242212124' or pos ==242212124:
            print("xxx")
        # if child is homozygous try to phase parental.
        if not c_rec.is_heterozygous():
            if f_rec != None and f_rec.is_heterozygous():
                p_hap_source(c_rec, f_rec)
                if f_rec.origin != -1 and not f_rec.phased():
                    if f_rec.origin == 1:
                        f_rec.flip_num = 1
                    f_phase_set.insert_record(f_rec)
            if m_rec != None and m_rec.is_heterozygous():
                p_hap_source(c_rec, m_rec)
                if m_rec.origin != -1 and not m_rec.phased():
                    if m_rec.origin == 1:
                        m_rec.flip_num = 1
                    m_phase_set.insert_record(m_rec)
            continue

        # fixme: more flexible logistic
        # if child is heterozygous but parent none
        if f_rec == None and m_rec == None:
            continue
        # if child is heterozygous and parent both not none
        if f_rec != None and m_rec != None:
            f_h = f_rec.is_heterozygous()
            m_h = m_rec.is_heterozygous()
            # f and m both heter
            if f_h and m_h:
                trio_heter_source(c_rec,f_rec,m_rec)
            # f and m both home
            if not f_h and not m_h and f_rec.hap0 == m_rec.hap0:
                pass
            elif not f_h:
                c_hap_source(c_rec, f_rec, 0)
            elif not m_h:
                c_hap_source(c_rec, m_rec, 1)
        elif f_rec !=None:
            if not f_rec.is_heterozygous():
                c_hap_source(c_rec, f_rec, 0)
        elif m_rec != None:
            if not m_rec.is_heterozygous():
                c_hap_source(c_rec, m_rec, 1)
        # if child haplotype can infered from parent
        if c_rec.origin != -1 and not c_rec.phased():
            # we set all rec origin 0, means hap0 from father
            if c_rec.origin == 1:
                c_rec.flip_num = 1
            c_phase_set.insert_record(c_rec)
        # if c_rec.origin != -1 and c_rec.phased():
        #     if c_rec.origin == 1:
        #         c_rec.flip()
    c_phase_set.finalize_phaseset_label()
    f_phase_set.finalize_phaseset_label()
    m_phase_set.finalize_phaseset_label()
    return [c_phase_set, f_phase_set, m_phase_set]


def write_chromosome(in_vcf, out_vcf, chromo_haplotype, contig):
    rec = None
    for rec in in_vcf.fetch(contig):
        het = rec.samples[0].gt_type
        if het != 1:  # not het loci
            out_vcf.write_record(rec)
        else:
            record = chromo_haplotype.chromo_record[rec.POS]
            record.finalize_record(rec)
            out_vcf.write_record(rec)


def main():
    parser = argparse.ArgumentParser("merge family-based phased vcf")
    parser.add_argument(
        '-f', '--father', help='Father VCF file(s) in gz format, indexed, seperated by space', required=False)
    parser.add_argument(
        '-m', '--mother', help='Mother VCF file(s) in gz format, indexed, seperated by space', required=False)
    parser.add_argument(
        '-c', '--child', help='Child VCF file(s) in gz format, indexed, seperated by space', required=True)
    parser.add_argument('-o', '--out', help='Out dir', required=True)
    args = parser.parse_args()
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    if not args.father and not args.mother:
        print("error: At least one parent is present")
        exit(1)
    if not args.father:
        m_vcf = vcf.Reader(filename=args.mother)
        c_vcf = vcf.Reader(filename=args.child)
        c1 = args.child.split("/")[-1].replace(".gz", "")
        m1 = args.mother.split("/")[-1].replace(".gz", "")
        child_out_vcf = vcf.Writer(open(os.path.join(args.out, c1), 'w'), c_vcf)
        m_out_vcf = vcf.Writer(open(os.path.join(args.out, m1), 'w'), m_vcf)
        chromos = m_vcf.contigs.keys()
        for chromo in chromos:
            try:
                m_vcf.fetch(chromo)
            except:
                continue
            m_chromo_info = ChromosomoHaplotype(m_vcf, str(chromo))
            c_chromo_info = ChromosomoHaplotype(c_vcf, str(chromo))

            child_haplotype(c_chromo_info, m_chromo_info)
            write_chromosome(c_vcf, child_out_vcf, c_chromo_info, str(chromo))
            write_chromosome(m_vcf, m_out_vcf, m_chromo_info, str(chromo))

        child_out_vcf.close()
        m_out_vcf.close()
        return
    if not args.mother:
        f_vcf = vcf.Reader(filename=args.father)
        c_vcf = vcf.Reader(filename=args.child)
        c1 = args.child.split("/")[-1].replace(".gz", "")
        f1 = args.father.split("/")[-1].replace(".gz", "")
        child_out_vcf = vcf.Writer(open(os.path.join(args.out, c1), 'w'), c_vcf)
        f_out_vcf = vcf.Writer(open(os.path.join(args.out, f1), 'w'), f_vcf)
        chromos = f_vcf.contigs.keys()
        for chromo in chromos:
            try:
                f_vcf.fetch(chromo)
            except:
                continue
            f_chromo_info = ChromosomoHaplotype(f_vcf, str(chromo))
            c_chromo_info = ChromosomoHaplotype(c_vcf, str(chromo))

            child_haplotype(c_chromo_info, f_chromo_info)
            write_chromosome(c_vcf, child_out_vcf, c_chromo_info, str(chromo))
            write_chromosome(f_vcf, f_out_vcf, f_chromo_info, str(chromo))

        child_out_vcf.close()
        f_out_vcf.close()
        return
    else:
        f_vcf = vcf.Reader(filename=args.father)
        m_vcf = vcf.Reader(filename=args.mother)
        c_vcf = vcf.Reader(filename=args.child)

        c1 = args.child.split("/")[-1].replace(".gz", ".final")
        f1 = args.father.split("/")[-1].replace(".gz", ".final")
        m1 = args.mother.split("/")[-1].replace(".gz", ".final")
        child_out_vcf = vcf.Writer(open(os.path.join(args.out, c1), 'w'), c_vcf)
        f_out_vcf = vcf.Writer(open(os.path.join(args.out, f1), 'w'), f_vcf)
        m_out_vcf = vcf.Writer(open(os.path.join(args.out, m1), 'w'), m_vcf)

        chromos = f_vcf.contigs.keys()
        for chromo in chromos:
            try:
                f_vcf.fetch(chromo)
            except:
                continue
            f_chromo_info = ChromosomoHaplotype(f_vcf, str(chromo))
            m_chromo_info = ChromosomoHaplotype(m_vcf, str(chromo))
            c_chromo_info = ChromosomoHaplotype(c_vcf, str(chromo))

            child_haplotype(c_chromo_info, f_chromo_info, m_chromo_info)
            write_chromosome(c_vcf, child_out_vcf, c_chromo_info, str(chromo))
            write_chromosome(f_vcf, f_out_vcf, f_chromo_info, str(chromo))
            write_chromosome(m_vcf, m_out_vcf, m_chromo_info, str(chromo))

        child_out_vcf.close()
        f_out_vcf.close()
        m_out_vcf.close()


if __name__ == "__main__":
    main()
