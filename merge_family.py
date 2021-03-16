import argparse
import vcf
from record import Record, PhaseSet, ChromosomoHaplotype
import sys
import os
from graph import Node, Graph
import math


def merge_chromo_haplotype(chromo_info1: ChromosomoHaplotype, chromo_info2: ChromosomoHaplotype):
    chromo_graph = create_phase_set_graph(chromo_info1, chromo_info2)
    # print(len(chromo_graph.connected_components))
    merge_graph(chromo_graph, chromo_info1, chromo_info2)
    print("merged")
    chromo_info1.finalize_chromosome_haplotype()


def child_haplotype(c_info: ChromosomoHaplotype, f_info: ChromosomoHaplotype, m_info: ChromosomoHaplotype):
    merge_phased_set = merge_unphased_snp(c_info, f_info, m_info)
    # print(len(merge_phased_set.records_idx))
    # phase_set_keys = list(c_info.chromo_phase_set.keys())
    # start_phase_set = c_info.chromo_phase_set[phase_set_keys[0]]
    # for phase_set_key in phase_set_keys:
    #     if phase_set_key == phase_set_keys[0]:
    #         continue
    #     phase_set = c_info.chromo_phase_set[phase_set_key]
    #     if phase_set.origin == 2:
    #         phase_set.build_origin()
    #     if phase_set.origin != 2:
    #         c_info.connect_phase_set(start_phase_set, phase_set)

def p_hap_source(c_r: Record, p_r: Record):
    if c_r.hap1 == p_r.hap1:
        p_r.origin = 0
    if c_r.hap1 == p_r.hap2:
        p_r.origin = 1

def c_hap_source(c_r: Record, p_r: Record, origin):
    if c_r.hap1 == p_r.hap1:
        c_r.origin = origin
    if c_r.hap2 == p_r.hap1:
        c_r.origin = math.abs(origin - 1)

def merge_unphased_snp(c_info: ChromosomoHaplotype, f_info: ChromosomoHaplotype, m_info: ChromosomoHaplotype):
    c_phase_set = PhaseSet(1)
    f_phase_set = PhaseSet(1)
    m_phase_set = PhaseSet(1)
    for pos in c_info.chromo_record.keys():
        c_rec = c_info.chromo_record[pos]
        f_rec: Record = None
        m_rec: Record = None
        if pos in f_info.chromo_record.keys():
            f_rec = f_info.chromo_record[pos]
        if pos in m_info.chromo_record.keys():
            m_rec = m_info.chromo_record[pos]
        
        # if child is homozygous try to phase parental.
        if not c_rec.is_heterozygous():
            if f_rec != None and f_rec.is_heterozygous():
                p_hap_source(c_rec, f_rec)
                if f_rec.origin != -1 and f_rec.origin != 2:
                    f_phase_set.insert_record(f_rec)
            if m_rec != None and m_rec.is_heterozygous():
                p_hap_source(c_rec, m_rec)
                if m_rec.origin != -1 and m_rec.origin != 2:
                    m_phase_set.insert_record(f_rec)
            continue
        
        # if child is heterozygous but parent none
        if f_rec == None and m_rec == None:
            continue
        # if child is heterozygous and parent both not none
        if f_rec != None and m_rec != None:
            f_h = f_rec.is_heterozygous()
            m_h = m_rec.is_heterozygous()
            # f and m both heter
            if f_h and m_h:
                continue
            # f and m both home
            if not f_h and not m_h and f_rec.hap1 == m_rec.hap1:
                continue
            if not f_h:
                c_hap_source(c_rec, f_rec, 0)
            if not m_h:
                c_hap_source(c_rec, m_rec, 1)
        elif f_rec:
            if f_rec.is_heterozygous():
                continue
            c_hap_source(c_rec, f_rec, 0)
        elif m_rec:
            m_h = m_rec.is_heterozygous()
            if m_rec.is_heterozygous():
                continue
            c_hap_source(c_rec, f_rec, 0)
        # if child haplotype can infered from parent
        if c_rec.origin != -1 and c_rec.origin != 2:
            c_phase_set.insert_record(c_rec)
    c_phase_set.set_origin(0)
    c_phase_set.finalize_phaseset_label()
    return c_phase_set


def write_chromosome(in_vcf: vcf.Reader, out_vcf: vcf.Writer, chromo_haplotype: ChromosomoHaplotype, contig: str):
    rec: vcf.model._Record
    for rec in in_vcf.fetch(contig):
        het = rec.samples[0].gt_type
        if het != 1:  # not het loci
            out_vcf.write_record(rec)
        else:
            record = chromo_haplotype.chromo_record[rec.POS]
            record.finalize_record(rec)
            out_vcf.write_record(rec)


def create_phase_set_graph(chromo_info1: ChromosomoHaplotype, chromo_info2: ChromosomoHaplotype):
    chromo_graph = Graph()
    for phase_set in chromo_info1.chromo_phase_set.values():
        chromo_graph.insert_node(Node(phase_set.starting_pos, 0))  # primary node
    for phase_set in chromo_info2.chromo_phase_set.values():
        chromo_graph.insert_node(Node(phase_set.starting_pos, 1))  # secondary Node

    phase_set: PhaseSet
    for phase_set in chromo_info1.chromo_phase_set.values():
        met_phase_set = set()
        for record_pos in phase_set.records_idx:
            if record_pos not in chromo_info2.chromo_record2phaseset_map.keys():  # not in secondary phase set
                continue
            ps_secondary = chromo_info2.chromo_record2phaseset_map[record_pos]
            if ps_secondary in met_phase_set:  # connection already found
                continue
            chromo_graph.add_edge_with_ps(phase_set.starting_pos, ps_secondary)
        met_phase_set.clear()

    chromo_graph.load_connected_components()
    return chromo_graph


def merge_graph(chromo_graph: Graph, chromo_info1: ChromosomoHaplotype, chromo_info2: ChromosomoHaplotype):
    """merge_chromo_haplotype [merge chromosome level haplotype based on the created graph structure]
        note that connected component are listed in such order: s - f -s or f-s-f
    Args:
        chromo_haplotype (ChromosomoHaplotype): [secondary chromosome level hap]
    """
    for connected_phase_sets in chromo_graph.connected_components:
        n_nodes = len(connected_phase_sets)
        if n_nodes == 1:
            continue
        else:
            start_node_id = chromo_graph.get_closeset_primary_node(connected_phase_sets[0])
            start_phase_set = chromo_info1.chromo_phase_set[chromo_graph.get_node(start_node_id).ps]
            visited = dict()
            chromo1_list = list()
            for id in connected_phase_sets:
                visited[id] = False
            queue = []
            queue.append(start_node_id)
            visited[start_node_id] = True

            while queue:
                s = queue.pop(0)
                if s != start_node_id:
                    node = chromo_graph.get_node(s)
                    c_phase_set: PhaseSet
                    node_phase_set: PhaseSet
                    if not node.is_primary():
                        node_phase_set = chromo_info2.chromo_phase_set[node.ps]
                        # add phased recodes that appears in chromo2 but not phased in chromo1
                        for r_idx in node_phase_set.records_idx:
                            if r_idx not in chromo_info1.chromo_record2phaseset_map.keys() and r_idx in chromo_info1.chromo_record.keys():
                                rec_1 = chromo_info1.chromo_record[r_idx]
                                if rec_1.hap == rec_1.hap1:
                                    continue
                                rec_2 = chromo_info2.chromo_record[r_idx]
                                rec_1.ps = rec_1.pos
                                if rec_1.hap != rec_2.hap:
                                    rec_1.flip()
                                start_phase_set.insert_record(rec_1)
                        # f_phase_set = chromo_info1.chromo_phase_set[node.ps]
                        # self.connect_phase_set(start_phase_set, r_phase_set)
                        for i in chromo_graph.adj_list[s]:
                            if i != start_node_id:
                                c_node = chromo_graph.get_node(i)
                                # all c_phase_set is primary node
                                c_phase_set = chromo_info1.chromo_phase_set[c_node.ps]
                                need_flip = node_phase_set.intersect_phase_set(c_phase_set)
                                if need_flip:
                                    if i in chromo1_list:
                                        for j in chromo1_list:
                                            j_phase_set = chromo_info1.chromo_phase_set[chromo_graph.get_node(j).ps]
                                            j_phase_set.flip()
                                    else:
                                        chromo_info1.chromo_phase_set[chromo_graph.get_node(i).ps].flip()
                                        chromo1_list.append(i)
                                else:
                                    if i not in chromo1_list:
                                        chromo1_list.append(i)
                            if visited[i] == False:
                                queue.append(i)
                                visited[i] = True
                    else:
                        node_phase_set = chromo_info1.chromo_phase_set[node.ps]
                        for i in chromo_graph.adj_list[s]:
                            sf_phase_set = chromo_info2.chromo_phase_set[chromo_graph.get_node(i).ps]
                            need_flip = node_phase_set.intersect_phase_set(sf_phase_set)
                            if need_flip:
                                sf_phase_set.flip()
                            if visited[i] == False:
                                queue.append(i)
                                visited[i] = True
                else:
                    for i in chromo_graph.adj_list[s]:
                        sf_phase_set = chromo_info2.chromo_phase_set[chromo_graph.get_node(i).ps]
                        need_flip = start_phase_set.intersect_phase_set(sf_phase_set)
                        if need_flip:
                            sf_phase_set.flip()
                        if visited[i] == False:
                            queue.append(i)
                            visited[i] = True
            # print(chromo1_list)
            # connect chromo1_list
            for i in chromo1_list:
                ps = chromo_graph.get_node(i).ps
                j_phase_set = chromo_info1.chromo_phase_set[ps]
                chromo_info1.connect_phase_set(start_phase_set, j_phase_set)


def main():
    parser = argparse.ArgumentParser("merge family-based phased vcf")
    parser.add_argument(
        '-f', '--father', help='Father VCF file(s) in gz format, indexed, seperated by space', required=True)
    parser.add_argument(
        '-m', '--mother', help='Mother VCF file(s) in gz format, indexed, seperated by space', required=True)
    parser.add_argument(
        '-c', '--child', help='Child VCF file(s) in gz format, indexed, seperated by space', required=True)
    parser.add_argument('-o', '--out', help='Out dir', required=True)
    args = parser.parse_args()
    f_vcf = vcf.Reader(filename=args.father)
    m_vcf = vcf.Reader(filename=args.mother)
    c_vcf = vcf.Reader(filename=args.child)
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    child_out_vcf = vcf.Writer(open(os.path.join(args.out, 'child.vcf'), 'w'), c_vcf)
    f_out_vcf = vcf.Writer(open(os.path.join(args.out, 'father.vcf'), 'w'), f_vcf)
    m_out_vcf = vcf.Writer(open(os.path.join(args.out, 'mother.vcf'), 'w'), m_vcf)

    # chromos = f_vcf.contigs.keys()
    chromos = ['chr1']
    for chromo in chromos:
        try:
            f_vcf.fetch(chromo)
        except:
            continue
        [f_chromo_info, m_chromo_info, c_chromo_info] = [ChromosomoHaplotype(f, str(chromo)) for f in
                                                         [f_vcf, m_vcf, c_vcf]]
        # merge_chromo_haplotype(c_chromo_info, f_chromo_info)
        # merge_chromo_haplotype(c_chromo_info, m_chromo_info)
        child_haplotype(c_chromo_info, f_chromo_info, m_chromo_info)
        # merge_chromo_haplotype(c_chromo_info, f_chromo_info)
        # merge_chromo_haplotype(c_chromo_info, m_chromo_info)
        write_chromosome(c_vcf, child_out_vcf, c_chromo_info, str(chromo))
    child_out_vcf.close()
    f_out_vcf.close()
    m_out_vcf.close()


if __name__ == "__main__":
    main()
