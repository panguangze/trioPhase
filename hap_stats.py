import argparse
import vcf
from record import Record, PhaseSet, ChromosomoHaplotype
from stats import PhaseSetStats, HapStats


def get_phase_set_stats(phase_set:PhaseSet):
    record: Record
    record_count = 0
    total_record = len(phase_set.records_idx)
    last_record_pos = 0
    last_record_idx = 0
    first_record_idx = 0
    prev_record_idx = 0
    prev_record_pos = 0
    curr_record_idx = 0
    curr_record_pos = 0
    positions = []
    for record_pos in phase_set.records.keys():
        record = phase_set.records[record_pos]
        record_count += 1
        curr_record_idx = record.idx
        curr_record_pos = record.pos
        if record_count == total_record:
            last_record_idx = record.idx
            last_record_pos = record.pos
            positions.append(str(last_record_pos))
        elif record_count == 1:
            first_record_idx = record.idx
            prev_record_idx = first_record_idx
            prev_record_pos = record.pos
            positions.append(str(prev_record_pos))
        else:
            if curr_record_idx != prev_record_idx + 1:
                positions.append(str(prev_record_pos))
                positions.append(str(curr_record_pos))
        prev_record_idx = curr_record_idx
        prev_record_pos = curr_record_pos
        
    S50 = total_record
    N50 = last_record_pos - phase_set.starting_pos
    spaned_record = last_record_idx - first_record_idx + 1
    AN50 = N50/spaned_record * S50
    sep = ";"
    poss = sep.join(positions)
    return AN50, S50, N50, spaned_record, poss


def get_haplotype_stats_chromo(in_chromo:ChromosomoHaplotype, out, contig, chromo_snp_count, chromo_span):
    phase_set : PhaseSet
    hap_stats = HapStats(chromo_snp_count, chromo_span)
    index = 0
    for phase_set in in_chromo.chromo_phase_set.values():
        AN50, S50, N50, spanned_snp, position_str = get_phase_set_stats(phase_set)
        phase_set_stats = PhaseSetStats(0, 0, S50, N50, AN50, spanned_snp)
        hap_stats.insert_phase_set_stats(0, phase_set_stats)
        index += 1
        out.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (contig, phase_set_stats.get_AN50(), phase_set_stats.get_N50(), phase_set_stats.get_phased_snp(), spanned_snp, position_str))
    out.write("%s\t%d\t%d\t%d\t%d\tNA\n" % (contig + "_total", hap_stats.get_AN50(), hap_stats.get_N50(), hap_stats.get_total_phased(), hap_stats.get_total_spanned()))
    return hap_stats


def count_contig_stats(in_vcf:vcf.Reader, contig: str):
    count = 0
    starting_loci = 0
    end_loci = 0 
    for rec in in_vcf.fetch(contig):
        count += 1
        end_loci = rec.POS
        if starting_loci == 0:
            starting_loci = end_loci
    return  count, end_loci - starting_loci


def get_haplotype_stats(in_vcf:vcf.Reader, out):
    contigs = in_vcf.contigs.keys()
    hap_stats = HapStats()
    for contig in contigs:
        try: 
            SNV_count, contig_span = count_contig_stats(in_vcf, contig)
            in_chromo = ChromosomoHaplotype(in_vcf, contig)
            chromo_hap_stats = get_haplotype_stats_chromo(in_chromo, out, contig, SNV_count, contig_span)
            hap_stats.insert_hap_stats(chromo_hap_stats)
        except:
            continue
    out.write("%s\t%d\t%d\t%d\t%d\tNA\n" % ("total", hap_stats.get_AN50(), hap_stats.get_N50(), hap_stats.get_total_phased(), hap_stats.get_total_spanned()))

def main():
    parser = argparse.ArgumentParser('hap_completeness.py')
    parser.add_argument('-v', '--vcf', help='input vcf, indexed', required=True)
    parser.add_argument('-o', '--out', help='output stats', required=True)
    options = parser.parse_args()

    inf = vcf.Reader(filename=options.vcf)
    out = open(options.out, "w")
    out.write("Chromosome\tAN50\tN50\tphased_snp\ttotal_snp\tbreak_point\n")
    get_haplotype_stats(inf, out)
    out.close()
    return


if __name__ == '__main__':
    main()

