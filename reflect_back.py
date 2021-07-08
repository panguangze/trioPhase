import argparse
import os
import vcf
def main():
    parser = argparse.ArgumentParser("reflect")
    parser.add_argument(
        '--hete', help='hete vcf', required=True)
    parser.add_argument(
        '--orig', help='original vcf', required=True)
    parser.add_argument(
        '--out', help='out vcf', required=True)
    args = parser.parse_args()
    c_vcf = vcf.Reader(filename=args.hete)
    # orig_vcf = vcf.Reader(filename=args.orig)
    orig_vcf = open(args.orig)

    out_vcf = vcf.Writer(open(args.out, 'w'), c_vcf)
    out_vcf.close()
    out_vcf = open(args.out,"a")
    #chromos=["1"]
    chromos = c_vcf.contigs.keys()
    hete_infos = {}
    for chromo in chromos:
        tmp= {}
        try:
            c_vcf.fetch(chromo)
        except:
            continue
        for rec in c_vcf.fetch(chromo):
            try:
                tmp[rec.POS]=[str(rec.samples[0]['GT']),str(rec.samples[0]['PS'])]
            except:
                tmp[rec.POS]=[str(rec.samples[0]['GT']),"0"]
        hete_infos[chromo] = tmp

    for line in orig_vcf.readlines():
        if line.startswith("#"):
            continue
        a = line.split("\t")
        chromo = a[0]
        pos = a[1]
        snp_format = a[-2]
        if "PS" not in snp_format:
            a[-2] = a[-2]+":PS"
        vs = a[-1].split(":")
        if int(pos) in hete_infos[chromo].keys():
            target = hete_infos[chromo][int(pos)]
            if len(vs) != 3:
                vs.append(target[1])
            else:
                vs[-1] = target[1]
            vs[0] = target[0]
        else:
            if len(vs) != 3:
                vs.append('0')
            else:
                vs[-1] = '0'
        a[-1] = ':'.join(vs)
        res = '\t'.join(a)
        res = res.replace('\r', '').replace('\n', '')
        out_vcf.write(res+"\n")

if __name__ == "__main__":
    main() 
