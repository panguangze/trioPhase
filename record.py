import vcf
import graph
import numpy as np

class Record:
    def set_origin(self,o):
        self.origin = o

    def __init__(self, rec: vcf.model._Record, PS: int, idx:int):
        allels=[rec.REF, rec.ALT]
        self.pos = rec.POS
        gt_str = rec.samples[0]['GT']
        gt0 = int(gt_str[0])
        gt1 = int(gt_str[2])
        if gt0 == 0:
            self.hap1 = rec.REF
        else:
            self.hap1 = rec.ALT[gt0 - 1]
        if gt1 == 0:
            self.hap2 = rec.REF
        else:
            self.hap2 = rec.ALT[gt1 - 1]
        self.ps = PS
        self.idx = idx
    
    def is_heterozygous(self):
        return self.hap1 != self.hap2

    def flip(self):
        t = hap2
        self.hap2 = self.hap1
        self.hap1 = t
    
    def phased(self):
        return self.ps != 0
    
    def finalize_record(self, rec: vcf.model._Record):
        if not self.phased():
            return

        gt_str = rec.samples[0]['GT']
        gt0 = int(gt_str[0])
        gt1 = int(gt_str[2])
        hap0 = self.hap
        hap1 = abs(self.hap - 1)
        if gt0 == 2 or gt1 == 2:
            hap0 += 1
            hap1 += 1
        
        gt_str = str(hap0) + '|' + str(hap1)

        if self.phased():
            if ('PS' in rec.FORMAT.split(':')):
                rec.samples[0].data = rec.samples[0].data._replace(GT=gt_str, PS=self.ps)
                
            else:
                rec.add_format('PS')
                samp_fmt = vcf.model.make_calldata_tuple(rec.FORMAT.split(':'))
                tmp = rec.samples[0].data._asdict()
                tmp['PS'] = self.ps
                tmp['GT'] = gt_str
                rec.samples[0].data = samp_fmt(**tmp)

            rec.samples[0].gt_nums = gt_str
        


class PhaseSet:
    def __init__(self, starting_pos):
        self.starting_pos = starting_pos
        self.records_idx = set()
        self.records = dict()
        self.origin = 2

    def contain_record(self, idx):
        return idx in self.records_idx
    
    def insert_record(self, record : Record):
        self.records_idx.add(record.pos)
        self.records[record.pos] = record

    def flip(self):
        for record in self.records.values():
            record.flip()

    def set_origin(self,o):
        self.origin = o

    def build_origin(self):
        f_support_count = 0
        m_support_count = 0
        m_supportintersect_phase_setcount = 0
        for rec in self.records.values():
            if rec.origin == 0:
                f_support_count += 1;
            if rec.origin ==1:
                m_support_count += 1
        if f_support_count >= m_support_count and f_support_count > 0:
            self.origin = 0
        if m_support_count > f_support_count and m_support_count > 0:
            self.origin = 1
        else:
            self.origin = 2
        
    def finalize_phaseset_label(self):
        self.starting_pos = min(self.records_idx) 
        for record in self.records.values():
            record.ps = self.starting_pos
    

    def intersect_phase_set(self, phase_set):
        intersection_length = 0
        hap0_hap0_support_length = 0
        hap0_hap1_support_length = 0
        hap1_hap0_support_length = 0
        hap1_hap1_support_length = 0
        # order: 0: 00 1:01, 2:10: 3:11
        # order = 0
        s_record: Record
        for s_record in phase_set.records.values():
            if s_record.pos in self.records_idx:
                intersection_length += 1
                record = self.records[s_record.pos]
                if record.hap == s_record.hap:
                    hap0_hap0_support_length += 1
                if record.hap == s_record.hap1:
                    hap0_hap1_support_length += 1
                if record.hap1 == s_record.hap:
                    hap1_hap0_support_length += 1
                if record.hap1 == s_record.hap1:
                    hap1_hap1_support_length += 1
        side = np.argmax([hap0_hap0_support_length, hap0_hap1_support_length, hap1_hap0_support_length, hap1_hap1_support_length])
        need_flip = False
        if side != 0 and side != 3:
            need_flip = True            
        # if intersection_length > 2:
        #     if hap0_support_length < hap1_support_length:
        #         need_flip = True
        return need_flip

class ChromosomoHaplotype:
    def __init__(self, in_vcf: vcf.Reader, chromo: str):
        self.chromo_record = dict()
        self.chromo_phase_set = dict()
        self.chromo_record2phaseset_map = dict()
        rec:vcf.model._Record
        ps_label_fix = dict()
        idx = 0
        
        for rec in in_vcf.fetch(chromo):
            het = rec.samples[0].gt_type
            if het != 1:        # not het loci
                continue
            PS_fix = 0
            if rec.samples[0].phased:
                fmt = rec.FORMAT.split(':')
                if 'PS' in fmt:
                    PS = rec.samples[0]['PS']
                    if PS in ps_label_fix.keys():
                        PS_fix = ps_label_fix[PS]
                    else:
                        ps_label_fix[PS] = rec.POS
                        PS_fix = rec.POS
                else:
                    PS_fix = 1
            record = Record()
            record.copy_from_rec(rec, PS_fix, idx)
            idx += 1
            self.chromo_record[record.pos] = record
            if record.ps != 0:
                PS = record.ps
                self.chromo_record2phaseset_map[record.pos] = PS
                phase_set : PhaseSet
                if PS in self.chromo_phase_set.keys():
                    phase_set = self.chromo_phase_set[PS]
                else:
                    phase_set = PhaseSet(record.ps)
                    self.chromo_phase_set[PS] = phase_set
                phase_set.insert_record(record)            

    def connect_phase_set(self, f_phase_set: PhaseSet, s_phase_set:PhaseSet):
        for record_pos, record in s_phase_set.records.items():
            f_phase_set.insert_record(record)
            self.chromo_record2phaseset_map[record_pos] = f_phase_set.starting_pos
        if s_phase_set.starting_pos in self.chromo_phase_set.keys():
            del self.chromo_phase_set[s_phase_set.starting_pos]
        s_phase_set.records.clear()
    
    def finalize_chromosome_haplotype(self):
        self.chromo_record2phaseset_map.clear()
        phase_set: PhaseSet
        temp = dict()
        for phase_set in self.chromo_phase_set.values():
            phase_set.finalize_phaseset_label()
            temp[phase_set.starting_pos] = phase_set
            for record_idx in phase_set.records_idx:
                self.chromo_record2phaseset_map[record_idx] = phase_set.starting_pos
        self.chromo_phase_set = temp

    