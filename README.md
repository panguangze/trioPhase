## require
python require, numpy,pysam,pyvcf


spechap https://github.com/deepomicslab/SpecHap


bcftools, tabix and bgzip


## Step by Step（Deprecate, please using One step cmd）
### Phase with bam 
#### child 
* step1 (If you have only father or mother data, cat is no need)
  
  extractHairs --vcf child.vcf --bam father.bam -o farhter_c.lst
  
  extractHairs --vcf child.vcf --bam mother.bam -o mother_c.lst
* step2

    cat farhter_c.lst mother_c.lst > parent.lst

    sort -n -k3 parent.lst > parent.sort.lst
* step3

    spechap -f parent.sort.lst -v child.vcf.gz -o child.out.vcf --keep_phasing_info

#### parent (mother or father)
* step1

  extractHairs --vcf father.vcf --bam child.bam -o child.lst
* step2

  sort -n -k3 child.lst > child.sort.lst
* step3

  spechap -f child.sort.lst -v father.vcf.gz -o father.out.vcf --keep_phasing_info

### Phase with only vcf
python merge_family.py -f xx.vcf -m xx.vcf -c xx.vcf -o out_dir
## One step
*\*_v means vcf file, \*_b means bam file*


python main.py --spechap xxx --extractHairs xxx --bcftools xxx --child_v xxx --child_b xxx --bgzip xxx --tabix xxx --out_dir xxx --ref xxx --script_root /home/caronkey/Documents/cityu/trio_phase --mother_v xx --mother_b xx --father_v xxx --father_b xxx


final output : *.phased.vcf.trio.vcf.final