### Raw phase with only vcf
python require, numpy,pysam,pyvcf
* step 1:

Install spechap https://github.com/deepomicslab/SpecHap
* step 2:

python main.py --father_v xx --mother_v xx --child_v xx --father_b xx --mother_b xx --child_b xx


### Further phase with bam (using vcf phased in raw phase step)
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

