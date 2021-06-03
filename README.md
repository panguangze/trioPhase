## require
python require, numpy,pysam,pyvcf

spechap https://github.com/deepomicslab/SpecHap

extractHairs https://github.com/panguangze/extractHairs


tabix and bgzip
## Step by Step
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

python main.py --spechap /home/caronkey/Documents/cityu/triophase/SpecHap/cmake-build-debug/SpecHap --extractHairs /home/caronkey/Documents/cityu/hap/extracthairs/cmake-build-debug/ExtractHAIRs --child_v ~/remote/dong_hpc/data_X101SC19050094-Z01-B10-21/2.cleandata/s0114-1_FDHG190451812-1a/chr1.vcf.gz --child_b ~/remote/dong_hpc/data_X101SC19050094-Z01-B10-21/2.cleandata/s0114-1_FDHG190451812-1a/chr1.bam -o test2 --mother_b ~/remote/dong_hpc/data_X101SC19050094-Z01-B10-21/2.cleandata/s0114-1_FDHG190451812-1a/chr1.bam --mother_v ~/remote/dong_hpc/data_X101SC19050094-Z01-B10-21/2.cleandata/s0114-1_FDHG190451812-1a/chr1.vcf.gz --bgzip bgzip --tabix tabix --father_b ~/remote/dong_hpc/data_X101SC19050094-Z01-B10-21/2.cleandata/s0114-1_FDHG190451812-1a/chr1.bam --father_v ~/remote/dong_hpc/data_X101SC19050094-Z01-B10-21/2.cleandata/s0114-1_FDHG190451812-1a/chr1.vcf.gz


final output : *.phased.vcf.gz.final.vcf