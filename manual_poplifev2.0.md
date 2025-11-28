# **PopLife "v2.0" manual**

PopLife "v2.0" is a bioinformatic pipeline designed to rapidly estimate genetic diversity and selection statistics for any species or population using a genome-wide population VCF file. The pipeline also generates ready-to-visualize files in .bigwig format, which can be easily converted to HTML for interactive exploration.


## **1. Input: VCF files** 

To run PopLife "v2.0", you will need Variant Call Format (VCF) files for your population of interest. Detailed VCF specifications can be found at *https://samtools.github.io/hts-specs/VCFv4.2.pdf*

An example is located in /data/poplife/scratch/PopLifev2.0/1.input:

```bash
input=/data/poplife/scratch/PopLifev2.0/1.input

cd $input

zcat 24negw.31dog.2badog.5medog.32ibgw.14itgw.16bagw.chr38.vcf.gz | grep -v "^##" | cut -f1-10 | head -4
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	gr_wolf.russia.SAMEA7038713
chr38	16881	.	A	C	192.6	PASS	AC=0;AF=0.006289;AN=164;BaseQRankSum=0;ClippingRankSum=0;DP=1873;ExcessHet=0.0137;FS=6.264;InbreedingCoeff=0.1384;MLEAC=2;MLEAF=0.006289;MQ=44.31;MQRankSum=0;QD=17.51;RAW_MQ=21600;ReadPosRankSum=-0.431;SOR=0.061	GT:AD:DP:GQ:PL	./.:12,0:12:19:0,19,458
chr38	17451	.	G	A	7300.15	PASS	AC=28;AF=0.089;AN=214;BaseQRankSum=0.347;ClippingRankSum=0;DP=2258;ExcessHet=0.0011;FS=3.96;InbreedingCoeff=0.2926;MLEAC=29;MLEAF=0.092;MQ=60;MQRankSum=0;QD=18.81;RAW_MQ=54000;ReadPosRankSum=0.022;SOR=0.888	GT:AD:DP:GQ:PL	0/0:12,0:12:30:0,30,450
chr38	17532	.	T	A	238.27	PASS	AC=0;AF=0.006289;AN=212;BaseQRankSum=0.915;ClippingRankSum=0;DP=2462;ExcessHet=0.0138;FS=0;InbreedingCoeff=0.1915;MLEAC=2;MLEAF=0.006289;MQ=44.16;MQRankSum=0;QD=9.93;RAW_MQ=46800;ReadPosRankSum=0.94;SOR=0.755	GT:AD:DP:GQ:PL	0/0:12,0:12:36:0,36,449
...

```

The file *24negw.31dog.2badog.5medog.32ibgw.14itgw.16bagw.chr38.vcf.gz* contains 5 populations:
- **24 Northern European gray wolves (negw)** from Scandinavia and Russia (Kardos et al. 2018 Nat Ecol Evol, doi: 10.1038/s41559-017-0375-4)
- **38 dogs** from several European breeds (multiple sources, primarily from the Dog10K Consortium)
- **32 Iberian gray wolves** (Sarabia, Salado et al. 2025 Mol Ecol, doi: 10.1111/mec.17639)
- **14 Italian gray wolves** (Battilani et al. 2025 J Hered, doi: 10.1093/jhered/esae041)
- **16 Balkan gray wolves** (Ciucani et al. 2023 iScience, doi: 10.1016/j.isci.2023)

### **VCF requirements: PL field**

To proceed, your VCF file must include a **PL (Phred-scaled Genotype Likelihoods) field**. This field can be generated from a BAM file using bcftools mpileup or GATK HaplotypeCaller. Alternatively, if your VCF contains GL fields, you can convert them to PL using a Python script with the cyvcf2 library, or with bcftools +tag2tag:

```bash
Usage: bcftools +tag2tag [General Options] -- [Plugin Options]
Options:
   run "bcftools plugin" for a list of common options

Plugin options:
       --gp-to-gl           convert FORMAT/GP to FORMAT/GL
       --gp-to-gt           convert FORMAT/GP to FORMAT/GT by taking argmax of GP
       --gl-to-pl           convert FORMAT/GL to FORMAT/PL
       --pl-to-gl           convert FORMAT/PL to FORMAT/GL
   -r, --replace            drop the source tag
   -t, --threshold <float>  threshold for GP to GT hard-call [0.1]

Example:
   bcftools +tag2tag in.vcf -- -r --gp-to-gl
```

**Why is the PL field necessary?** ANGSD produces more accurate results when using the --vcf-pl option, as we will see in the next step.


### **Extracting populations**

You can easily extract specific populations using bcftools:

```bash
lists=/data/poplife/scratch/PopLifev2.0/0.lists
input=24negw.31dog.2badog.5medog.32ibgw.14itgw.16bagw.chr38.vcf.gz
path=/data/poplife/scratch/PopLifev2.0/1.input

bcftools view -S $lists/32ibgw.list $path/$input | bgzip > 32ibgw.chr38.vcf.gz
tabix 32ibgw.chr38.vcf.gz
```


## **2. Site Frequency Spectrum (SFS) with ANGSD** 

ANGSD (Korneliussen et al. 2014 BMC Bioinformatics, doi: 10.1186/s12859-014-0356-4) is a software package designed for analyzing low-coverage genomic data. It utilizes genotype likelihoods rather than called genotypes, offering greater robustness. A suite of complementary tools (ngsTools) is available (*https://github.com/mfumagalli/ngsTools*). 

This tutorial focuses on calculating genetic diversity statistics (e.g., Tajima's D, pairwise nucleotide diversity π, Watterson's θ) derived from the Site Frequency Spectrum (SFS) using ANGSD and realSFS. For advanced applications, refer to the ANGSD tutorial: *https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests*

Within the PopLifev2.0 environment, activate the ANGSD conda environment:

```bash
conda init && conda activate angsd
```

Alternatively, you can use the installation at /data/poplife/scratch/programs/angsd_def/angsd or download ANGSD from the official repository: *https://github.com/ANGSD/angsd*. 

### **2.1. Generate Sample Allele Frequency (SAF) Files**

You will need a reference genome for your species. Example usage:

```bash
angsd=/data/poplife/scratch/programs/angsd_def/angsd/angsd
path=/data/poplife/scratch/PopLifev2.0

# If we want to only calculate SFS and statistics for the 32ibgw population:
file=32ibgw.chr38

$angsd -P 8 -vcf-pl $path/1.input/$file.vcf.gz -anc $path/0.reference/canfam31.chr38.fa \
   -ref $path/0.reference/canfam31.chr38.fa \
	-minMapQ 20 -minQ 20 -minInd 14 \
	-doSaf 1 -doMajorMinor 1 -out $path/2.saf/$file
```

Parameters:

- *P 8 -minMapQ 20 -minQ 20*: Use 8 parallel threads, minimum mapping quality of 20 and minimum
base quality of 20 
- *doSaf 1*: Calculate site allele frequency likelihoods
- *doMajorMinor 1*: Infer major and minor alleles from read counts. 
- *minInd 14*: Include only sites with data from at least 14 individuals. Adjust this based on your population size and desired missing data threshold (e.g., -minInd 32 for 32 individuals would require no missing data).


### **2.2. Estimate the Site Frequency Spectrum (SFS)**

Use realSFS to compute the SFS from the SAF file:

```bash
realsfs=/data/poplife/scratch/programs/angsd_def/angsd/misc/realSFS
path=/data/poplife/scratch/PopLifev2.0

$realsfs $path/2.saf/32ibgw.chr38.saf.idx -P 10 -fold 1 > $path/2.sfs/32ibgw.chr38.sfs
```

This command computes the folded (-fold 1) site frequency spectrum using 10 threads (-P 10). 


### **2.3. Calculate Diversity Statistics (θ) based on SFS**

Convert the SAF file to the theta format and compute window-based statistics using the saf2theta function in realSFS:

```bash
# Convert SAF to Theta format using SFS
$realsfs saf2theta $path/2.saf/32ibgw.chr38.saf.idx -P 2 -fold 1 \
-sfs $path/2.sfs/32ibgw.chr38.sfs -outname $path/3.thetas/32ibgw.chr38

thstat=/data/poplife/scratch/programs/angsd_def/angsd/misc/thetaStat

# Calculate statistics per 50000bp window (using a step of 50000bp too):
$thstat do_stat $path/3.thetas/32ibgw.chr38.thetas.idx -win 50000 -step 50000 \
-outnames $path/3.thetas/32ibgw.chr38.w50k.st50k
```

You can specify sliding windows by using a -step value smaller than -win (e.g., -win 50000 -step 10000). 

Inspect the output:

```bash
# Check the values for different statistics:
head $path/3.thetas/32ibgw.chr38.w50k.st50k.pestPG
```

**Output Column Descriptions:**

| Column   | Description   |
| -------- | ------------  |
| (indexStart,indexStop)   | Variant indices included in the window   |
| (firstPos_withData,lastPos_withData)   | Genomic coordinates of the first and last site with data in the window   |
| (WinStart,WinStop)   | Start and end positions of the window   |
| Chr   | Chromosome   |
| WinCenter   | Genomic position of the window center   |
| tW   | Watterson's theta (θ) per window  |
| tP   | Pairwise theta (π) per window   |
| tF   | Fu and Li's tF estimator for theta per window   |
| tH   | Fay and Wu's H per window   |
| tL   | Tajima's thetaL per window   |
| Tajima   | Tajima's D   |
| fuf   | Fu and Li's F   |
| fud   | Fu and Li's D   |
| fayh   | Fay and Wu's H   |
| zeng   | Zeng's E   |
| nSites   | Number of sites in window   |

**Important Note**: The raw values for tW, tP, tF, tH, and tL are summed across the entire window and are not normalized per site. To obtain per-site estimates, divide these values by nSites.

**Best Practice**: Filter out windows with a very low number of sites (nSites), as they can produce extreme and unreliable diversity estimates. Excluding windows in the bottom 5% of nSites is recommended


## **3. F<sub>ST</sub> Calculation** 

Suppose you have a VCF file containing three populations:

```bash
pop1.pop2.pop3.vcf.gz
```

Create individual list files for each population:

```bash
head pop1.list
ind1pop1
ind2pop1
ind3pop1
...

head pop2.list
ind1pop2
ind2pop2
ind3pop2
...

head pop3.list
ind1pop3
ind2pop3
ind3pop3
...

```

### **3.1. Prepare Population Pairs and Phenotype Files**

For pairwise F<sub>ST</sub> analysis, create combined sample lists and PLINK phenotype files:

```bash
# Combine sample lists for each pair
cat pop1.list pop2.list > pop1.pop2.txt
cat pop1.list pop3.list > pop1.pop3.txt
cat pop2.list pop3.list > pop2.pop3.txt

# Create PLINK phenotype files (IMPORTANT: POP cannot start with a number)
echo -e "FID\tIID\tPOP" > head.list
cat pop1.list | awk -v OFS='\t' '{print "0", $1, "pop1"}' > pop1.pheno 
cat pop2.list | awk -v OFS='\t' '{print "0", $1, "pop2"}' > pop2.pheno
cat pop3.list | awk -v OFS='\t' '{print "0", $1, "pop3"}' > pop3.pheno

# Assemble pheno files
cat head.list pop1.pheno pop2.pheno > pop1.pop2.pheno
cat head.list pop1.pheno pop3.pheno > pop1.pop3.pheno
cat head.list pop2.pheno pop3.pheno > pop2.pop3.pheno
```

### **3.2. Extract Population Pairs and Split by Chromosome**

Extract VCFs for each population pair to ensure consistent sites, and split them by chromosome:

```bash
# Extract pairs (no allele count [AC] filtering to keep all sites consistent)
bcftools view -S pop1.pop2.list -Oz -o pop1.pop2.vcf.gz pop1.pop2.pop3.vcf.gz && tabix pop1.pop2.vcf.gz
bcftools view -S pop1.pop3.list -Oz -o pop1.pop3.vcf.gz pop1.pop2.pop3.vcf.gz && tabix pop1.pop3.vcf.gz
bcftools view -S pop2.pop3.list -Oz -o pop2.pop3.vcf.gz pop1.pop2.pop3.vcf.gz && tabix pop2.pop3.vcf.gz

# Split each pair VCF by chromosome (assuming n chromosomes):
for pop in pop1.pop2 pop1.pop3 pop2.pop3; do
   for j in {1..n}; do echo chr$j
      bcftools view -r chr$j $pop.vcf.gz | bgzip > $pop.chr$j.vcf.gz
      tabix $pop.chr$j.vcf.gz 
   done
done
```

### **3.3. Run F<sub>ST</sub> with PLINK2**

PLINK2 can calculate both Weir & Cockerham's (1984) and Hudson's (1992) F<sub>ST</sub>.

```bash
plink=/data/poplife/scratch/programs/plink2
vcf_dir=/path_to_your_VCF
out_dir=/path_to_fst_files
pheno_dir=/path_to_pheno_lists

for pop_pair in pop1.pop2 pop1.pop3 pop2.pop3; do
  for chr in {1..n}; do
    echo "Processing $pop_pair, chromosome $chr"

    # Create PLINK binary files
    $plink --vcf $vcf_dir/$pop_pair.chr$chr.vcf.gz --make-bed \
           --out $out_dir/$pop_pair.chr$chr.plink --chr-set 38 no-xy
    # Note: --chr-set 38 no-xy specifies 38 autosomes, excluding sex chromosomes.
    # Adjust this for your species (e.g., use 'no-x' for Drosophila, 'no-xy' for mammals, 'no-wz' for birds).

    # Calculate F_ST using both methods
    for method in hudson wc; do
      $plink --bfile $out_dir/$pop_pair.chr$chr.plink --chr-set 38 no-xy \
             --pheno $pheno_dir/$pop_pair.pheno \
             --fst POP method=$method report-variants \
             --out $out_dir/$pop_pair.chr$chr.$method
    done
  done
done

```

## **4. Population Branch Statistic (PBS)** 

PBS measures lineage-specific divergence based on F<sub>ST</sub> values among three populations. The divergence time (T<sub>div</sub>) between two populations is calculated as T<sub>pop1_pop2</sub> = -ln(1 - F<sub>ST_pop1_pop2</sub>).

A custom Python script is available to compute PBS at /data/poplife/scratch/programs/homescripts. 

You will probably have to modify the names of the output files from the previous PLINK2 step, since the output fst.var names will include the name of each population. I am using hudson as an example name. 

```bash
calc_script=/data/poplife/scratch/programs/homescripts/calculate_pbs.py
fst_dir=/path/to/fst/files

for chr in {1..n}; do
  python $calc_script \
    $fst_dir/pop1.pop2.chr$chr.hudson.fst.var \
    $fst_dir/pop1.pop3.chr$chr.hudson.fst.var \
    $fst_dir/pop2.pop3.chr$chr.hudson.fst.var \
    pop1 pop2 pop3 chr$chr
done

```

This script generates three files, each centered on one population:

- PBS.pop1.vs.pop2.vs.pop3.chr$chr.tsv (centered around pop1)
- PBS.pop2.vs.pop1.vs.pop3.chr$chr.tsv (centered around pop2)
- PBS.pop3.vs.pop1.vs.pop2.chr$chr.tsv (centered around pop3)


## **5. Phasing** 

Phasing infers the haplotype phase of genotypes. While not strictly required for selscan v2.0, it can improve accuracy. The following example uses Beagle to phase a two-population Bos taurus dataset across 29 chromosomes; you may need to adjust for your species.


```bash
chr="$1" # Chromosome passed as a script argument

input_vcf="/path/to/unphased/vcf/15bsbt.15htbt.chr${chr}.vcf.gz"
output_prefix="/path/to/phased/vcf/15bsbt.15htbt.chr${chr}.phased" # Beagle automatically adds ".vcf.gz"
beagle_jar="/data/poplife/scratch/programs/beagle5.5/beagle.27Feb25.75f.jar"

echo "Phasing chromosome ${chr}"

java -jar ${beagle_jar} gt=${input_vcf} out=${output_prefix} nthreads=20 iterations=12
```

Optional Performance Tuning:

```bash
temp_dir=/path/to/temporary_files
java -Xmx512m -Xss4m -XX:ParallelGCThreads=14 -Djava.io.tmpdir=$temp_dir \
     -jar ${beagle_jar} gt=${input_vcf} out=${output_prefix} nthreads=20 iterations=12 > beagle_chr${chr}.log
```

- -Xmx512m: Sets the maximum Java heap size to 512 MB. This controls how much memory the JVM can use for object allocation (i.e., main working memory for Java objects).
- -Xss4m: Sets the thread stack size to 4 MB. This affects how deep recursive method calls can go before causing a stack overflow.
- -XX:ParallelGCThreads=14: Tells the JVM to use 14 threads for garbage collection (GC). This is useful for multi-core machines to speed up GC pauses.
- -Djava.io.tmpdir=$temp: Sets the temporary directory. This is where Java will write temp files (used, for example, in file buffering or decompression).


## **6. Selection Scans with selscan**

selscan identifies genomic regions under selection using various haplotype-based statistics:

| Statistic   | Full Name   | Key Reference  |
| -------- | ------------  | ------------  | 
| **EHH**   | Extended Haplotype Homozygosity   | Sabeti et al. 2002 *Nature*   |
| **iHS**   | integrated Haplotype Score   | Voight et al. 2006 *PLoS Biol*   |
| **nSL**   | number of Segregating sites by Length   | Ferrer-Admetlla et al. 2014 *MBE*   |
| **XPEHH**   | Cross-Population EHH   | Sabeti et al. 2007 *Nature*   |
| **XPnSL**   | Cross-Population nSL   | Szpiech et al. 2021 *Evol Lett*   |

selscan incorporates pi too, but the window size cannot be determined (we can do it in ANGSD). 

Normalization is required after running selscan per chromosome. norm is an executable incorporated in the selscan package to normalize files. 

### **6.1. Create a Genetic Map File**

selscan requires a map file for each chromosome. First, extract physical positions from your VCF:

```bash
vcf_dir=/path/to/vcf/files
zcat $vcf_dir/pop1.chr1.vcf.gz | grep -v "^#" | awk '{print $1"\t"$3"\t"$2"\t"$2}' \
 > $vcf_dir/pop1.chr1.map
```

.map files are genetic map files that selscan need to use to calculate the statistics. If you want to have **physical map** files, you must find them in literature. There exists a number of maps for a diverse variety of species, especially for model species (mice, rats, C. elegans, Drosophila...) for humans (hg19, hg38...) and for domestic species (dogs, cows, chicken, goats...). An example for the physical map of human chromosome 22 in PLINK format is:

```bash
# PLINK format: Chr SNP_Name Genetic_Dist(cM) Physical_Pos(bp)
22	.	1.457757	15287922
22	.	1.585395	16370978
22	.	1.585931	16372046
22	.	1.586426	16373044
22	.	1.586781	16373740
```

Where the first column is the chromosome name, the second column is the name of the SNP, the third is the physical distance (in cM) and the fourth column is the genetic distance (in nucleotides from the beginning of the chromosome)

Interpolating a PLINK physical map is easy with the use of an in-home script that draws a regression line between genetic and physical distances and interpolates it to the SNP dataset that you have in your VCF, making 0 those distances drawn before the first position of the PLINK map and your dataset. There is a script at /data/poplife/scratch/programs/homescripts that can do it. Its usage is:

```bash
# Provide interpolate_genmap_ends.py with a 1) physical map (PLINK format), 2) genetic map (extracted from VCF in previous step), 3) a name for the output map.
python3.12 interpolate_genmap_ends.py PLINK.physical.map chr1.map interpolated_plink_genmap.chr1.map
```

If no genetic map is available, selscan can use the physical map extracted from the VCF, though results may be less precise.

### **6.2. Single-Population Statistics (iHS, nSL)**

```bash
selscan=/data/poplife/scratch/programs/selscan/bin/linux/selscan-2.0.3
vcf_dir=/path/to/vcf/files
map_dir=/path/to/map/files

# For nSL
$selscan --nsl --vcf $vcf_dir/pop1.chr1.vcf.gz --map $map_dir/pop1.chr1.map \
         --maf 0.001 --threads 20 --out $vcf_dir/pop1.chr1

# For iHS
$selscan --ihs --vcf $vcf_dir/pop1.chr1.vcf.gz --map $map_dir/pop1.chr1.map \
         --maf 0.001 --threads 20 --out $vcf_dir/pop1.chr1
```

Use flags like --trunc-ok, --max-extend, or --keep-low-freq for datasets with large populations but few SNPs. Output files will be *.nsl.out and *.ihs.out. 

### **6.3. Cross-Population Statistics (XP-EHH, XP-nSL)**

Populations must be differentiated in pop1, pop2, pop3. As with the single population statistics, a .map file is needed. The first population in the --vcf argument is treated as the "test" population, and the second (--vcf-ref) is the population against which we are running the statistic. This scripts detects what SNPs are putatively under selection in pop1 that are not found in pop2.

```bash
selscan=/data/poplife/scratch/programs/selscan/bin/linux/selscan-2.0.3
vcf_dir=/path/to/vcf/files
map_dir=/path/to/map/files

# For XP-EHH
$selscan --xpehh --vcf $vcf_dir/pop1.chr1.vcf.gz --vcf-ref $vcf_dir/pop2.chr1.vcf.gz \
         --map $map_dir/pop1.chr1.map --maf 0.001 --threads 20 --out $vcf_dir/pop1.vs.pop2.chr1

# For XP-nSL
$selscan --xpnsl --vcf $vcf_dir/pop1.chr1.vcf.gz --vcf-ref $vcf_dir/pop2.chr1.vcf.gz \
         --map $map_dir/pop1.chr1.map --maf 0.001 --threads 20 --out $vcf_dir/pop1.vs.pop2.chr1
```

### **6.4. Normalization**

Normalize scores across the genome in 100 Kb windows, binning by allele frequency (20 bins = 5% increments):

```bash
norm=/data/poplife/scratch/programs/selscan/bin/linux/norm

selscan_dir=/path/to/selscan/output
pop_prefix=pop1
stat_suffix=xpehh.out # Can be nsl.out, ihs.out, xpnsl.out, etc.

# Normalize all chromosomes simultaneously
$norm --xpehh --files $in/$pref.chr1.$suff $in/$pref.chr2.$suff $in/$pref.chr3.$suff $in/$pref.chr4.$suff $in/$pref.chr5.$suff $in/$pref.chr6.$suff $in/$pref.chr7.$suff $in/$pref.chr8.$suff $in/$pref.chr9.$suff $in/$pref.chr10.$suff $in/$pref.chr11.$suff $in/$pref.chr12.$suff $in/$pref.chr13.$suff $in/$pref.chr14.$suff $in/$pref.chr15.$suff $in/$pref.chr16.$suff $in/$pref.chr17.$suff $in/$pref.chr18.$suff $in/$pref.chr19.$suff $in/$pref.chr20.$suff $in/$pref.chr21.$suff $in/$pref.chr22.$suff  --winsize 100000 --bins 20
```

## **7. Convert Statistics to WIG Format**

For visualization, statistics need to be converted to WIG format, which will be converted to the compressed version, bigwig. Scripts are provided for different data types. All scripts are in folder /data/poplife/scratch/programs/homescripts

### **7.1. Genetic Diversity (from ANGSD)**

The script accepts: 1) A file name, 2) The output directory, 3) The output file name, 4) The chromosome number, 5) The name of the statistic. The names of the statistics are:

| tW   | Watterson's theta (θ) per window  |
| tP   | Pairwise theta (π) per window   |
| tF   | Fu and Li's tF estimator for theta per window   |
| tH   | Fay and Wu's H per window   |
| tL   | Tajima's thetaL per window   |
| Tajima   | Tajima's D   |
| fuf   | Fu and Li's F   |
| fud   | Fu and Li's D   |
| fayh   | Fay and Wu's H   |
| zeng   | Zeng's E   |

An usage example is:


```bash
# Arguments: FileName OutputDir OutputFileName Chromosome_number Statistic 
python $script $root/3.thetas/32ibgw.chr38.w50k.st50k.pestPG $root/4.wig/tajimaD tajimaD.32ibgw.chr38.wig 38 Tajima

```

### **7.2. PBS**

This script accepts: 1) Input file, 2) Directory where to write the output, and 3) Chromosome number. An usage example is:


```bash
script=/data/poplife/scratch/programs/homescripts/pbs2wig.py
root=/data/poplife/scratch/PopLifev2.0

# Arguments: InputFile OutputDir Chromosome_number
python pbs2wig.py $root/3.pbs/PBS.32ibgw.vs.14itgw.vs.16bagw.chr38.tsv $root/4.wig/pbs chr38
```

### **7.3. Single-Population Selection Statistics**

This script accepts: 1) The input file name, 2) The output directory, 3) The chromosome name. An usage example is:

```bash
script=/data/poplife/scratch/programs/homescripts/singlepopselection2wig.py
root=/data/poplife/scratch/PopLifev2.0

# Arguments: InputFile OutputDir Chromosome_number
python $script $root/2.singlescan/32ibgw.chr38.ihs.out.20bins.norm $root/4.wig/ihs 38
```

### **7.4. Cross-Population Selection Statistics**

This script accepts: 1) The input file name, 2) The output directory, 3) The chromosome name. An usage example is:

```bash
script=/data/poplife/scratch/programs/homescripts/crosspopselection2wig.py
root=/data/poplife/scratch/PopLifev2.0

# Arguments: InputFile OutputDir Chromosome_number
python $script $root/2.crosspopscan/32ibgw.vs.14itgw.chr38.xpehh.out.norm $root/4.wig/xpehh 38
```

## **8. Convert WIG to BigWIG Format**

Compress WIG files to BigWIG format for efficient storage and visualization using the UCSC tool wigToBigWig. This conversion requires a chromosome sizes file, which can be generated from .fa.fai index files. Below is an example for CanFam3.1:

```bash
# View .fai files containing chromosome lengths
head /data/poplife/scratch/csarabia/data/Canis_lupus/reference/canfam31.chr0*.fa.fai 
```

**Output:**

```bash
==> /data/poplife/scratch/csarabia/data/Canis_lupus/reference/canfam31.chr01.fa.fai <==
chr01	122678785	7	60	61

==> /data/poplife/scratch/csarabia/data/Canis_lupus/reference/canfam31.chr02.fa.fai <==
chr02	85426708	7	60	61

==> /data/poplife/scratch/csarabia/data/Canis_lupus/reference/canfam31.chr03.fa.fai <==
chr03	91889043	7	60	61
...
```

```bash
# Create a chromosome sizes file:
cat /path_to_reference/species.chr*.fa.fai | cut -f1,2 > /path_to_reference/chrom.sizes

# Example of CanFam3.1:
head /data/poplife/scratch/csarabia/data/Canis_lupus/reference/chrom.sizes
chr01	122678785
chr02	85426708
chr03	91889043
chr04	88276631
chr05	88915250
```

wigToBigWig requires three arguments: 1) An input WIG file, 2) A file with the names and sizes of all chromosomes, and 3) An output BigWIG file. An example usage:

```bash
ref=/data/poplife/scratch/csarabia/data/Canis_lupus/reference/chrom.sizes
inputfile=$root/4.wig/32ibgw.chr38.ihs.out.20bins.wig
outputfile=$root/5.bigwig/32ibgw.chr38.ihs.out.20bins.bigwigç
wig2bigwig=/data/poplife/scratch/programs/wigToBigWig

$wig2bigwig $inputfile $ref $outputfile
```

BigWIG files can be easily transferred from the server and visualized using PopLife.