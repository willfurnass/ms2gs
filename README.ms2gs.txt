# ms2gs
ms2gs reads genome data from ms file - format and perform forward simulation to evaluate genomic selection scenarios

# 2014-10-02 13:53:05 
! comments can be included as this
# or as this
# the command line is
#   ./ms2gs [-i parfile] [-h|-help] [-ms msfile] [-snp snpfile] [-iter niter] [-seed seed] [-quiet]
# results are printed to STDOUT
# sections in ms2gs.par file follow
NOBLUP    !--> does not perform blup
NOCAUSAL  !--> does not perform gblup with causal mutations
NOSEQ     !--> does not perform gblup with sequence data
NOCHIP    !--> does not perform gblup with snp data
NOFIT     !--> does not perform fit analysis
NITER    !--> no. of iterates, overriden by command line -iter
niter
MSFILE    !--> file with ms format data
msfile
MSCOMMAND !--> ms command that can be called from system, overrides MSFILE
mscommand
PLINKFILE  !--> if input is read from a pfile.ped and pfile.map file, either msfile or plinkfile
pfile
OUTFILE    !--> outfile prefix, overriden by command line -out
outfile
PEDFILE    !--> pedfile to carry out gene dropping w/o selection
pedfile
SNPFILE    !--> snpfile with snps to be used in chip, overriden by -snp in command line
snpfile    !    overrides ASCERTAINMENT section
SNPREGFILE !--> SNPREG FILE, if present, overrides ASCERTAINMENT section: chr, pos1, pos
snpregfile
QTLREGFILE !--> if present, overrides QTL_POS section: chr, pos1, pos
snpregfile
MAPFILE    !--> mapfile chr pos xrate : xrate from pos-1 to pos, default otherwise
mapfile
PRINT_SOL   !--> if present, prints snp solutions for every replicate and method
PRINT_PLINK !--> if present, prints plink file
PRINT_YG    !--> if present, prints phenotype, genotype and ebv file
NCHR        !--> no. of chrs specifies genome features, lengths in bp 
nchr 
CM2MB       !--> cm2Mb ratio (1 by default), can be changed with MAPFILE
cm2mb
CHRLENGHTS  !--> chr lengths in bp, all in one line
chr1_length chr2_length ...
H2          !--> heritability of the trait
h2
ADJUST_VE     !--> if present, var e is adjusted to match exactly expected h2, given observed haplotypes
BINARY_TRAIT  !--> if a threshold trait simulated, incidence; h2 provided in obs scale
incidence
KILLINDS      !--> delete these individuals for GBLUP
idk1 idk2
NQTL          !--> no. of qtls
nqtl
QTL_INDS      !--> QTL parameters computed using ind ind1 to ind2
qtl_ind1 qtl_ind2
QTL_SIGN      !--> probability of derived allele being deleterious [0.5]
pderived
SIMULATE_QTL_SNP_R2 !--> simulate qtl snps with r2 with nearest next left marker                  
r2 
SIMULATE_QTL_SNP_FNEUTRAL !--> simulate qtl snps with neutral frequency distribution
SIMULATE_QTL_SNP_FUNIFORM !--> simulate qtl snps with uniform[b1,b2]
b1 b2
QTL_POS !--> fixed qtl positions, ipos is snp order within chr
qtl1_chr qtl1_ipos
qtl2_chr qtl2_ipos
...
QTL_FREQ_RANGE !--> freq range on which potentially causal snps are chosen from
qfreq1, qfreq2
QTL_VA !--> fixed additive variance per qtl 
qtl1_va, qtl2_va, ...
QTL_VD  !--> fixed dominant variance per qtl 
qtl1_vd, qtl2_vd, ...
QTL_EFFECT_A !--> fixed additive effects per qtl 
qtl1_a, qtl2_a, ... qtl_nqtl_a
QTL_EFFECT_D  !--> fixed dominant effects per qtl 
qtl1_d, qtl2_d, ... 
QTL_DISTRIBUTION_A  !--> QTL add effects are sampled from a distribution: u(niform), g(amma), n(ormal)
[u, l_bound, u_bound], [n, mu, var], [g, s, b]
QTL_DISTRIBUTION_D !--> QTL dom effects are sampled from a distribution
[u, l_bound, u_bound], [n, mu, var], [g, s, b]
QTL_DISTRIBUTION_VA  !--> QTL add variances are sampled from a distribution
[u, l_bound, u_bound], [n, mu, var], [g, s, b]
QTL_DISTRIBUTION_VD  !--> QTL dom variances are sampled from a distribution
[u, l_bound, u_bound], [n, mu, var], [g, s, b]
ASCERTAINMENT_INDS  !--> snps are ascertained froms ind 1 to ind 2, a fraction p of snps with minimum maf are selected
asc_ind1, asc_ind2
ASCERTAINMENT_MAF  !--> snps are ascertained froms ind 1 to ind 2, a fraction p of snps with minimum maf are selected
maf
ASCERTAINMENT_P  !--> snps are ascertained froms ind 1 to ind 2, a fraction p of snps with minimum maf are selected
p
NXVAL     !-->no of xvalidation iterates
nxval
PXVAL     !--> % of samples idx1-idx2 deleted for xvalidation
p
IDXVAL    !--> sample ids range to perform xvalidation
idx1 idx2
SEED  !--> random seed can be overriden in command line by -seed (WARNING: does not affect ms program)
seed
KMIN !--> min allele count for a snp to be considered [1]
kmin
P_ERROR_CHIP !--> genotyping error [0]
perr_chip
P_ERROR_NGS  !--> base sequencing error [0]
perr_ngs
P_ERROR_IMP  !--> imputation error [0]
perr_imp
# the minimum command line is
#   ./ms2gs -i parfile

