#!/usr/bin/env nextflow

params.version=false
params.help=false
version='0.0.1'
timestamp='2019-06-04'

/**
	Prints version when asked for
*/
if (params.version) {
	System.out.println("")
	System.out.println("UK Biobank Genotype Quality Control - Version: $version ($timestamp)")
	exit 1
}


if (params.help) {
	System.out.println("")
	System.out.println("UK Biobank Genotype Quality Control -  Version: $version ($timestamp)")
	System.out.println("This pipeline is distributed in the hope that it will be useful")
	System.out.println("but WITHOUT ANY WARRANTY. See the MIT License for more details.")
	System.out.println("")
	System.out.println("Please report comments and bugs to https://github.com/choishingwan/uk-biobank-piplines")
	System.out.println("")
	System.out.println("Usage: ")
	System.out.println("   nextflow run perform_qc.nf [options]")
	System.out.println("")
	System.out.println("Mandatory arguments:")
	System.out.println("    --greed     Location of the GreedyRelated program")
	System.out.println("    --out       Prefix of the output files")
    System.out.println("    --sqc       Location of the UK biobank SQC file")
    System.out.println("    --rel       Location of the relatedness data file")
    System.out.println("    --bfile     Prefix of the genotype file. When provided, doesn't require")
    System.out.println("                the --bed --bim and --fam parameter")
    System.out.println("    --bed       Location of the bed file. Must use together with --bim and --fam")
    System.out.println("    --bim       Location of the bim file. Must use together with --bed and --fam")
    System.out.println("    --fam       Location of the fam file. Must use together with --bim and --bed")
	System.out.println("Options:")
	System.out.println("    --maf       Minor allele frequency filtering. Default=0.01")
	System.out.println("    --geno      Genotype missingness filtering. Default=0.02")
	System.out.println("    --mind      Individual missingness filtering. Default=0.02")
	System.out.println("    --hwe       HWE filtering. Default=1e-8")
	System.out.println("    --windSize  Window size for prunning. Default=200")
	System.out.println("    --windStep  Step size for prunning. Default=50")
    System.out.println("    --r2        Threshold for prunning. Default=0.2")
    System.out.println("    --maxSize   Maximum number of samples used for prunning. Default=10000")
    System.out.println("    --relThres  Threshold for removing related samples. Default=0.044")
    System.out.println("    --pheno     Phenotype file to guide greedy related")
    System.out.println("    --kmean     Number of cluster to form when defining ancestry. Default=4")
    System.out.println("    --sex       Can either be sd or fix. If sd, will exclude samples N sd away")
    System.out.println("                from mean, as defined by --sdm; whereas for fix, will exclude ")
    System.out.println("                samples based on --male and --female parameters. Default=sd")
    System.out.println("    --sdm       Define the number of SD away from the mean for sex filtering.")
    System.out.println("                Default=3")
    System.out.println("    --male      Define the F stat threshold for male. Male with F stat lower")
    System.out.println("                than this number will be removed. Default=0.8")
    System.out.println("    --male      Define the F stat threshold for female. Female with F stat higher")
    System.out.println("                than this number will be removed. Default=0.2")
    System.out.println("    --build     Genome build. Can either be grch37 or grch38. Use to define")
    System.out.println("                long LD regions. Default=grch37")
	System.out.println("")
    exit 1
}
params.maf = 0.01
params.geno = 0.02
params.mind = 0.02
params.hwe = 1e-8
params.windSize=200
params.windStep=50
params.r2=0.2
params.maxSize=10000
params.relThres=0.044
params.greed=false
params.bfile=false
params.bed=false
params.bim=false
params.fam=false
params.pheno=false
params.out=false
params.sqc=false
params.rel=false
params.kmean=4
params.dir="."
params.seed=false
params.male=0.8
params.female=0.2
params.build="grch37"
params.sex="sd"
params.sdm=3
// programs
greed=params.greed
// binary files
bfile=params.bfile
bed=params.bed
bim=params.bim
fam=params.fam
// Assign variables
build=params.build
maf=params.maf 
geno=params.geno
mind=params.mind
hwe=params.hwe
sdm=params.sdm
out=params.out
male=params.male
female=params.female
wind_size=params.windSize
wind_step=params.windStep
wind_r2=params.r2
max_size=params.maxSize
seed=params.seed
rel_thres=params.relThres
kmean=params.kmean
sexFilter=params.sex
if(!params.greed){
    error("GreedyRelated not provided!\n")
}
if(!params.bfile){
    if(!params.bed || !params.bim || !params.fam){
        error("Error: You must provide the bed bim and fam files!\n")
    }else if(params.bed && params.bim && params.fam){
    }else{
        error("Error: You must provide the binary file prefix!\n")
    }
}
if(!params.out){
    error("Error: You must provide an output prefix\n")
}
if(!params.sqc){
    error("Error: You must provide the SQC file\n")
}
if(!params.rel){
    error("Errorr: You must provide the relatedness file\n")
}
if(!params.seed){
    seed=1234
}
fileExists = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}
if(params.sex!="sd" && params.sex!="fix"){
    error("Error: Only support two type of sex based filtering, sd based or fixed\nsd based will filter samples 3 SD away from the mean, and fixed based filtering will filter male with F stat lower than 0.8 and Female with F stat higher than 0.2\n")
}
firstPass=Channel.create();
if(params.bfile){
    firstPass = Channel
        .fromFilePairs("${bfile}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  
        .ifEmpty { error "No matching plink files" }        
        // apply the fileExists function to the three input (bed bim fam)
        .map { a -> [fileExists(a[1]), fileExists(a[2]), fileExists(a[3])] } 
}else{
    firstPass = Channel
        .fromFilePairs("{${bed},${bim},${fam}}",size:3, flat: true)
        .ifEmpty { error "No matching plink files" }
        // apply the fileExists function to the three input (bed bim fam)
        .map { a -> [fileExists(a[1]), fileExists(a[2]), fileExists(a[3])] } 
}
      
sqcFile = Channel.fromPath(params.sqc)
relFile = Channel.fromPath(params.rel)
process firstPassGeno{
    cpus 12
    module 'bioinformatics/plink2/1.90b3.38'
    input:
    set file(bed), file(bim), file(fam) from firstPass

    output:
    file ("${out}-geno${geno}.snplist") into (fgeno)
    set file(bed), file(bim), file(fam) into basicQC

    script:

    """
    plink   --bed ${bed} \
            --bim ${bim} \
            --fam ${fam} \
            --geno ${geno} \
            --write-snplist \
            --out ${out}-geno${geno}
    """
}

process extractSQC{
    cpus 1
    module 'general/R/3.3.3'
    input:
    file(sqc) from sqcFile
    publishDir params.dir, mode:'copy', pattern:'*.covariates' 
    output:
    file("${out}.invalid") into hetProblem
    file("${out}.sex") into sexInfo
    file("${out}.covariates") into covariates  
    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    list.files()
    dat <- fread("${sqc}", data.table=F)
    cov <- dat[,c("FID", "IID", "Batch", paste("PC",1:40,sep=""))]
    write.table(cov, "${out}.covariates", quote=F, row.names=F)
    sex <- dat[,c("FID", "IID", "Submitted.Gender")]
    write.table(sex, "${out}.sex", quote=F, row.names=F)
    problem <- subset(dat,dat\$het.missing.outliers==1 | dat\$excess.relative==1 | FID<0)
    write.table(problem[,c("FID","IID")], "${out}.invalid", quote=F, row.names=F)
    """
}

process getEUR{
    cpus 1
    module 'general/R/3.3.3'
    input:
    file(qc) from covariates
    publishDir params.dir, mode:'copy', pattern:'*meanEUR'
    publishDir params.dir, mode:'copy', pattern:'*png'
    output:
    file("${out}.${kmean}meanEUR") into eur
    file("${out}.pca.png") into pcaPlot 
    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(ggplot2)
    cov <- fread("${qc}", data.table=F)
    # Set the seed for the kmean clustering
    if("${seed}"!="false"){
        set.seed("${seed}")
    }
    pc1k<-kmeans(cov\$PC1, ${kmean})
    pc2k<-kmeans(cov\$PC2, ${kmean})
    cov\$clusters<-as.factor(paste(pc1k\$cluster,pc2k\$cluster,sep="."))
    table(pc1k\$cluster, pc2k\$cluster)
    g <- ggplot(cov, aes(x=PC1,y=PC2, color=clusters))+\
            geom_point()+\
            theme_classic()
    ggsave("${out}.pca.png", plot=g, height=7,width=7)
    max.cluster <- names(which.max(table(cov\$clusters)))
    # We assume the largest cluster is the EUR cluster
    write.table(subset(cov, as.character(cov\$clusters)==max.cluster)[,c("FID", "IID")], "${out}.${kmean}meanEUR", quote=F, row.names=F)
    """
}



process basicQC{
    cpus 12
    module 'bioinformatics/plink2/1.90b3.38'
    input:
    set file(bed), file(bim), file(fam) from basicQC
    file(snplist) from (fgeno)
    file(eur) from eur
    file(invalid) from hetProblem 
    output:
    set file("${out}-qc.fam"), file("${out}-qc.snplist") into (qced)
    set file(bed), file(bim), file(fam) into (high_ld)
    script:
    base=bed.baseName
    """
    plink   --keep ${eur} \
            --bed ${bed} \
            --bim ${bim} \
            --fam ${fam} \
            --geno ${geno} \
            --maf ${maf} \
            --hwe ${hwe} \
            --write-snplist \
            --make-just-fam \
            --out ${out}-qc \
            --remove ${invalid} 
    """
}

process highLD{
    cpus 12
    module 'bioinformatics/plink2/1.90b3.38'
    input:
    set file(qc_fam), file(qc_snp) from qced
    set file(bed), file(bim), file(fam) from high_ld
    output:
    set file(qc_fam), file(qc_snp) into ldQC
    set file(bed), file(bim), file(fam) into prune
    file("${out}.set") into ld_set
    script:
    base=bed.baseName
    """
    echo "1     48000000     52000000   High_LD
2     86000000     100500000    High_LD
2     134500000     138000000   High_LD
2     183000000     190000000   High_LD
3     47500000     50000000 High_LD
3     83500000     87000000 High_LD
3     89000000     97500000 High_LD
5     44500000     50500000 High_LD
5     98000000     100500000    High_LD
5     129000000     132000000   High_LD
5     135500000     138500000   High_LD
6     25000000     35000000 High_LD
6     57000000     64000000 High_LD
6     140000000     142500000   High_LD
7     55000000     66000000 High_LD
8     7000000     13000000  High_LD
8     43000000     50000000 High_LD
8     112000000     115000000   High_LD
10     37000000     43000000    High_LD
11     46000000     57000000    High_LD
11     87500000     90500000    High_LD
12     33000000     40000000    High_LD
12     109500000     112000000  High_LD
20     32000000     34500000 High_LD" > high_ld_37
    echo "1     48060567     52060567     hild
2     85941853     100407914     hild
2     134382738     137882738     hild
2     182882739     189882739     hild
3     47500000     50000000     hild
3     83500000     87000000     hild
3     89000000     97500000     hild
5     44500000     50500000     hild
5     98000000     100500000     hild
5     129000000     132000000     hild
5     135500000     138500000     hild
6     25500000     33500000     hild
6     57000000     64000000     hild
6     140000000     142500000     hild
7     55193285     66193285     hild
8     8000000     12000000     hild
8     43000000     50000000     hild
8     112000000     115000000     hild
10     37000000     43000000     hild
11     46000000     57000000     hild
11     87500000     90500000     hild
12     33000000     40000000     hild
12     109521663     112021663     hild
20     32000000     34500000     hild
X     14150264     16650264     hild
X     25650264     28650264     hild
X     33150264     35650264     hild
X     55133704     60500000     hild
X     65133704     67633704     hild
X     71633704     77580511     hild
X     80080511     86080511     hild
X     100580511     103080511     hild
X     125602146     128102146     hild
X     129102146     131602146     hild" > high_ld_38
    ldFile=high_ld_37
    if [[ "${build}" != "grch37" ]];
    then
        ldFile=high_ld_38
    fi
    echo \${ldFile}
    plink \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --extract ${qc_snp} \
        --keep ${qc_fam} \
        --make-set \${ldFile} \
        --write-set \
        --out ${out}
    """
}
process prunning{
    cpus 12
    module 'bioinformatics/plink2/1.90b3.38'
    input: 
    set file(bed), file(bim), file(fam) from prune
    set file(qc_fam), file(qc_snp) from ldQC
    file(high_ld) from ld_set
    output:
    set file(bed), file(bim), file(fam) into sexCheck
    set file(qc_fam), file(qc_snp) into sexQC
    file("${out}-qc.prune.in") into pruned
    script:
    base=bed.baseName
    """
    if [[ \$(wc -l < ${qc_fam}) -ge ${max_size} ]];
    then
        plink \
            --bed ${bed} \
            --bim ${bim} \
            --fam ${fam} \
            --extract ${qc_snp} \
            --keep ${qc_fam} \
            --indep-pairwise ${wind_size} ${wind_step} ${wind_r2} \
            --out ${out}-qc \
            --thin-indiv-count ${max_size} \
            --seed ${seed} \
            --exclude ${high_ld}
    else
        plink \
            --bed ${bed} \
            --bim ${bim} \
            --fam ${fam} \
            --extract ${qc_snp} \
            --keep ${qc_fam} \
            --indep-pairwise ${wind_size} ${wind_step} ${wind_r2} \
            --out ${out}-qc \
            --exclude ${high_ld}
    fi 
    """
    }  

// currently, we assume the sex chromosome is always available
process calculateSexF{
    cpus 12
    module 'bioinformatics/plink2/1.90b3.38'
    input:
    set file(bed), file(bim), file(fam) from sexCheck
    file(prune) from pruned
    set file(qc_fam), file(qc_snp) from sexQC
    output:
    set file(bed), file(bim), file(fam) into endFilter
    set file(qc_fam), file(qc_snp) into sexFilterQC
    file("${out}.sexcheck") into sexF
    script:
    base=bed.baseName
    """
    plink \
        --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --extract ${prune} \
        --keep ${qc_fam} \
        --check-sex \
        --out ${out}
    """
}

process sexFiltering{
    cpus 1
    module 'general/R/3.3.3'
    input:
    file(sex) from sexInfo
    set file(qc_fam), file(qc_snp) from sexFilterQC
    file(sexF) from sexF
    output:
    file("${out}.invalid.sex") into sexFiltered
    file( "${out}.valid.rel") into relFilter
    set file(qc_fam), file(qc_snp) into postSexQC
    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    sex.info <- fread("${sex}")
    sex.fstat <- fread("${sexF}")
    fam <- fread("${qc_fam}")
    sex <- merge(sex.info, sex.fstat, by=c("FID", "IID")) 
    # we want to remove samples with negative FID (dropout) and samples 
    # with abnormal F statistics
    sex.summary <- data.table(Submitted.Gender=c("M","F"), lower=c(${male},-999),upper=c(999,${female}))
    if("${sexFilter}"=="sd"){
        sex.summary <- sex[FID>=0,list(avg_f = mean(F), sd_f = sd(F)), 'Submitted.Gender']
        sex.summary <- sex.summary[,.(lower=avg_f-${sdm}*sd_f, upper=avg_f+${sdm}*sd_f, Submitted.Gender)]    
    }
    sex.final <- merge(sex, sex.summary, by="Submitted.Gender")
    sex.invalid <- sex.final[(F<lower | F > upper)]
    fam.keep <- fam[!V1%in%sex.invalid\$FID & V1>=0]
    write.table(sex.invalid, "${out}.invalid.sex", quote=F, row.names=F)
    write.table(fam.keep, "${out}.valid.rel", quote=F, row.names=F)
    """
}

process relFiltering{
    cpus 1
    input: 
    file(rel) from relFile
    file(valid) from relFilter
    output:
    file("${out}.invalid.rel") into relFiltered

    script:
    """
    ${greed} \
        -r ${rel} \
        -i ID1 \
        -I ID2 \
        -f Kinship \
        -k ${valid} \
        -o ${out}.invalid.rel \
        -t ${rel_thres} \
        -s ${seed}
    """
}

process finalData{
    cpus 12
    module 'bioinformatics/plink2/1.90b3.38'
    input:
    file(sex) from sexFiltered
    set file(qc_fam), file(qc_snp) from postSexQC
    file(rel) from relFiltered
    set file(bed), file(bim), file(fam) from endFilter

    publishDir params.dir, mode:'copy', pattern:'*invalid*'
    publishDir params.dir, mode:'copy', pattern:'*-qc*'
    output:
    file("${out}-qc.snplist") into finalSNP
    file("${out}-qc.fam") into finalFAM

    script:
    base=bed.baseName
    """
    cat ${rel} ${sex} > ${out}.invalid.samples
    plink --bfile ${base} \
        --extract ${qc_snp} \
        --keep ${qc_fam} \
        --remove ${out}.invalid.samples \
        --make-just-fam \
        --write-snplist \
        --out ${out}-qc
    """


}