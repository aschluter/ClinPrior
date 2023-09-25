#' Title
#'
#' @param variants object of class vcfR-class with the patient variants.
#' @param sampleName character name with the patient code as is written in the VCF file.
#' @param filter character name with the desired filter to select the variants. It should be as is written in the FILTER column in the vcf file. 
#' @param isCodingVar logical indicating whether to consider only variants in coding regions. Default TRUE.
#' @param frequenceAR numerical value indicating the MAF(minor allele frequence) to filter out the variants in autosomal recessive and X-linked inheritance. Default 0.01
#' @param frequenceAD numerical value indicating the MAF(minor allele frequence) to filter out the variants in autosomal dominant recessive inheritance. Default 0.00005
#' @param GlobalPhenotypicScore matrix with the phenotypic metrics obtained from the MatrixPropagation function.
#' @param assembly Genome assembly used. Default assembly human GRCh37.
#' @param processors Default=1.
#'
#' @return data frame with the patient's variants and their associated information classified from most to least likely to be the cause of the patient's phenotype.
#' @export
#'
#' @examples
#' library(vcfR)
#' patientHPOsFile <- paste(system.file("extdata/example", package = "ClinPrior"),"HPOpatient.txt",sep="/")
#' HPOpatient <- unique(read.csv(patientHPOsFile, header = FALSE, sep = "\t")[, 1])
#' Y<-proteinScore(HPOpatient)
#' ClinPriorGeneScore<-MatrixPropagation(Y,alpha=0.2)
#' 
#' vcfFile = paste(system.file("extdata/example", package = "ClinPrior"),"HG001_GRCh37_1_22_v4.2.1_benchmark.vep01.KCNQ2Met546Thr.vcf.gz",sep="/")
#' variants <- read.vcfR(vcfFile)
#' variantsFiltered <- readVCF(sampleName = "HG001",variants=variants)
#' result = priorBestVariantVcfR(variants = variantsFiltered, sampleName = "HG001",GlobalPhenotypicScore = ClinPriorGeneScore)
priorBestVariantVcfR<-function(variants, sampleName, filter = "",isCodingVar = TRUE, frequenceAR = 0.01,frequenceAD = 0.00005, GlobalPhenotypicScore, assembly="assembly37",processors=1)
{
  variants <- variants
  assembly <<- assembly
#load files
  dest = system.file(paste("extdata",assembly,sep = "/"), package = "ClinPrior")
  ClinPriorfiles<-list.files(path = dest, full.names = TRUE)
  CCRautosomesFile = ClinPriorfiles[grep("ccrs.autosomes.*gz$",ClinPriorfiles)]
  CCRxchromFile = ClinPriorfiles[grep("ccrs.xchrom.*gz$",ClinPriorfiles, perl = TRUE)]
  blacklistFile = ClinPriorfiles[grep("blacklist",ClinPriorfiles)]
  #CADDsnvFile = ClinPriorfiles[grep("CADD_SNV.*gz$",ClinPriorfiles, perl = TRUE)]
  #CADDindelFile = ClinPriorfiles[grep("CADD_Indels.*gz$",ClinPriorfiles, perl = TRUE)]
  #exomeFile = ClinPriorfiles[grep("exomes.*Below001.*gz$",ClinPriorfiles, perl = TRUE)]
  #genomeFile = ClinPriorfiles[grep("genomes.*Below001.*gz$",ClinPriorfiles, perl = TRUE)]
  #clinvarFile = ClinPriorfiles[grep("clinvar.*gz$",ClinPriorfiles, perl = TRUE)]

  loadVariants<-function(variants)
  {
    library(vcfR)
    ObjectVcfR<-read.vcfR(file=vcfFile)
      gt = variants@gt
      pos = match(sampleName, colnames(gt))
      gtPOS = grep("0/0",gt[,pos],invert=TRUE)
      variantsObject = variants[gtPOS]
   return(variantsObject)
  }


  #functions
  getAlt_forMutliAllelic<-function(variants)
  {

    gt = variants@gt
    pos = match(sampleName, colnames(gt))
    options(warn=-1)
    genotype = do.call(rbind, strsplit(as.character(gt[,pos]), ":"))[,1]
    options(warn=0)

    m<-regexpr("(\\d{1,2}\\/)|(\\d{1,2}\\|)",genotype,perl = TRUE)
    genotype1<-regmatches(genotype, m)
    genotype1<-gsub("\\||\\/", "",genotype1,perl=TRUE)

    m<-regexpr("(\\/\\d{1,2})|(\\|\\d{1,2})",genotype,perl = TRUE)
    genotype2<-regmatches(genotype, m)
    genotype2<-gsub("\\||\\/", "",genotype2,perl=TRUE)


    alternative = getFIX(variants)[,5]

    m<-regexpr("(\\w+\\,|\\w+)",alternative,perl = TRUE)
    allel1<-regmatches(alternative, m)
    allel1<-gsub("\\,", "",allel1,perl=TRUE)

    m<-regexpr("(\\,\\w+)",alternative,perl = TRUE)
    posAllel2<-grep("(\\,\\w+)",alternative,perl = TRUE)
    allel2<-regmatches(alternative, m)
    allel2<-gsub("\\,", "",allel2,perl=TRUE)

    #GENOTYPE2
    alt<-matrix(0,ncol = 1, nrow = dim(variants)[1])
    #0/1
    if(length(grep("1",genotype2))>0){alt[grep("1",genotype2)]<-allel1[grep("1",genotype2)]}
    #0/2
    if(length(grep("2",genotype2))>0){alt[posAllel2[grep("2",genotype2[posAllel2])]]<-allel2[grep("2",genotype2[posAllel2])]}
    #0/3
    if(length(grep("3",genotype2))>0){alt[posAllel2[grep("3",genotype2[posAllel2])]]<-allel2[grep("3",genotype2[posAllel2])]}
    if(length(grep("4",genotype2))>0){alt[posAllel2[grep("4",genotype2[posAllel2])]]<-allel2[grep("4",genotype2[posAllel2])]}
    if(length(grep("5",genotype2))>0){alt[posAllel2[grep("5",genotype2[posAllel2])]]<-allel2[grep("5",genotype2[posAllel2])]}
    if(length(grep("6",genotype2))>0){alt[posAllel2[grep("6",genotype2[posAllel2])]]<-allel2[grep("6",genotype2[posAllel2])]}
    if(length(grep("7",genotype2))>0){alt[posAllel2[grep("7",genotype2[posAllel2])]]<-allel2[grep("7",genotype2[posAllel2])]}
    if(length(grep("8",genotype2))>0){alt[posAllel2[grep("8",genotype2[posAllel2])]]<-allel2[grep("8",genotype2[posAllel2])]}
    if(length(grep("9",genotype2))>0){alt[posAllel2[grep("9",genotype2[posAllel2])]]<-allel2[grep("9",genotype2[posAllel2])]}
    if(length(grep("10",genotype2))>0){alt[posAllel2[grep("10",genotype2[posAllel2])]]<-allel2[grep("10",genotype2[posAllel2])]}


    #GENOTYPE1
    #alt<-matrix(0,ncol = 1, nrow = dim(variants)[1])
    #1/0
    if(length(grep("1",genotype1))>0){alt[grep("1",genotype1)]<-allel1[grep("1",genotype1)]}
    #2/0
    if(length(grep("2",genotype1))>0){alt[posAllel2[grep("2",genotype1[posAllel2])]]<-allel2[grep("2",genotype1[posAllel2])]}
    #3/0
    if(length(grep("3",genotype1))>0){alt[posAllel2[grep("3",genotype1[posAllel2])]]<-allel2[grep("3",genotype1[posAllel2])]}
    if(length(grep("4",genotype1))>0){alt[posAllel2[grep("4",genotype1[posAllel2])]]<-allel2[grep("4",genotype1[posAllel2])]}
    if(length(grep("5",genotype1))>0){alt[posAllel2[grep("5",genotype1[posAllel2])]]<-allel2[grep("5",genotype1[posAllel2])]}
    if(length(grep("6",genotype1))>0){alt[posAllel2[grep("6",genotype1[posAllel2])]]<-allel2[grep("6",genotype1[posAllel2])]}
    if(length(grep("7",genotype1))>0){alt[posAllel2[grep("7",genotype1[posAllel2])]]<-allel2[grep("7",genotype1[posAllel2])]}
    if(length(grep("8",genotype1))>0){alt[posAllel2[grep("8",genotype1[posAllel2])]]<-allel2[grep("8",genotype1[posAllel2])]}
    if(length(grep("9",genotype1))>0){alt[posAllel2[grep("9",genotype1[posAllel2])]]<-allel2[grep("9",genotype1[posAllel2])]}
    if(length(grep("10",genotype1))>0){alt[posAllel2[grep("10",genotype1[posAllel2])]]<-allel2[grep("10",genotype1[posAllel2])]}

    variants@fix[,5]<-alt

    return(variants)
  }

  deleteGnomadBlacklistVariants<-function(variants)
  {
    gt = variants@gt
    pos = match(sampleName, colnames(gt))
    options(warn=-1)
    genotype = do.call(rbind, strsplit(as.character(gt[,pos]), ":"))[,1]
    options(warn=0)

    m<-regexpr("(\\d{1,2}\\/)|(\\d{1,2}\\|)",genotype,perl = TRUE)
    genotype1<-regmatches(genotype, m)
    genotype1<-gsub("\\||\\/", "",genotype1,perl=TRUE)

    m<-regexpr("(\\/\\d{1,2})|(\\|\\d{1,2})",genotype,perl = TRUE)
    genotype2<-regmatches(genotype, m)
    genotype2<-gsub("\\||\\/", "",genotype2,perl=TRUE)


    alternative = getFIX(variants)[,5]

    m<-regexpr("(\\w+\\,|\\w+)",alternative,perl = TRUE)
    allel1<-regmatches(alternative, m)
    allel1<-gsub("\\,", "",allel1,perl=TRUE)

    m<-regexpr("(\\,\\w+)",alternative,perl = TRUE)
    posAllel2<-grep("(\\,\\w+)",alternative,perl = TRUE)
    allel2<-regmatches(alternative, m)
    allel2<-gsub("\\,", "",allel2,perl=TRUE)

    #shortVariant = do.call(paste, c(longVariant, sep="-"))
    #longVariant = gsub(" ", "", longVariant, fixed = TRUE)


    #GENOTYPE1
    alt<-matrix(0,ncol = 1, nrow = dim(variants)[1])
    #1/0
    if(length(grep("1",genotype1))>0){alt[grep("1",genotype1)]<-allel1[grep("1",genotype1)]}
    #2/0
    if(length(grep("2",genotype1))>0){alt[posAllel2[grep("2",genotype1[posAllel2])]]<-allel2[grep("2",genotype1[posAllel2])]}
    #3/0
    if(length(grep("3",genotype1))>0){alt[posAllel2[grep("3",genotype1[posAllel2])]]<-allel2[grep("3",genotype1[posAllel2])]}
    if(length(grep("4",genotype1))>0){alt[posAllel2[grep("4",genotype1[posAllel2])]]<-allel2[grep("4",genotype1[posAllel2])]}
    if(length(grep("5",genotype1))>0){alt[posAllel2[grep("5",genotype1[posAllel2])]]<-allel2[grep("5",genotype1[posAllel2])]}
    if(length(grep("6",genotype1))>0){alt[posAllel2[grep("6",genotype1[posAllel2])]]<-allel2[grep("6",genotype1[posAllel2])]}
    if(length(grep("7",genotype1))>0){alt[posAllel2[grep("7",genotype1[posAllel2])]]<-allel2[grep("7",genotype1[posAllel2])]}
    if(length(grep("8",genotype1))>0){alt[posAllel2[grep("8",genotype1[posAllel2])]]<-allel2[grep("8",genotype1[posAllel2])]}
    if(length(grep("9",genotype1))>0){alt[posAllel2[grep("9",genotype1[posAllel2])]]<-allel2[grep("9",genotype1[posAllel2])]}
    if(length(grep("10",genotype1))>0){alt[posAllel2[grep("10",genotype1[posAllel2])]]<-allel2[grep("10",genotype1[posAllel2])]}


    longVariant = cbind.data.frame(variants@fix[,c(1,2,4)],alt)
    longVariant = do.call(paste, c(longVariant, sep="-"))
    longVariant = gsub(" ", "", longVariant, fixed = TRUE)
    deleteVariants1<-match(ensemble,longVariant)[!is.na(match(ensemble,longVariant))]



    #GENOTYPE2
    alt<-matrix(0,ncol = 1, nrow = dim(variants)[1])
    #0/1
    if(length(grep("1",genotype2))>0){alt[grep("1",genotype2)]<-allel1[grep("1",genotype2)]}
    #0/2
    if(length(grep("2",genotype2))>0){alt[posAllel2[grep("2",genotype2[posAllel2])]]<-allel2[grep("2",genotype2[posAllel2])]}
    #0/3
    if(length(grep("3",genotype2))>0){alt[posAllel2[grep("3",genotype2[posAllel2])]]<-allel2[grep("3",genotype2[posAllel2])]}
    if(length(grep("4",genotype2))>0){alt[posAllel2[grep("4",genotype2[posAllel2])]]<-allel2[grep("4",genotype2[posAllel2])]}
    if(length(grep("5",genotype2))>0){alt[posAllel2[grep("5",genotype2[posAllel2])]]<-allel2[grep("5",genotype2[posAllel2])]}
    if(length(grep("6",genotype2))>0){alt[posAllel2[grep("6",genotype2[posAllel2])]]<-allel2[grep("6",genotype2[posAllel2])]}
    if(length(grep("7",genotype2))>0){alt[posAllel2[grep("7",genotype2[posAllel2])]]<-allel2[grep("7",genotype2[posAllel2])]}
    if(length(grep("8",genotype2))>0){alt[posAllel2[grep("8",genotype2[posAllel2])]]<-allel2[grep("8",genotype2[posAllel2])]}
    if(length(grep("9",genotype2))>0){alt[posAllel2[grep("9",genotype2[posAllel2])]]<-allel2[grep("9",genotype2[posAllel2])]}
    if(length(grep("10",genotype2))>0){alt[posAllel2[grep("10",genotype2[posAllel2])]]<-allel2[grep("10",genotype2[posAllel2])]}


    longVariant = cbind.data.frame(variants@fix[,c(1,2,4)],alt)
    longVariant = do.call(paste, c(longVariant, sep="-"))
    longVariant = gsub(" ", "", longVariant, fixed = TRUE)
    deleteVariants2<-match(ensemble,longVariant)[!is.na(match(ensemble,longVariant))]

    deleteVariants<-union(deleteVariants1,deleteVariants2)


    return(deleteVariants)
  }

  deleteIndelVariantsSamePos<-function(variants)
  {
    #delete indel variants in same position

    shortVariant = getFIX(variants)[,c(1,2)]
    shortVariant = apply(shortVariant[,1:2],1,function(x) paste(x,collapse="-"))

    m<-regexpr("(^\\S+\\-\\d+)",ensemble,perl = TRUE)
    shortEnsemble<-regmatches(ensemble, m)
    m<-regexpr("(\\w+\\-\\w+$)",ensemble,perl = TRUE)
    regexEnsembl<-regmatches(ensemble, m)
    refEnsembl = gsub("(\\w+)(-)(\\w+$)", "\\1", regexEnsembl, perl = TRUE)
    AltEnsembl = gsub("(\\w+)(-)(\\w+$)", "\\3", regexEnsembl, perl = TRUE)


    pos1<-match(shortEnsemble,shortVariant)[!is.na(match(shortEnsemble,shortVariant))]
    posInShortEnsembl<-c(1:length(shortEnsemble))[!is.na(match(shortEnsemble,shortVariant[pos1]))]

    # ref ==, Alt !=
    posRef = c(1:length(posInShortEnsembl))[nchar(variants@fix[pos1,4]) == nchar(refEnsembl[posInShortEnsembl])]
    posAlt = c(1:length(posInShortEnsembl))[nchar(variants@fix[pos1,5]) != nchar(AltEnsembl[posInShortEnsembl])]
    posIntersect1  = intersect(posRef,posAlt)
    # ref !=, Alt ==
    posRef = c(1:length(posInShortEnsembl))[nchar(variants@fix[pos1,4]) != nchar(refEnsembl[posInShortEnsembl])]
    posAlt = c(1:length(posInShortEnsembl))[nchar(variants@fix[pos1,5]) == nchar(AltEnsembl[posInShortEnsembl])]
    posIntersect2  = intersect(posRef,posAlt)

    # ref !=, Alt !=
    posRef = c(1:length(posInShortEnsembl))[nchar(variants@fix[pos1,4]) != nchar(refEnsembl[posInShortEnsembl])]
    posAlt = c(1:length(posInShortEnsembl))[nchar(variants@fix[pos1,5]) != nchar(AltEnsembl[posInShortEnsembl])]
    posIntersect3  = intersect(posRef,posAlt)


    deleteVariants<-pos1[unique(c(posIntersect1,posIntersect2,posIntersect3))]

    return(deleteVariants)
  }



  qualVariants<-function(qualityThreshold,variants)
  {
    posDiscard = 0
    posGQDiscard <- integer(0)
    posGQ<-match("GQ",strsplit(as.character(variants[1,9]),":")[[1]])
    if(!is.na(posGQ))
    {
      options(warn=-1)
      GeneQuality = do.call(rbind, strsplit(as.character(variants[,10]), ":"))[,posGQ]
      posGQDiscard<-c(1:length(GeneQuality))[as.numeric(GeneQuality)<=20]
      options(warn=0)
    }
    posAllelicDepth<-match("DP",strsplit(as.character(variants[1,9]),":")[[1]])
    options(warn=-1)
    deep = do.call(rbind, strsplit(as.character(variants[,10]), ":"))[,posAllelicDepth]
    options(warn=0)
    posDeepDiscard0<-grep("\\d+",deep,perl=TRUE,invert = TRUE)
    #deep<-sapply(variants[,10],function(x) strsplit(as.character(x),":")[[1]][posAllelicDepth])
    #deep<-sapply(deep,function(x) sum(as.numeric(strsplit(x,",")[[1]])))
    deep[posDeepDiscard0]<-0
    posDeepDiscard<-c(1:length(deep))[as.numeric(deep)<=10]

    #Discard bias <25% and >75% in HET variants
    posADvariants<-getADvariants()
    posAllelicDepth<-match("AD",strsplit(as.character(variants[1,9]),":")[[1]])
    options(warn=-1)
    deepAD = do.call(rbind, strsplit(as.character(variants[,10]), ":"))[,posAllelicDepth]
    options(warn=0)
    HETvariants<-deepAD[posADvariants]

    m<-regexpr("^(\\d{1,2}\\/\\d{1,2})|^(\\d{1,2}\\|\\d{1,2})",variants[,10],perl = TRUE)
    genotype<-regmatches(variants[,10], m)

    m<-regexpr("(\\d{1,2}\\/)|(\\d{1,2}\\|)",genotype,perl = TRUE)
    genotype1<-regmatches(genotype, m)
    genotype1<-gsub("\\||\\/", "",genotype1,perl=TRUE)
    genotype1<-genotype1[posADvariants]

    m<-regexpr("(\\/\\d)|(\\|\\d)",genotype,perl = TRUE)
    genotype2<-regmatches(genotype, m)
    genotype2<-gsub("\\||\\/", "",genotype2,perl=TRUE)
    genotype2<-genotype2[posADvariants]

    #GENOTYPE1
    freqBias<-matrix(0.5,ncol = 1, nrow = length(posADvariants))
    options(warn=-1)
    #1/0
    if(length(grep("1",genotype1))>0){freqBias[grep("1",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("1",genotype1)], ","))[,2])/as.numeric(deep)[posADvariants][grep("1",genotype1)]}
    #2/0
    if(length(grep("2",genotype1))>0){freqBias[grep("2",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("2",genotype1)], ","))[,3])/as.numeric(deep)[posADvariants][grep("2",genotype1)]}
    #3/0..
    if(length(grep("3",genotype1))>0){freqBias[grep("3",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("3",genotype1)], ","))[,4])/as.numeric(deep)[posADvariants][grep("3",genotype1)]}
    if(length(grep("4",genotype1))>0){freqBias[grep("4",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("4",genotype1)], ","))[,5])/as.numeric(deep)[posADvariants][grep("4",genotype1)]}
    if(length(grep("5",genotype1))>0){freqBias[grep("5",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("5",genotype1)], ","))[,6])/as.numeric(deep)[posADvariants][grep("5",genotype1)]}
    if(length(grep("6",genotype1))>0){freqBias[grep("6",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("6",genotype1)], ","))[,7])/as.numeric(deep)[posADvariants][grep("6",genotype1)]}
    if(length(grep("7",genotype1))>0){freqBias[grep("7",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("7",genotype1)], ","))[,8])/as.numeric(deep)[posADvariants][grep("7",genotype1)]}
    if(length(grep("8",genotype1))>0){freqBias[grep("8",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("8",genotype1)], ","))[,9])/as.numeric(deep)[posADvariants][grep("8",genotype1)]}
    if(length(grep("9",genotype1))>0){freqBias[grep("8",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("9",genotype1)], ","))[,10])/as.numeric(deep)[posADvariants][grep("9",genotype1)]}
    if(length(grep("10",genotype1))>0){freqBias[grep("8",genotype1)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("10",genotype1)], ","))[,11])/as.numeric(deep)[posADvariants][grep("10",genotype1)]}


    #GENOTYPE2
    freqBias2<-matrix(0.5,ncol = 1, nrow = length(posADvariants))
    #1/0
    if(length(grep("1",genotype2))>0){freqBias2[grep("1",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("1",genotype2)], ","))[,2])/as.numeric(deep)[posADvariants][grep("1",genotype2)]}
    #2/0
    if(length(grep("2",genotype2))>0){freqBias2[grep("2",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("2",genotype2)], ","))[,3])/as.numeric(deep)[posADvariants][grep("2",genotype2)]}
    #3/0..
    if(length(grep("3",genotype2))>0){freqBias2[grep("3",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("3",genotype2)], ","))[,4])/as.numeric(deep)[posADvariants][grep("3",genotype2)]}
    if(length(grep("4",genotype2))>0){freqBias2[grep("4",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("4",genotype2)], ","))[,5])/as.numeric(deep)[posADvariants][grep("4",genotype2)]}
    if(length(grep("5",genotype2))>0){freqBias2[grep("5",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("5",genotype2)], ","))[,6])/as.numeric(deep)[posADvariants][grep("5",genotype2)]}
    if(length(grep("6",genotype2))>0){freqBias2[grep("6",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("6",genotype2)], ","))[,7])/as.numeric(deep)[posADvariants][grep("6",genotype2)]}
    if(length(grep("7",genotype2))>0){freqBias2[grep("7",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("7",genotype2)], ","))[,8])/as.numeric(deep)[posADvariants][grep("7",genotype2)]}
    if(length(grep("8",genotype2))>0){freqBias2[grep("8",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("8",genotype2)], ","))[,9])/as.numeric(deep)[posADvariants][grep("8",genotype2)]}
    if(length(grep("9",genotype2))>0){freqBias2[grep("9",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("9",genotype2)], ","))[,10])/as.numeric(deep)[posADvariants][grep("9",genotype2)]}
    if(length(grep("10",genotype2))>0){freqBias2[grep("10",genotype2)]<-as.numeric(do.call(rbind, strsplit(as.character(HETvariants)[grep("10",genotype2)], ","))[,11])/as.numeric(deep)[posADvariants][grep("10",genotype2)]}

    #freqBias = as.numeric(do.call(rbind, strsplit(as.character(HETvariants), ","))[,1])/as.numeric(deep[posADvariants])
    #freqBias<-sapply(c(1:length(HETvariants)),function(x) as.numeric(strsplit(names(HETvariants[x]),",")[[1]][1])/HETvariants[x])

    options(warn=0)
    posnoNA<-c(1:length(freqBias))[!is.na(freqBias)]
    freqBias<-freqBias[posnoNA]
    posDiscard<-union(c(1:length(freqBias))[freqBias>0.75],c(1:length(freqBias))[freqBias<0.25])
    posBiasDiscard1<-c(1:length(deep))[posADvariants][posnoNA][posDiscard]

    posnoNA<-c(1:length(freqBias2))[!is.na(freqBias2)]
    freqBias2<-freqBias2[posnoNA]
    posDiscard<-union(c(1:length(freqBias2))[freqBias2>0.75],c(1:length(freqBias2))[freqBias2<0.25])
    posBiasDiscard2<-c(1:length(deep))[posADvariants][posnoNA][posDiscard]

    lengthMajor1 <- sapply(as.character(HETvariants), function(x) length(c(1:length(as.numeric(strsplit(x, ",")[[1]])))[as.numeric(strsplit(x, ",")[[1]])>1]))
    pos<-c(1:length(lengthMajor1))[lengthMajor1>2]
    posBiasDiscard3<-c(1:length(deep))[posADvariants][pos]


    posDiscard<-unique(c(posGQDiscard,posDeepDiscard,posBiasDiscard1,posBiasDiscard2,posBiasDiscard3))
    if(qualityThreshold>0)
    {
      pos100<-c(1:length(variants[,6]))[as.numeric(gsub(",",".",variants[,6]))<qualityThreshold]
      posDiscard<-union(posDiscard,pos100[!is.na(pos100)])
    }
    return(posDiscard)
  }

  getHOMvariants<-function(variants)
  {
    gt = variants@gt
    pos = match(sampleName, colnames(gt))
    posHOMvariants<-grep("1/1|2/2|3/3",gt[,pos])
    posHOMvariants<-posHOMvariants[grep("\\bX\\b",variants@fix[posHOMvariants,1],invert=TRUE)]
    return(posHOMvariants)
  }

  getPhysicalPhasing<-function(selectedVariants,variants)
  {
    gt = variants@gt
    pos = match(sampleName, colnames(gt))
    phaseVariants = matrix(".", ncol = 1, nrow = length(selectedVariants))
    posPhaseVariants<-grep("\\|",gt[selectedVariants,pos],perl =TRUE)

    if(length(posPhaseVariants)>0)
    {
      m<-regexpr("\\d+\\_\\S+\\_[A-Z]{1,300}",gt[selectedVariants,pos][posPhaseVariants],perl = TRUE)
      phaseVariants[posPhaseVariants] = regmatches(gt[selectedVariants,pos][posPhaseVariants], m)

    }
    return(phaseVariants)
  }

  adjustPhysicalPhasing<-function(posCmpHETvariants,variants)
  {
    gt = variants@gt
    pos = match(sampleName, colnames(gt))

    posPhaseVariants<-grep("\\|",gt[posCmpHETvariants,pos],perl =TRUE)
    posPhaseVariantsNOT<-grep("\\|",gt[posCmpHETvariants,pos],invert = TRUE,perl =TRUE)

    if(length(posPhaseVariants)>0)
    {
      m<-regexpr("\\d+\\_\\S+\\_[A-Z]{1,300}",gt[posCmpHETvariants,pos][posPhaseVariants],perl = TRUE)
      phaseVariants = regmatches(gt[posCmpHETvariants,pos][posPhaseVariants], m)

      posDups = match(unique(phaseVariants),phaseVariants)

      posCmpHETvariants<-sort(union(posCmpHETvariants[posPhaseVariantsNOT],posCmpHETvariants[posPhaseVariants][posDups]))
    }
    return(posCmpHETvariants)
  }

  getGTsamples<-function(selectedVariants,variants)
  {
    gt = variants@gt
    GTsamples=0
    for(x in 2:length(colnames(gt)))
    {
      options(warn=-1)
      GTsamples = cbind(GTsamples,do.call(rbind, strsplit(as.character(gt[selectedVariants,x]), ":"))[,1])
      options(warn=0)
    }
    #GTsamples = GTsamples[,-1]
    return(GTsamples)
  }

  getCmpHETvariants<-function(variants)
  {
    posADvariants<-getADvariants(variants)
    VariantAnnotation<-getVariant(posADvariants,variants)
    geneCmpHET<-getGeneSymbol(VariantAnnotation,variants)
    posUniqueRemove = match(names(table(geneCmpHET))[table(geneCmpHET)==1],geneCmpHET)
    posCmpHETvariants<-posADvariants[-posUniqueRemove]
    posADvariants<-adjustPhysicalPhasing(posCmpHETvariants,variants)

    VariantAnnotation<-getVariant(posADvariants,variants)
    geneCmpHET<-getGeneSymbol(VariantAnnotation,variants)

    posCmpHETvariants = posADvariants

    posNULL = c(1:length(geneCmpHET))[grep(".",geneCmpHET,invert = TRUE)]
    if(length(posNULL)>0){geneCmpHET = geneCmpHET[-posNULL]
    posCmpHETvariants = posADvariants[-posNULL]}

    posUniqueRemove = match(names(table(geneCmpHET))[table(geneCmpHET)==1],geneCmpHET)

    if(length(posUniqueRemove)>0){posCmpHETvariants<-posCmpHETvariants[-posUniqueRemove]}
    posCmpHETvariants<-posCmpHETvariants[grep("\\bX\\b",variants@fix[posCmpHETvariants,1],invert=TRUE)]


    #geneCmpHET<-getGeneName(VariantAnnotation)
    #dups<-duplicated(geneCmpHET)
    #geneDuplicated<-unique(geneCmpHET[dups])
    #options(warn=-1)
    #matchItem<-do.call(paste, c(cbind.data.frame("\\b",geneDuplicated,"\\b"), sep=""))
    #options(warn=0)
    #cl <- parallel::makeCluster(processors)
    #parallel::clusterExport(cl,c("geneCmpHET"))
    #parallel::clusterEvalQ(cl, library(igraph))

    #posGenes<-parSapply(cl,matchItem, function(x) grep(x,geneCmpHET))
    #stopCluster(cl)
    #print(Sys.time() - strt)
    #quant<-sapply(c(1:length(posGenes)), function(x) length(posGenes[[x]]))
    #posCmpHETvariants<-posADvariants[unlist(posGenes[quant<4])]

    return(posCmpHETvariants)
  }

  getADvariants<-function(variants)
  {
    gt = variants@gt
    pos = match(sampleName, colnames(gt))
    posADvariants<-grep("0/1|0/2|0/3|1/2|2/3|1/3",gt[,pos])
    posADvariants<-posADvariants[grep("\\bX\\b",variants@fix[posADvariants,1],invert=TRUE)]
    return(posADvariants)
  }

  getXLvariants<-function(variants)
  {
    gt = variants@gt
    pos = match(sampleName, colnames(gt))
    posXLvariants<-grep("0/0",gt[,pos],invert = TRUE)
    posXLvariants<-posXLvariants[grep("\\bX\\b",variants@fix[posXLvariants,1])]
    return(posXLvariants)
  }

  getVariant<-function(selectedVariants,variants)
  {
    VariantAnnotation = variants@fix[selectedVariants,8]
    return(VariantAnnotation)
  }

  getProtein<-function(VariantAnnotation,variants)
  {

    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    ProteinPOS = match("HGVSp",meta)

    options(warn=-1)
    ProtAnnotation = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,ProteinPOS]
    options(warn=0)
    #ProtAnnotation<-sapply(VariantAnnotation,function(x) strsplit(as.character(x),"[|]")[[1]][11])
    return(ProtAnnotation)
  }

  getcDNA<-function(VariantAnnotation,variants)
  {

    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    cDNAPOS = match("HGVSc",meta)

    options(warn=-1)
    cDNAannotation = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,cDNAPOS]
    options(warn=0)

    return(cDNAannotation)
  }

  getEnsemblGene<-function(VariantAnnotation,variants)
  {

    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    EnsemblPOS = match("Gene",meta)

    options(warn=-1)
    Ensemblannotation = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,EnsemblPOS]
    options(warn=0)

    return(Ensemblannotation)
  }

  getConsequence<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    ConsequencePOS = match("Consequence",meta)

    options(warn=-1)
    consequence = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,ConsequencePOS]
    options(warn=0)

    return(consequence)
  }

  getGeneName<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    symbolPOS = match("Gene",meta)

    options(warn=-1)
    geneName = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,symbolPOS]
    options(warn=0)


    geneName<-totalEnsembl[match(geneName,totalEnsembl[,3]),2]

    return(geneName)
  }

  getGeneSymbol<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    symbolPOS = match("Gene",meta)

    options(warn=-1)
    geneName = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,symbolPOS]
    options(warn=0)


    geneName<-totalEnsembl[match(geneName,totalEnsembl[,3]),1]

    return(geneName)
  }


  isCoding<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]##!!!!!!!
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    posCoding = match("BIOTYPE",meta)

    options(warn=-1)
    coding = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,posCoding]
    options(warn=0)

    posCoding = c(1:length(VariantAnnotation))[grep("\\bprotein_coding\\b",coding)]
    return(posCoding)
  }

  getisCoding<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    posCoding = match("BIOTYPE",meta)

    options(warn=-1)
    coding = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,posCoding]
    options(warn=0)

    return(coding)
  }

  isModerate<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    impactPOS = match("IMPACT",meta)

    options(warn=-1)
    impact = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,impactPOS]
    options(warn=0)

    posModerate<-grep("MODERATE",impact)
    ModerateMatrix<-matrix(".",ncol = 1,nrow = length(VariantAnnotation))
    if(length(posModerate)>0)
    {
      ModeratecDNA<-getcDNA(VariantAnnotation[posModerate],variants)
      ModerateProtein<-getProtein(VariantAnnotation[posModerate],variants)
      options(warn=-1)
      ModerateAnnotation<-do.call(paste, c(cbind.data.frame(ModeratecDNA,ModerateProtein), sep="|"))
      options(warn=0)
      #ModerateAnnotation<-sapply(c(1:length(posModerate)),function(x) paste(ModeratecDNA[x],ModerateProtein[x],sep = "|"))
      ModerateMatrix[posModerate]<-ModerateAnnotation
    }
    return(ModerateMatrix)
  }

  isHigh<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    impactPOS = match("IMPACT",meta)

    options(warn=-1)
    impact = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,impactPOS]
    options(warn=0)

    posHigh<-grep("HIGH",impact)

    #pos5UTR<-c(1:length(VariantAnnotation))[grep("5_prime_UTR",VariantAnnotation)]
    #posHigh<-setdiff(posHigh,pos5UTR)
    HighMatrix<-matrix(".",ncol = 1,nrow = length(VariantAnnotation))
    if(length(posHigh)>0)
    {
      HIGHcDNA<-getcDNA(VariantAnnotation[posHigh],variants)
      HighProtein<-getProtein(VariantAnnotation[posHigh],variants)
      options(warn=-1)
      HighAnnotation<-do.call(paste, c(cbind.data.frame(HIGHcDNA,HighProtein), sep="|"))
      options(warn=0)
      #HighAnnotation<-sapply(c(1:length(posHigh)),function(x) paste(HIGHcDNA[x],HighProtein[x],sep = "|"))
      HighMatrix[posHigh]<-HighAnnotation
    }
    return(HighMatrix)
  }

  isSplicingOld<-function(VariantAnnotation,shortVariant,longVariant, Threshold,variants)
  {
    SplicingMatrix<-matrix(".",ncol = 1,nrow = length(VariantAnnotation))
    consequence<-getConsequence(VariantAnnotation,variants)
    posSplicing1<-grep("splice",consequence,perl=TRUE)
    posSplicing2<-grep("intron_variant",consequence,perl=TRUE)
    posSplicing3<-grep("synonymous_variant",consequence,perl=TRUE)

    if(length(posSplicing2)>0)
    {
      branchpoint<-ClinPrior::branchpointer(longVariant[posSplicing2],assembly,0.2)
      posBranchPoint<-grep("branchpoint",branchpoint)
      if(length(posBranchPoint)>0){SplicingMatrix[posSplicing2][posBranchPoint]<-branchpoint[posBranchPoint]}
      cDNA<-getcDNA(VariantAnnotation[posSplicing2],variants)
      names(cDNA)<-c()
      m <- regexpr("\\+[\\d+]{1,8}|\\-[\\d+]{1,8}", cDNA,perl=TRUE)
      posMatch <- grep("\\+[\\d+]{1,8}|\\-[\\d+]{1,8}", cDNA,perl=TRUE)
      match<-regmatches(cDNA, m)
      distSplic<-as.numeric(gsub("\\+|\\-", "",match,perl=TRUE))
      posSplicThreshold<-c(1:length(distSplic))[distSplic<Threshold]
      posSplicing2<-posSplicing2[posMatch][posSplicThreshold]

      posSplicing<<-unique(c(posSplicing1,posSplicing2,posSplicing3))

      if(length(posSplicing)>0)
      {
        #SplicingcDNA<-getcDNA(VariantAnnotation[posSplicing])
        #strt <- Sys.time()
        posMT<<-grep("MT",shortVariant[posSplicing], invert = TRUE)
        if(length(posMT)>0)
        {
          VariantAnnotation<<-VariantAnnotation
          shortVariant<<-shortVariant
          longVariant<<-longVariant
          assembly<<-assembly
          if(length(posSplicing[posMT])>20)
          {
            cl <- parallel::makeCluster(processors)
            parallel::clusterExport(cl,c("VariantAnnotation","posSplicing","shortVariant","posMT","longVariant","assembly"))
            parallel::clusterEvalQ(cl, library(ClinPrior))

            SplicingMaxentScan<-parSapply(cl,c(1:length(posSplicing[posMT])), function(x){
              getGeneName<-function(VariantAnnotation,variants)
              {
                options(warn=-1)
                geneName = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,4]
                options(warn=0)
                geneName<-total_unique[match(geneName,total_unique[,1]),2]
                #geneName<-unname(sapply(VariantAnnotation,function(x) strsplit(as.character(x),"[|]")[[1]][4]))
                return(geneName)
              }

              #SplicingMaxentScan<-sapply(c(1:length(posSplicing[posMT])),function(x){
              geneName<-getGeneName(VariantAnnotation[posSplicing][posMT][x],variants)
              chr<-paste("chr",strsplit(shortVariant[posSplicing][posMT][x],":")[[1]][1],sep="")
              start<-as.numeric(strsplit(shortVariant[posSplicing][posMT][x],"-")[[1]][2])
              Ref<-strsplit(longVariant[posSplicing][posMT][x],"-")[[1]][2]
              Alt<-strsplit(longVariant[posSplicing][posMT][x],"-")[[1]][3]
              ClinPrior::MaxentScanClinPrior(geneName,chr,start,Ref,Alt,assembly,0.2)
            })
            stopCluster(cl)
            gc()
          }else{
            SplicingMaxentScan<-sapply(c(1:length(posSplicing[posMT])),function(x){
              geneName<-getGeneName(VariantAnnotation[posSplicing][posMT][x],variants)
              chr<-paste("chr",strsplit(shortVariant[posSplicing][posMT][x],":")[[1]][1],sep="")
              start<-as.numeric(strsplit(shortVariant[posSplicing][posMT][x],"-")[[1]][2])
              Ref<-strsplit(longVariant[posSplicing][posMT][x],"-")[[1]][2]
              Alt<-strsplit(longVariant[posSplicing][posMT][x],"-")[[1]][3]
              ClinPrior::MaxentScanClinPrior(geneName,chr,start,Ref,Alt,assembly,0.2)
            })}
          #print(Sys.time() - strt)
          discard<-c(1:length(posSplicing[posMT]))[colSums(SplicingMaxentScan)==0]
          SplicingMaxentScan<-SplicingMaxentScan[,-discard]
          options(warn=-1)
          SplicingMatrix[posSplicing][posMT][-discard]<-do.call(paste, c(as.data.frame(t(SplicingMaxentScan))))
          options(warn=0)

        }
      }
    }
    return(SplicingMatrix)
  }

  isSplicing<-function(VariantAnnotation,variants)
   {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    MaxEntAlt = match("MaxEntScan_alt",meta)
    MaxEntDiff = match("MaxEntScan_diff",meta)
    MaxEntRef = match("MaxEntScan_ref",meta)


    VariantAnnotation = paste(VariantAnnotation,"|",sep="")

      options(warn=-1)
      MaxEnt = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,c(MaxEntDiff,MaxEntRef,MaxEntAlt)]
      options(warn=0)

      MaxEnt = rbind(c("","",""),MaxEnt)
      percentImp = as.numeric(MaxEnt[,1])/as.numeric(MaxEnt[,2])
      pos = c(1:length(percentImp))[abs(percentImp)>=0.2][!is.na(c(1:length(percentImp))[abs(percentImp)>=0.2])]
      MaxEnt2 = matrix("",ncol = 3,nrow = length(percentImp))
      MaxEnt2[pos,] = MaxEnt[pos,]

      SplicingMatrix = apply(MaxEnt2,1,function(x) paste(x,collapse="|"))
      SplicingMatrix = SplicingMatrix[-1]

    return(SplicingMatrix)
  }

  clinvarPathoScore<-function(VariantScore,clinvar,phenoScore)
  {

    pos<-intersect(grep("Pathogenic|Likely_pathogenic",clinvar), c(1:length(VariantScore))[phenoScore[,1] <= 2000][!is.na(c(1:length(VariantScore))[phenoScore[,1] <= 2000])])
    #    if(length(pos)>0){VariantScore[pos] = VariantScore[pos] + 1}
    if(length(pos)>0){VariantScore[pos] = VariantScore[pos] + 0} #treure clinvar PATH score
    return(VariantScore)
  }

  #Calculte Variant score
  calculateVariantScore<-function(HIGH,Splicing,Moderate,CADD,CCR,clinvar,metrics,inheritance,gtPOS)
  {
    VariantScore = 0

    if(HIGH !=".")
    {
      VariantScore = VariantScore + 0.4
      if(!is.na(metrics$pLI)){if(metrics$pLI>0.9){VariantScore = VariantScore + 0.2}}
      if(!is.na(metrics$lof_z)){if(metrics$lof_z>3){VariantScore = VariantScore + 0.1}}
      if(CADD >=20){VariantScore = VariantScore + 0.05}
      if(CADD >=25){VariantScore = VariantScore + 0.05}
      if(CADD >=30){VariantScore = VariantScore + 0.1}
    }

    if((Moderate !=".") && (Splicing =="||") && (HIGH =="."))
    {
      VariantScore = VariantScore + 0.2
      if(!is.na(metrics$pLI)){if(metrics$pLI>0.9){VariantScore = VariantScore + 0.1}}
      if(!is.na(metrics$mis_z)){if(metrics$mis_z>3){VariantScore = VariantScore + 0.05}}
      if(CADD >=20){VariantScore = VariantScore + 0.05}
      if(CADD >=25){VariantScore = VariantScore + 0.05}
      if(CADD >=30){VariantScore = VariantScore + 0.1}
      if(CCR >=85){VariantScore = VariantScore + 0.1}
      if(CCR >=95){VariantScore = VariantScore + 0.2}
    }

    if((Splicing !="||") &&  (HIGH =="."))
    {
      VariantScore = VariantScore + 0.2
      if(!is.na(metrics$pLI)){if(metrics$pLI>0.9){VariantScore = VariantScore + 0.1}}
      if(!is.na(metrics$lof_z)){if(metrics$lof_z>3){VariantScore = VariantScore + 0.1}}
      if(CADD >=20){VariantScore = VariantScore + 0.1}
      if(CADD >=25){VariantScore = VariantScore + 0.1}
      if(CADD >=30){VariantScore = VariantScore + 0.1}

    }

    if(inheritance == "HOM")
    {
      VariantScore = VariantScore + 0.5
      if(!is.na(metrics$pRec)){if(metrics$pRec>0.7){VariantScore = VariantScore + 0.1}}
    }

    if(inheritance == "CmpHET")
    {
      VariantScore = VariantScore + 0.2
      if(!is.na(metrics$pRec)){if(metrics$pRec>0.7){VariantScore = VariantScore + 0.05}}
    }

    if((inheritance == "HET") && (Splicing =="||") && (HIGH ==".") && (CCR ==0)) {VariantScore = VariantScore/10}

    if((Splicing =="||") && (HIGH ==".") &&  (Moderate ==".")){VariantScore = 0}

    if((inheritance =="XLink") && (gtPOS =="0/1")){VariantScore = 0}

    #if(length(grep("benign",clinvar,ignore.case = TRUE))>0){VariantScore = 0}

    return(round(VariantScore, digits = 4))
  }

  createVCFshortVariantCCR<-function(selectedVariants,variants)
  {
    shortVariant <- getFIX(variants[selectedVariants])[,1:2]
    position <- as.numeric(getFIX(variants[selectedVariants])[,2])
    alt <- getFIX(variants[selectedVariants])[,5]

    #shortVariant = variants[selectedVariants,c(1,2)]
    shortVariant = apply(shortVariant,1,function(x) paste(x,collapse=":"))
    lengthALT<-as.numeric(sapply(alt,function(x) nchar(as.character(x))))
    shortVariant = apply(cbind(shortVariant,(position+lengthALT-1)),1,function(x) paste(x,collapse="-"))

    #process <- function(x) {
    #  shortVariant <- getFIX(variants[selectedVariants[x]])[1:2]
    #  position <- as.numeric(getFIX(variants[selectedVariants[x]])[2])
    #  alt <- getFIX(variants[selectedVariants[x]])[5]
    #  shortVariant <- paste(shortVariant,collapse=":")
    #  lengthALT <- as.numeric(nchar(as.character(alt)))
    #  shortVariant <- paste(cbind(shortVariant,(position+lengthALT-1)),collapse="-")
    #}

    #shortVariant <-do.call(rbind, lapply(c(1:length(selectedVariants)), process))
    shortVariant = gsub(" ", "", shortVariant, fixed = TRUE)
    return(shortVariant)
  }

  createVCFshortVariant<-function(selectedVariants,variants)
  {
    shortVariant = getFIX(variants[selectedVariants])[,c(1,2,2)]
    shortVariant1 = apply(shortVariant[,1:2],1,function(x) paste(x,collapse=":"))
    shortVariant = apply(cbind(shortVariant1,shortVariant[,2]),1,function(x) paste(x,collapse="-"))
    shortVariant = gsub(" ", "", shortVariant, fixed = TRUE)
    return(shortVariant)
  }

  createVCFlongVariant<-function(selectedVariants,variants)
  {
    longVariant = getFIX(variants[selectedVariants])[,c(1,2,4,5)]
    longVariant1 = apply(longVariant[,1:2],1,function(x) paste(x,collapse=":"))
    longVariant = apply(cbind(longVariant1,longVariant[,c(3,4)]),1,function(x) paste(x,collapse="-"))
    longVariant = gsub(" ", "", longVariant, fixed = TRUE)
    return(longVariant)
  }

  ClinVarConstructor<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    clinvarPOS = match("CLIN_SIG",meta)

    options(warn=-1)
    clinvar = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,clinvarPOS]
    options(warn=0)

    return(clinvar)
  }

  CCRconstructor<-function(VariantAnnotation,selectedVariants,variants)
  {
    consequence<-getConsequence(VariantAnnotation,variants)
    posCCR<-grep("missense_variant|inframe",consequence,perl=TRUE)
    shortVariantCCR<-createVCFshortVariantCCR(selectedVariants,variants)

    outputCCR<-matrix(".",ncol = 1,nrow = length(selectedVariants))

    if(length(posCCR)>0)
    {
      resCCR<-sapply(shortVariantCCR[posCCR],function(x) seqminer::tabix.read.table(CCRautosomesFile, x))
      if(length(resCCR)==length(posCCR))
      {
        posSelectedVariantCCR<-c(1:length(resCCR))[sapply(c(1:length(resCCR)),function(x) length(resCCR[[x]]))>0]
        TabixCCRscore<-unlist(sapply(posSelectedVariantCCR,function(x) mean(resCCR[[x]][[4]])))
        #TabixCCRscore<-resCCR[4,]
        outputCCR[posCCR][posSelectedVariantCCR]<-round(as.numeric(TabixCCRscore), digits = 4)
      }
      if(length(resCCR)>length(posCCR))
      {
        #TabixCCRscore<-resCCR[4,]
        TabixCCRscore<-unlist(sapply(resCCR[4,],function(x) max(x)))
        outputCCR[posCCR]<-round(as.numeric(TabixCCRscore), digits = 4)
      }

      if(length(grep("X",shortVariantCCR[posCCR]))>0)
      {
        resCCR<-sapply(shortVariantCCR[posCCR],function(x) seqminer::tabix.read.table(CCRxchromFile, x))

        if(length(resCCR)==length(posCCR))
        {
          posSelectedVariantCCR<-c(1:length(resCCR))[sapply(c(1:length(resCCR)),function(x) length(resCCR[[x]]))>0]
          TabixCCRscore<-unlist(sapply(posSelectedVariantCCR,function(x) mean(resCCR[[x]][[4]])))
          #TabixCCRscore<-resCCR[4,]
          outputCCR[posCCR][posSelectedVariantCCR]<-round(as.numeric(TabixCCRscore), digits = 4)
        }
        if(length(resCCR)>length(posCCR))
        {
          TabixCCRscore<-unlist(sapply(resCCR[4,],function(x) max(x)))
          outputCCR[posCCR]<-round(as.numeric(TabixCCRscore), digits = 4)

        }
      }
    }
    return(outputCCR)
  }

  CADDconstructor<-function(VariantAnnotation,variants)
  {
    #outputCADD<-matrix(".",ncol = 1,nrow = length(shortVariant))
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    CADD_PHRED = match("CADD_PHRED",meta)

    options(warn=-1)
    outputCADD = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,CADD_PHRED]
    options(warn=0)


    return(outputCADD)
  }

  gnomADConstructor<-function(VariantAnnotation,variants)
  {
    #meta = queryMETA(variants,"CSQ")[[3]][4]
    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    gnomADPOS = match("AF",meta)
    gnomADgenPOS = match("gnomADg_AF",meta)
    gnomADexoPOS = match("gnomADe_AF",meta)


    process <- function(x) {
      MAF = strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation[x])), "[|]")[[1]][c(gnomADPOS)]
      gen = strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation[x])), "[|]")[[1]][c(gnomADgenPOS)]
      exome = strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation[x])), "[|]")[[1]][c(gnomADexoPOS)]
      c(MAF,gen,exome)[order(c(MAF,gen,exome),decreasing = TRUE)][1]
    }

    gnomAD <-do.call(c, lapply(c(1:length(VariantAnnotation)), process))


    return(gnomAD)
  }


  gnomADexomeConstructor<-function(shortVariant,longVariant)
  {
    #gnomadExome
    gnomADexome<-matrix(".",ncol = 1,nrow = length(shortVariant))

    gnomADTabix <- seqminer::tabix.read.table(exomeFile, shortVariant)


    if(length(gnomADTabix)>0)
    {
      options(warn=-1)
      gnomadlongVariant = do.call(paste, c(gnomADTabix[,c(1,2,4,5)], sep="-"))
      options(warn=0)
      gnomadlongVariant = gsub("^(\\d+)-(\\d+)-(\\w+)-(\\w+)", "\\1:\\2-\\3-\\4", gnomadlongVariant, perl = TRUE)

      #gnomadlongVariant = apply(gnomADTabix[,1:2],1,function(x) paste(x,collapse=":"))
      #gnomadlongVariant = apply(cbind(gnomadlongVariant,gnomADTabix[,c(4,5)]),1,function(x) paste(x,collapse="-"))
      #gnomadlongVariant = gsub(" ", "", gnomadlongVariant, fixed = TRUE)

      posGnomadTabix<-c(1:length(gnomadlongVariant))[!is.na(match(gnomadlongVariant,longVariant))]
      posSelectedVariants<-match(gnomadlongVariant,longVariant)[!is.na(match(gnomadlongVariant,longVariant))]

      m<-regexpr("controls_AF_popmax=\\S+?;", gnomADTabix[posGnomadTabix,8],perl=TRUE)
      match<-regmatches(gnomADTabix[posGnomadTabix,8], m)
      match<-gsub("controls_AF_popmax=", "",match,perl=TRUE)
      matchGnomAD<-gsub("\\;", "",match,perl=TRUE)

      gnomADexome[posSelectedVariants]<-matchGnomAD
    }
    return(gnomADexome)
  }

  getMetrics<-function(VariantAnnotation,variants)
  {
    options(warn=-1)
    genesList = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",VariantAnnotation)), "[|]"))[,4]
    options(warn=0)

    genes<-getGeneName(VariantAnnotation,variants)
    pos<-match(genes,geneConstraints$gene)
    geneMetrics<-geneConstraints[pos,2:6]
    #genesList<-as.matrix(geneConstraints[pos,1])
    #genesList[is.na(pos),1]<-genes[is.na(pos)]
    #genesList<-total_unique[match(genesList,total_unique[,2]),1]

    geneMetrics<-cbind(genesList,geneMetrics)
    return(geneMetrics)
  }

  getPhenoScore<-function(VariantAnnotation,GlobalPhenotypicScore,variants)
  {
    genes<-getGeneName(VariantAnnotation,variants)

    FunctPos<-match(genes,GlobalPhenotypicScore[,2])
    PhysPos<-match(genes,GlobalPhenotypicScore[,5])

    FunctScore<-GlobalPhenotypicScore[FunctPos,3]
    PhysScore<-GlobalPhenotypicScore[PhysPos,6]

    FunctScore[is.na(FunctScore)]<-0.0001
    PhysScore[is.na(PhysScore)]<-0.0001

    FunctPos[is.na(FunctPos)]<-dim(GlobalPhenotypicScore)[1]/2
    PhysPos[is.na(PhysPos)]<-dim(GlobalPhenotypicScore)[1]/2

    FunctScore<-round(as.numeric(FunctScore),digits = 4)
    PhysScore<-round(as.numeric(PhysScore),digits = 4)

    PhenoScore<-cbind(FunctPos,FunctScore,PhysPos,PhysScore)
    return(PhenoScore)
  }

  getKnownDisease<-function(VariantAnnotation,variants)
  {
    genes<-c("")
    genes<-c(genes,getGeneName(VariantAnnotation,variants))
    #indices<-sapply(unique(genes), function(x) which(HPO2genes[,1] == x))
    genes<-genes[-1]
    #indices<-indices[-1]
    HPOrecessive<-"HP:0000007|HP:0001416|HP:0001526"
    HPOdominant<-"HP:0000006|HP:0001415|HP:0001447|HP:0001448|HP:0001451|HP:0001455|HP:0001456|HP:0001463"
    HPOXlink<-"HP:0001417|HP:0001418|HP:0001419|HP:0001423"

    #recessive<-unlist(sapply(c(1:length(indices)), function(x) grep(HPOrecessive,HPO2genes[indices[[x]],2])[1]))
    #dominant<-unlist(sapply(c(1:length(indices)), function(x) grep(HPOdominant,HPO2genes[indices[[x]],2])[1]))
    #xlink<-unlist(sapply(c(1:length(indices)), function(x) grep(HPOXlink,HPO2genes[indices[[x]],2])[1]))
    recessive<-match(genes,HPO2genes[grep(HPOrecessive,HPO2genes[,2]),1])
    dominant<-match(genes,HPO2genes[grep(HPOdominant,HPO2genes[,2]),1])
    xlink<-match(genes,HPO2genes[grep(HPOXlink,HPO2genes[,2]),1])

    knownDisease <- matrix("",nrow = length(genes), ncol = 1)
    knownDisease[!is.na(recessive)]<-paste(knownDisease[!is.na(recessive)],"|AR|",sep="")
    knownDisease[!is.na(dominant)]<-paste(knownDisease[!is.na(dominant)],"|AD|",sep="")
    knownDisease[!is.na(xlink)]<-paste(knownDisease[!is.na(xlink)],"|XLINK|",sep="")
    return(knownDisease)
  }

  getClinPriorScore<-function(VariantScore,phenoScore)
  {
    #ClinPriorScore<-(0.2*as.numeric(VariantScore)) + (10*(as.numeric(phenoScore[,2])/as.numeric(phenoScore[,1])))
    #ClinPriorScore<-(0.2*as.numeric(VariantScore)) + 0.8*(log10(1 + exp(-800*((as.numeric(phenoScore[,1])/20146)-0.2)))*(as.numeric(phenoScore[,2])))
    position = do.call(pmin,as.data.frame(phenoScore[,c(1,3)]))
    score = do.call(pmax,as.data.frame(phenoScore[,c(2,4)]))
    #ClinPriorScore<-0.9*(log10(1 + exp(-80*((as.numeric(as.numeric(min(phenoScore[,c(1,3)])/dim(total_unique)[1])-0.2)))*(as.numeric(max(phenoScore[,c(2,4)])))))+0.1*as.numeric(VariantScore))
    ClinPriorScore<-10*(log10(1 + exp(-30*(as.numeric(position/dim(total_unique)[1])-0.3)))*score)*as.numeric(VariantScore)
    ClinPriorScore[VariantScore==0]<-ClinPriorScore[VariantScore==0]/20
    return(round(ClinPriorScore,digits = 4))
  }

  filterGnomAD<-function(gnomAD,freqThreshold)
  {
    posNa = c(1:length(gnomAD))[is.na(as.numeric(gnomAD))]
    posThreshold = c(1:length(gnomAD))[!is.na(as.numeric(gnomAD))][as.numeric(gnomAD)[!is.na(as.numeric(gnomAD))]<=as.numeric(freqThreshold)]
    selectedVariants = sort(union(posNa,posThreshold))
    return(selectedVariants)
  }

  orderTable<-function(tableRaw, inheritance)
  {
    tableOrder<-tableRaw
    if(inheritance == "HOM" || inheritance == "HET" || inheritance == "XLink")
    {
      tableRaw<-tableRaw[order(tableRaw$ClinPriorScore, decreasing=TRUE),]
      #inheritanceRaw<-sapply(c(1:length(tableRaw$ClinPriorScore)),function(x) paste(x, inheritance,sep = "-"))
      remove<-tableRaw$ClinPriorScore
      #tableOrder<-cbind(tableRaw,inheritanceRaw,remove)
      tableOrder<-cbind(tableRaw,inheritance,remove)
    }

    if(inheritance == "CmpHET")
    {
      geneUnique<-unique(tableRaw$genesList[duplicated(tableRaw$genesList)])

      process <- function(x) {
        posGU =  grep(paste("\\b",geneUnique[x],"\\b",sep=""),tableRaw$genesList)
        score = matrix(mean(tableRaw$ClinPriorScore[posGU][order(tableRaw$ClinPriorScore[posGU],decreasing=TRUE)][1:2]),ncol=1,nrow = length(posGU))
        posGU
      }

      pos <-do.call(c, lapply(c(1:length(geneUnique)), process))

      process <- function(x) {
        posGU =  grep(paste("\\b",geneUnique[x],"\\b",sep=""),tableRaw$genesList)
        score = matrix(mean(tableRaw$ClinPriorScore[posGU][order(tableRaw$ClinPriorScore[posGU],decreasing=TRUE)][1:2]),ncol=1,nrow = length(posGU))
        score
      }

      remove <-do.call(c, lapply(c(1:length(geneUnique)), process))

      tableOrder<-cbind(tableRaw[pos,],inheritance,remove)
      tableOrder = tableOrder[order(tableOrder$remove,decreasing = TRUE),]


    }
    return(tableOrder)
  }

  constructTable<-function(selectedVariants,inheritance, freqThreshold, isCodingVar,variants)
  {
    Ordertable<-c()
    gt = variants@gt
    pos = match(sampleName, colnames(gt))


    VariantAnnotation<-getVariant(selectedVariants,variants)
    if(isCodingVar){selectedVariants2 = isCoding(VariantAnnotation,variants)
    selectedVariants = selectedVariants[selectedVariants2]
    }
    VariantAnnotation<-getVariant(selectedVariants,variants)
    gnomAD<-gnomADConstructor(VariantAnnotation,variants)
    selectedVariants2 = filterGnomAD(gnomAD,freqThreshold)
    selectedVariants = selectedVariants[selectedVariants2]

    if(length(selectedVariants)>0)
    {
        selectedVariants = c(1,2,selectedVariants)
        VariantAnnotation<-getVariant(selectedVariants,variants)
        gnomAD<-gnomADConstructor(VariantAnnotation,variants)
        Moderate<-isModerate(VariantAnnotation,variants)
        HIGH<-isHigh(VariantAnnotation,variants)
        Splicing<-isSplicing(VariantAnnotation,variants)
        metrics<-getMetrics(VariantAnnotation,variants)
        cDNA<-getcDNA(VariantAnnotation,variants)
        Protein<-getProtein(VariantAnnotation,variants)
        Consequence<-getConsequence(VariantAnnotation,variants)
        clinvar<-ClinVarConstructor(VariantAnnotation,variants)
        CADD<-CADDconstructor(VariantAnnotation,variants)
        CCR<-CCRconstructor(VariantAnnotation,selectedVariants,variants)
        GTsamples = getGTsamples(selectedVariants,variants)
        VariantScore<-sapply(c(1:length(selectedVariants)),function(x) calculateVariantScore(HIGH[x], Splicing[x], Moderate[x], CADD[x], CCR[x], clinvar[x], metrics[x,], inheritance,GTsamples[x,pos]))
        phenoScore<-getPhenoScore(VariantAnnotation,GlobalPhenotypicScore,variants)
        ClinPriorScore<-getClinPriorScore(VariantScore,phenoScore)

        knownDisease<-getKnownDisease(VariantAnnotation,variants)

        phaseGT = getPhysicalPhasing(selectedVariants,variants)
        coding = getisCoding(VariantAnnotation,variants)
        GTsamplesFull = gt[selectedVariants,pos]
        table<-cbind(getFIX(variants[selectedVariants])[,c(1:5)],metrics,clinvar,knownDisease,coding,Consequence,cDNA,Protein,HIGH,Splicing,Moderate,CADD,CCR,gnomAD,phenoScore,VariantScore,ClinPriorScore,GTsamples,GTsamplesFull,phaseGT)
        Ordertable<-orderTable(table, inheritance)

        #delete extra
        posEXTRA = match(variants@fix[selectedVariants[c(1:2)],2],Ordertable$POS)[!is.na(match(variants@fix[selectedVariants[c(1:2)],2],Ordertable$POS))]
        if(length(posEXTRA)>0){Ordertable = Ordertable[-posEXTRA,]}

    }
    return(Ordertable)
  }


#variants<-filterVariantsVcfR(variants,filter)
#Blacklist
 # blacklist = read.csv(blacklistFile, sep="\t", header=FALSE)
 # ensemble = do.call(paste, c(blacklist, sep="-"))

  #delete ensemble exome,genome Gnomad and blacklist variants
  #ensemble<-unique(c(shortSNPexome,shortSNPgenome,blacklist))
  # deleteVariants<-deleteGnomadBlacklistVariants()
  #if(length(deleteVariants)>0){variants<-variants[-deleteVariants,]}

  #variants<-getAlt_forMutliAllelic(variants)

  #delete indels in same genomic position
  #deleteVariants<-deleteIndelVariantsSamePos()
  #if(length(deleteVariants)>0){variants<-variants[-deleteVariants,]}



#HOM variants
HOMvariants<-getHOMvariants(variants)
tableHOM<-c()
if(length(HOMvariants)>0){tableHOM<-constructTable(HOMvariants,"HOM", frequenceAR,isCodingVar,variants)}

#AD variants
ADvariants<-getADvariants(variants)
tableAD<-c()
if(length(ADvariants)>0){tableAD<-constructTable(ADvariants,"HET", frequenceAD,isCodingVar,variants)}

#Compound heterozygous variants
CMPHETvariants<-getCmpHETvariants(variants)
tableCMPHET<-c()
if(length(CMPHETvariants)>0){tableCMPHET<-constructTable(CMPHETvariants,"CmpHET", frequenceAR,isCodingVar,variants)}

#XLINK variants
XLvariants<-getXLvariants(variants)
tableXL<-c()
if(length(XLvariants)>0){tableXL<-constructTable(XLvariants,"XLink",frequenceAR,isCodingVar,variants)}

tableEnsemble<-rbind(tableHOM,tableAD,tableCMPHET,tableXL)
tableEnsemble<-tableEnsemble[order(tableEnsemble$remove,tableEnsemble$genesList,decreasing=TRUE),]
del<-match("remove",colnames(tableEnsemble))
tableEnsemble<-tableEnsemble[,-del]
ClinPriorPosition <- c(1:dim(tableEnsemble)[1])
tableEnsemble<-cbind(ClinPriorPosition,tableEnsemble)

gc()
return(tableEnsemble)
}
