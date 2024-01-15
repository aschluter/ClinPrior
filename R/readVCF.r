
#' Filtering the patient variants in VCF format annotated by VEP
#'
#' @param sampleName character name with the patient code as is written in the VCF file.
#' @param filter character name with the desired filter to select the variants. It should be as is written in the FILTER column in the vcf file. 
#' @param geneQuality numeric value indicating the desired GQ threshold. Default 20.
#' @param readDepth numeric value indicating the desired DP threshold. Default 10.
#' @param variants object of class vcfR-class with the patient variants using the read.vcfR function from vcfR package.
#' @param assembly Genome assembly used. Default assembly human GRCh37.
#' @param distSplicThreshold integer indicating the maximum distance in base pairs (bp) allowed between intronic variants and the intron-exon boundary. Default 1000.
#' @param synonymous logical indicating whether to include the synonymous variants. Default TRUE.
#' 
#' @return returns an object of class vcfR-class. 
#' @export
#'
#' @examples
#' library(vcfR)
#' vcfFile = paste(system.file("extdata/example", package = "ClinPrior"),"HG001_GRCh37_1_22_v4.2.1_benchmark.vep01.KCNQ2Met546Thr.vcf.gz",sep="/")
#' variants <- read.vcfR(vcfFile)
#' variantsFiltered <- readVCF(sampleName = "HG001",variants=variants)
readVCF<-function(sampleName,filter = "",geneQuality = 20,readDepth = 10,variants,assembly = "assembly37",distSplicThreshold = 1000,synonymous = TRUE)
{
  library(vcfR)

  dest = system.file(paste("extdata",assembly,sep = "/"), package = "ClinPrior")
  ClinPriorfiles<-list.files(path = dest, full.names = TRUE)
  blacklistFile = ClinPriorfiles[grep("blacklist",ClinPriorfiles)]
  #functions

  getADvariants<-function(variants)
  {
    gt = variants@gt
    pos = match(sampleName, colnames(gt))
    posADvariants<-grep("0/1|0/2|0/3|1/2|2/3|1/3",gt[,pos])
    posADvariants<-posADvariants[grep("\\bX\\b",variants@fix[posADvariants,1],invert=TRUE)]
    return(posADvariants)
  }

  filterVariantsVcfR<-function(variants,filter)
  {
    if(filter != "")
    {
      filterColumn<-variants@fix[,7]
      pos = grep(filter,filterColumn)
      if(length(pos)>0){variants<-variants[pos,]}
    }
    return(variants)
  }

  deleteGnomadBlacklistVariants<-function(variants,ensemble)
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

  qualVariants<-function(geneQuality,readDepth,variants)
  {

    gt = variants@gt
    pos = match(sampleName, colnames(gt))

    #GQ
    posGQDiscard = integer(0)

    lGQ = length(grep("GQ",as.character(variants@gt[,1]),perl = TRUE))

    if(length(lGQ)>0)
    {

      process <- function(x) {
        posGQ <-  match("GQ",strsplit(as.character(variants@gt[x,1]), ":")[[1]])
        GQ <- strsplit(as.character(variants@gt[x,pos]), ":")[[1]][posGQ]
        GQ
      }

      GQ <-do.call(rbind, lapply(c(1:dim(variants)[1]), process))

      posGQDiscard1<-  c(1:length(GQ))[!is.na(as.numeric(GQ)<=geneQuality)]
      posGQDiscard2<-  c(1:length(posGQDiscard1))[as.numeric(GQ)[posGQDiscard1]<as.numeric(geneQuality)]
      posGQDiscard<- c(1:length(GQ))[posGQDiscard1][posGQDiscard2]
    }


    #AD
    posDeepDiscard = integer(0)

    lDP = length(grep("DP",as.character(variants@gt[,1]),perl = TRUE))

    if(length(lDP)>0)
    {

      process <- function(x) {
        posDP <-  grep("^DP$|^DPI$",perl=TRUE,strsplit(as.character(variants@gt[x,1]), ":")[[1]])[1]
        deep <- strsplit(as.character(variants@gt[x,pos]), ":")[[1]][posDP]
        deep
      }


      deep <-do.call(rbind, lapply(c(1:dim(variants)[1]), process))

        posDPDiscard1<-  c(1:length(deep))[!is.na(as.numeric(deep)<=readDepth)]
        posDPDiscard2<-  c(1:length(posDPDiscard1))[as.numeric(deep)[posDPDiscard1]<as.numeric(readDepth)]
        posDeepDiscard<- c(1:length(deep))[posDPDiscard1][posDPDiscard2]
     }


    #Discard bias <25% and >75% in HET variants

    variants2<- variants
    posVariants2 = c(1:as.numeric(dim(variants)[1]))

    posDiscard<-unique(c(posGQDiscard,posDeepDiscard))

    if(length(posDiscard)>0)
    {

      posVariants2 = c(1:as.numeric(dim(variants)[1]))[-posDiscard]
      variants2<- variants[-posDiscard,]
    }

    posBiasDiscard = integer(0)

    posADvariants<-getADvariants(variants2)

    #GT
    lGT = length(grep("GT",as.character(variants2@gt[,1]),perl = TRUE))

    if(length(lGT)>0)
    {

      process <- function(x) {
        posGT <-  match("GT",strsplit(as.character(variants2@gt[x,1]), ":")[[1]])
        GT <- strsplit(as.character(variants2@gt[x,pos]), ":")[[1]][posGT]
        genotype1 = strsplit(GT,"/")[[1]][1]
        genotype2 = strsplit(GT,"/")[[1]][2]
        list(genotype1,genotype2)
      }

      genotype <-do.call(rbind, lapply(c(1:dim(variants2)[1]), process))



    #AD
    lAD = length(grep("AD",as.character(variants2@gt[,1]),perl = TRUE))

    if(length(lAD)>0)
    {

     HETvariants<-variants2@gt[posADvariants,]
          genotype1<- unlist(genotype[posADvariants,1])
          genotype2<- unlist(genotype[posADvariants,2])

         process <- function(x) {
            numG1 = as.numeric(genotype1[x])
            numG2 = as.numeric(genotype2[x])

            posAD <-  match("AD",strsplit(as.character(HETvariants[x,1]), ":")[[1]])
            deepAD <- strsplit(as.character(HETvariants[x,pos]), ":")[[1]][posAD]
            allel1 = as.numeric(strsplit(deepAD,",")[[1]][(numG1+1)])
            allel2 = as.numeric(strsplit(deepAD,",")[[1]][(numG2+1)])
           freqBias = allel1/(allel1+allel2)
          }

          freqBias <-do.call(rbind, lapply(c(1:length(posADvariants)), process))

          posnoNA<-c(1:length(freqBias))[!is.na(freqBias)]
          freqBias<-freqBias[posnoNA]
          posDiscard2<-union(c(1:length(freqBias))[freqBias>0.75],c(1:length(freqBias))[freqBias<0.25])
          posBiasDiscard<-c(1:length(deep))[posADvariants][posnoNA][posDiscard2]
          posBiasDiscard<- posVariants2[posBiasDiscard]
      }

    }

    posDiscard<-unique(c(posGQDiscard,posDeepDiscard,posBiasDiscard))
    #posDiscard<-unique(c(posGQDiscard,posDeepDiscard,posBiasDiscard1,posBiasDiscard2,posBiasDiscard3))

    #if(qualityThreshold>0)
    #{
    #  pos100<-c(1:length(variants[,6]))[as.numeric(gsub(",",".",variants[,6]))<qualityThreshold]
    #  posDiscard<-union(posDiscard,pos100[!is.na(pos100)])
    #}

    return(posDiscard)
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

  getIMPACT<-function(VariantAnnotation,variants)
  {

    meta = queryMETA(variants,"\\bCSQ\\b")[[1]][4]
    meta = unlist(strsplit(meta, "[|]"))
    impactPOS = match("IMPACT",meta)

    options(warn=-1)
    impact = do.call(rbind, strsplit(as.character(gsub(".*CSQ=","",variants@fix[,8])), "[|]"))[,impactPOS]
    options(warn=0)

    return(impact)
  }

  ####
  #remove Chr                       
  if(length(grep("^chr",variants@fix[,1]))>0){variants@fix[,1]<-gsub("^chr","",variants@fix[,1])}

                         
  #ObjectVcfR<-read.vcfR(vcfFile)
  gt = variants@gt
  pos = match(sampleName, colnames(gt))
  gtPOS = grep("0/0",gt[,pos],invert=TRUE)
  variants = variants[gtPOS,]



  consequence = getConsequence(variants@fix[,8],variants)
  posCONSE<-grep("intron",consequence)


  cDNA = getcDNA(variants@fix[,8],variants)
  m <- regexpr("\\+[\\d+]{1,8}|\\-[\\d+]{1,8}", cDNA,perl=TRUE)
  posMatch <- grep("\\+[\\d+]{1,8}|\\-[\\d+]{1,8}", cDNA,perl=TRUE)
  match<-regmatches(cDNA, m)
  distSplic<-as.numeric(gsub("\\+|\\-", "",match,perl=TRUE))

  posSplicThreshold<-c(1:length(distSplic))[distSplic<=distSplicThreshold]
  pos30<-c(1:length(cDNA))[posMatch][posSplicThreshold]
  posKEEP = intersect(pos30,posCONSE)


  impact = getIMPACT(variants@fix[,8],variants)
  posIMPACT<-grep("HIGH|MODERATE|LOW",impact)

  posSYNOM<-grep("synonymous",consequence)
  posSPLICE<-grep("splice",consequence)

  posSYNOM_SPLICE = intersect(posSYNOM,posSPLICE)

  posSYNOM_NO_SPLICE = posSYNOM[is.na(match(posSYNOM,posSYNOM_SPLICE))]


  posKEEP = sort(union(posIMPACT,posKEEP))
  if(synonymous==FALSE){posKEEP = setdiff(posKEEP,posSYNOM_NO_SPLICE)}

  if(length(posKEEP)>0){variants<-variants[posKEEP,]}



  posDISCARD<-qualVariants(geneQuality,readDepth,variants)
  if(length(posDISCARD)>0){variants<-variants[-posDISCARD,]}


  variants<-filterVariantsVcfR(variants,filter)

  #delete ensemble exome,genome Gnomad and blacklist variants
  #Blacklist
  blacklist = read.csv(blacklistFile, sep="\t", header=FALSE)
  ensemble = do.call(paste, c(blacklist, sep="-"))
  deleteVariants<-deleteGnomadBlacklistVariants(variants,ensemble)
  if(length(deleteVariants)>0){variants<-variants[-deleteVariants,]}

  variants<-getAlt_forMutliAllelic(variants)

  #delete indels in same genomic position
  deleteVariants<-deleteIndelVariantsSamePos(variants)
  if(length(deleteVariants)>0){variants<-variants[-deleteVariants,]}


  return(variants)
}
