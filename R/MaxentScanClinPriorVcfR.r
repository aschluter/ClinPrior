#' Title
#'
#' @param Threshold
#' @param assembly
#' @param matrixVariants
#'
#' @return
#' @export
#'
#' @examples
MaxentScanClinPriorVcfR<-function(matrixVariants, assembly="assembly37",Threshold=0.2)
{

  # matrixVariants
  #CHROM = c(22,22,22,22)
  #POS = c(21064114,21083617,21119035,21152931)
  #REF = c("T","C","C","A")
  #ALT = c("C","T","T","G")
  #genesList = c("PI4KA","PI4KA","PI4KA","PI4KA")
  #matrixVariants = data.frame(cbind(CHROM,POS,REF,ALT,genesList))

  if(assembly == "assembly37"){
    library(BSgenome.Hsapiens.UCSC.hg19)
    BSgenome_assembly = BSgenome.Hsapiens.UCSC.hg19}
  if(assembly == "assembly38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    #BSgenome_assembly = deparse(substitute(BSgenome.Hsapiens.UCSC.hg38))}
    BSgenome_assembly = BSgenome.Hsapiens.UCSC.hg38}
  library(reticulate)


  #currentDir<-getwd()
  maxentpyPath<-system.file("extdata", package = "ClinPrior")
  #setwd(maxentpyPath)
  maxentpy<-import_from_path("maxentpy", path = maxentpyPath, convert = TRUE)
  maxentpy$maxent$dir_path<-maxentpyPath
  #setwd(currentDir)



  getStrand<-function(gene)
  {
    strand<-gene2strand[match(gene,gene2strand[,1]),2]
    return(strand)
  }

  Ref2Alt3<-function(seq1,seq2,change1,change2)
  {
    scores<-0
    m<-gregexpr("[A|C|G|T]{18}AG[A|C|G|T]{3}",as.character(seq1),perl=TRUE)
    seqs<-regmatches(as.character(seq1), m)[[1]]


    foreachseq<-function(seqOne)
    {
      score<-0
        m<-gregexpr(seqOne,as.character(seq1),perl=TRUE)
        l<-m[[1]][1]
        #lpost<-l+as.numeric(attributes(m[[1]])[1])

          if(l>3 && l<(27+nchar(change1)))
          {
            #cmd<-paste("perl score3seq.pl", seqOne, sep=" ")
            #res<-system(cmd,intern = TRUE)
            #score<-as.numeric(strsplit(res,"\t")[[1]][2])
            score<-maxentpy$maxent$score3(seqOne)

            if(score>=3)
            {
                #difference<-nchar(change2)-nchar(change1)
                altSeq<-seq2[l:(l+22)]
                #if(l>26){altSeq<-seq2[27:(27+22)]}
                #cmd<-paste("perl score3seq.pl", altSeq, sep=" ")
                #res<-system(cmd,intern = TRUE)
                score2<-maxentpy$maxent$score3(as.character(altSeq))
                score<-if(score2<score){abs(score-score2)/score}else{score =0}
                if(score<Threshold){score = 0}
                #score<-score/as.numeric(strsplit(res,"\t")[[1]][2])
            }else{score =0}
          }
      return(score)
    }

    scores<-if(length(seqs)>0){lapply(seqs,function(x) foreachseq(x))}
    if(length(scores)==0){scores=0}

  return(scores)
  }

  Ref2Alt5<-function(seq1,seq2,change1,change2)
  {
    scores<-0
    m<-gregexpr("[A|C|G|T]{3}GT[A|C|G|T]{4}",as.character(seq1),perl=TRUE)
    seqs<-regmatches(as.character(seq1), m)[[1]]

    foreachseq<-function(seqOne)
    {
      score<-0
      m<-gregexpr(seqOne,as.character(seq1),perl=TRUE)
      l<-m[[1]][1]
      #lpost<-l+as.numeric(attributes(m[[1]])[1])

      if(l>18 && l<(27+nchar(change1)))
      {
        #cmd<-paste("perl score5seq.pl", seqOne, sep=" ")
        #res<-system(cmd,intern = TRUE)
        #score<-as.numeric(strsplit(res,"\t")[[1]][2])
        score<-maxentpy$maxent$score5(seqOne)

        if(score>=3)
        {
          #difference<-nchar(change2)-nchar(change1)
          altSeq<-seq2[l:(l+8)]
          #if(l>26){altSeq<-seq2[27:(27+8)]}
          #cmd<-paste("perl score5seq.pl", altSeq, sep=" ")
          #res<-system(cmd,intern = TRUE)
          score2<-maxentpy$maxent$score5(as.character(altSeq))
          score<-if(score2<score){abs(score-score2)/score}else{score =0}
          if(score<Threshold){score = 0}
          #score<-score/as.numeric(strsplit(res,"\t")[[1]][2])
        }else{score=0}
      }
      return(score)
    }

    scores<-if(length(seqs)>0){lapply(seqs,function(x) foreachseq(x))}
    if(length(scores)==0){scores=0}
    return(scores)
  }

  ##################################

  MaxentScan<-function(CHROMpos)
  {

    chr = paste("chr",CHROMpos$CHROM,sep="")
    start = CHROMpos$POS
    Ref = as.character(CHROMpos$REF)
    Alt = as.character(CHROMpos$ALT)
    query = as.character(CHROMpos$genesList)

    gene<-total_unique[match(query,total_unique[,1]),2]

    SplicingScores<-c(0,0,0,0)
    pos1 = as.numeric(start) -26
    pos2 = as.numeric(start) + 25 + nchar(Ref)

    if(!is.na(getStrand(gene)))
    {
          seq<-getSeq(BSgenome_assembly, chr,	start = as.integer(pos1),end= as.integer(pos2))
          #############
          ##ALt to Ref#
          ############
          pos1 = as.integer(start)+nchar(Ref)
          pos2 = pos1+26
          seq2<-getSeq(BSgenome_assembly, chr,	start = as.integer(pos1),end= as.integer(pos2))
          seqALT<-seq
          seqALT<-paste(seqALT[1:26],Alt,seq2,sep="")
          seqALT<-DNAString(seqALT)


          if(getStrand(gene)=="+")
          {

            ##Ref to ALt
            scoresRef2Alt3<-Ref2Alt3(seq,seqALT,Alt,Ref)
            scoresRef2Alt5<-Ref2Alt5(seq,seqALT,Alt,Ref)

            ##ALt to Ref
            scoresAlt2Ref3<-Ref2Alt3(seqALT,seq,Ref,Alt)
            scoresAlt2Ref5<-Ref2Alt5(seqALT,seq,Ref,Alt)

            SplicingScores[1]<-(-unlist(scoresRef2Alt3)[order(unlist(scoresRef2Alt3),decreasing = TRUE)[1]])
            SplicingScores[3]<-(-unlist(scoresRef2Alt5)[order(unlist(scoresRef2Alt5),decreasing = TRUE)[1]])
            SplicingScores[2]<-unlist(scoresAlt2Ref3)[order(unlist(scoresAlt2Ref3),decreasing = TRUE)[1]]
            SplicingScores[4]<-unlist(scoresAlt2Ref5)[order(unlist(scoresAlt2Ref5),decreasing = TRUE)[1]]

          }

          if(getStrand(gene)=="-")
          {
              seqComplement<-reverseComplement(seq)
              seqALTComplement<-reverseComplement(seqALT)

              RefCompl<-as.character(reverseComplement(DNAString(Ref)))
              AltCompl<-as.character(reverseComplement(DNAString(Alt)))



              ##Ref to ALt complementary
              scoresRef2Alt3_AS<-Ref2Alt3(seqComplement,seqALTComplement,AltCompl,RefCompl)
              scoresRef2Alt5_AS<-Ref2Alt5(seqComplement,seqALTComplement,AltCompl,RefCompl)

              ##ALt to Ref complementary
              scoresAlt2Ref3_AS<-Ref2Alt3(seqALTComplement,seqComplement,RefCompl,AltCompl)
              scoresAlt2Ref5_AS<-Ref2Alt5(seqALTComplement,seqComplement,RefCompl,AltCompl)

              SplicingScores[1]<-(-unlist(scoresRef2Alt3_AS)[order(unlist(scoresRef2Alt3_AS),decreasing = TRUE)[1]])
              SplicingScores[3]<-(-unlist(scoresRef2Alt5_AS)[order(unlist(scoresRef2Alt5_AS),decreasing = TRUE)[1]])
              SplicingScores[2]<-unlist(scoresAlt2Ref3_AS)[order(unlist(scoresAlt2Ref3_AS),decreasing = TRUE)[1]]
              SplicingScores[4]<-unlist(scoresAlt2Ref5_AS)[order(unlist(scoresAlt2Ref5_AS),decreasing = TRUE)[1]]

          }
    }
    gc()
    SplicingScores = round(SplicingScores,digits = 2)
    return(matrix(SplicingScores,nrow = 1,ncol = 4))
  }

  acumulat = 0
  for(x in 1:dim(matrixVariants)[1])
  {
    print(x)
    acumulat = rbind(acumulat,cbind(matrixVariants[x,],MaxentScan(matrixVariants[x,])))
  }

  maxentScanOutput = acumulat[-1,]

  #branchpoint
  #longVariant =  cbind.data.frame(do.call(paste, c(matrixVariants[,1:2], sep=":")),matrixVariants[,3:4])
  #longVariant = do.call(paste, c(longVariant[,1:3], sep="-"))
  #longVariant = gsub(" ", "", longVariant, fixed = TRUE)


  #branchpoint<-ClinPrior::branchpointer(longVariant,assembly,0.2)
  #maxentScanOutput<-cbind(maxentScanOutput,branchpoint)

return(maxentScanOutput)
}
