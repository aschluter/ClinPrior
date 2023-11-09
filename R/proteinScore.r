#' Phenotypic score for each gene based on the description of the patient in HPO terms
#'
#' @param HPOpatient vector with the patient HPO terms
#'
#' @return vector with the gene scores
#' @export
#'
#' @examples
#' HPOpatient = c("HP:0004481","HP:0002376","HP:0001257","HP:0001250","HP:0000238","HP:0002922","HP:0000365")
#' Y<-proteinScore(HPOpatient)
proteinScore<-function(HPOpatient)
{
    library(igraph)
    library(Matrix)

    HPO2genes<-HPO2genes 
    treureHPO<-treureHPO
    HPOqueryGene<-HPOqueryGene
    total_unique<-total_unique

    g1<-as.undirected(g1)
    HPOpatient<-as.matrix(unique(HPOpatient))

    treureHPO <- paste(treureHPO,collapse="|")
    posTreure<-grep(treureHPO, HPO2genes[,2])
    if(length(posTreure)>0) HPO2genes<-HPO2genes[-posTreure,]

    # Delete some irrelevant HPOs
    posTreure<-grep(treureHPO, HPOpatient)
    if(length(posTreure)>0) HPOpatient<-HPOpatient[-posTreure]

    #create a subgraph with the expanded (1 order) patient HPOs
    HPOorig_expanded<-unique(unlist(sapply(HPOpatient, function (x) rownames(as.matrix(igraph::ego(g1, order = 1, x)[[1]])))))
    g.sub <- induced.subgraph(graph = g1, HPOorig_expanded)
    res<-cluster_edge_betweenness(g.sub)
    HPOorig_expanded<-cbind(res$names,res$membership)
    HPOorigGroups <<- HPOorig_expanded[match(HPOpatient,HPOorig_expanded[,1]),]

    genes <- unique(HPO2genes[, 1]) #genes with HPO
    posGenes<-match(genes,rownames(HPOadj))
    acumulat<-Matrix(matrix(0,ncol = 1,nrow =length(rownames(HPOadj)[posGenes])))
    acumulatFreq<-acumulat

    for(z in 1:length(HPOpatient))
    {
      pos<-match(HPOpatient[z],colnames(HPOadj))
      HPOPatientItem = HPOdistance[pos,]
      column<-HPOadj[posGenes,] %*% HPOPatientItem
      acumulat<-cbind(acumulat,column)
    }
    acumulat<-acumulat[,-1]

    if(length(HPOorigGroups[,2])!=length(unique(HPOorigGroups[,2])))
    {
      memo = ""
      dupli =  unique(HPOorigGroups[duplicated(HPOorigGroups[,2]),2])
      for(i in 1:length(dupli))
      {
        pos<-which(HPOorigGroups[,2] == dupli[i])
        new<-Matrix::rowSums(acumulat[,pos])
        acumulat[,pos[1]]<-new
        memo<-c(memo,pos[-1])
      }
      memo <- memo[-1]
      acumulat<-acumulat[,-as.numeric(memo)]

    }



    acumulat = acumulat / acumulat
    acumulat[is.na(acumulat)]=0
    acumulat = Matrix(acumulat)
    HPOmatch_quant<-Matrix::rowSums(acumulat)


    q=HPOmatch_quant-1
    m=length(unique(HPOorigGroups[,2]))
    n=length(V(g1)$name)-length(unique(HPOorigGroups[,2]))
    k=as.matrix(HPOqueryGene)*10
        stats<-phyper(q,m,n,k,lower.tail = FALSE, log.p = FALSE)
    gc()
    stats[stats == -Inf] = 1
    stats[stats == 0] = 10^(log10(min(stats[stats != 0]))-1)
    testResult<-stats

    D = abs(log10(abs(testResult)))
    Dred <- as.numeric(D)
    pos <- match(genes, total_unique[,2])
    D <- matrix(0, nrow = length(total_unique[,2]))
    D[pos] <- Dred
    DNormed <- (D - min(D, na.rm = TRUE))/(max(D, na.rm = TRUE) - min(D, na.rm = TRUE))
    Y = 1/(1 + exp((DNormed * (-12)) + log(9999)))
    Y = as.numeric(Y)



   return(Y)
}
