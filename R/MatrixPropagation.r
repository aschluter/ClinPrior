#' Propagation of the phenotypic score within of a multilayer network with physical and functional interactions
#'
#' @param alpha numeric value [0-1] that ponderates the propagation (alpha) and the previous knowledge (Y) (1-alpha) contribution in the final phenotypic score in the iterative propagation process.
#' @param Y vector with the phenotypic score obtained from the proteinScore function.
#'
#' @return matrix with final phenotypic scores after iterative propagation in a physical protein-protein interaction network and in a functional one.
#' @export
#'
#' @examples
#' HPOpatient = c("HP:0004481","HP:0002376","HP:0001257","HP:0001250","HP:0000238","HP:0002922","HP:0000365")
#' Y<-proteinScore(HPOpatient)
#' ClinPriorGeneScore<-MatrixPropagation(Y,alpha=0.2)
MatrixPropagation<-function(Y,alpha=0.2)
{
    library(Matrix)

    #Physical network

    pos <- match(total_unique_Conn_Physical, total_unique[,2])
    #alpha = c(0.01,0.1,0.2,0.3,0.4,0.5,0.6)

    F <- Y
    F2 <- F[pos]
    memo<-F2
    l= 1
    memoL = matrix(c(1:20),ncol = 20)
    bol = 1


    while(bol)
    {
        WF = normPhysical %*% F2
        F2 = (alpha * WF) + (1 - alpha) * F[pos]
        l = length(c(1:length(total_unique_Conn_Physical))[total_unique_Conn_Physical[order(as.matrix(memo))]==total_unique_Conn_Physical[order(as.matrix(F2))]])
        memoL = c(memoL,l)
        bol = mean(memoL[c((length(memoL)-19):length(memoL))])!=l
        memo = F2
    }
    F[pos] <- F2

    resultatPhysical<-cbind(total_unique[order(F,decreasing=T),],F[order(F,decreasing=T)])


    #Functional network

    pos <- match(total_unique_Conn_Func, total_unique[,2])
    #alpha = c(0.01,0.1,0.2,0.3,0.4,0.5,0.6)

    F <- Y
    F2 <- F[pos]
    memo<-F2
    l= 1
    memoL = matrix(c(1:20),ncol = 20)
    bol = 1

    while(bol)
    {
        WF = normFunc %*% F2
        F2 = (alpha * WF) + (1 - alpha) * F[pos]
        l = length(c(1:length(total_unique_Conn_Func))[total_unique_Conn_Func[order(as.matrix(memo))]==total_unique_Conn_Func[order(as.matrix(F2))]])
        memoL = c(memoL,l)
        bol = mean(memoL[c((length(memoL)-19):length(memoL))])!=l
        memo = F2
    }
    F[pos] <- F2

    resultatFunctional<-cbind(total_unique[order(F,decreasing=T),],F[order(F,decreasing=T)])
    ClinPriorScore <- cbind(resultatFunctional,resultatPhysical)
    colnames(ClinPriorScore) <-c("Symbol", "geneID", "PriorFunct", "Symbol", "geneID", "PriorPhys")

    gc()
    return(ClinPriorScore)
}
