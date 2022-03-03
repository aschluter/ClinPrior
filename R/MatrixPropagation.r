#' Title
#'
#' @param alpha
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
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
