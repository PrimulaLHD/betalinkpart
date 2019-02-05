###########################################################################################################################
###The second function used for comparing more than two Sites (networks)


beta.sor.dist = function(W,...){
    dWN = matrix(NA,ncol=length(W),nrow=length(W))
    colnames(dWN) = names(W)
    rownames(dWN) = names(W)
    dBWN = dWN
    dBPA = dWN
    dBP = dWN
    dBA = dWN
    dBOS = dWN
    for(i in c(1:(length(W)-1))){
        for(j in c((i+1):(length(W)))){
            partition = beta.sor(W[[i]],W[[j]],...)
            dBWN[j,i]	= partition$BWN
            dBPA[j,i]	= partition$BPA
            dBP[j,i]	= partition$BP
            dBA[j,i]	= partition$BA
            dBOS[j,i]	= partition$BOS
        }
    }
    distances = list(BWN=dBWN, BPA=dBPA, BP=dBP, BA=dBA, BOS=dBOS)
    distances = lapply(distances,as.dist)
    return(distances)
}

