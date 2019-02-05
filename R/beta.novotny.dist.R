###########################################################################################################################
###The furth function used for comparing more than two (Site) matrices


beta.novotny.dist = function(W,...){
    dWN = matrix(NA,ncol=length(W),nrow=length(W))
    colnames(dWN) = names(W)
    rownames(dWN) = names(W)
    dBcc = dWN
    dBPH = dWN
    dBP = dWN
    dBH = dWN
    dB0 = dWN
    for(i in c(1:(length(W)-1))){
        for(j in c((i+1):(length(W)))){
            partition = beta.novotny(W[[i]],W[[j]],...)
            dBcc[j,i]	= partition$Bcc
            dBPH[j,i]	= partition$BPH
            dBP[j,i]	= partition$BP
            dBH[j,i]	= partition$BH
            dB0[j,i]	= partition$B0
        }
    }
    distances = list(Bcc=dBcc, BPH=dBPH, BP=dBP, BH=dBH, B0=dB0)
    distances = lapply(distances,as.dist)
    return(distances)
}


