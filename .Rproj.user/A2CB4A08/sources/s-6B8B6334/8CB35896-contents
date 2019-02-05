#R function modified from Poisot et al. (2012) and Simanonok and  Burkle 2014 to calculate partitions of interaction turnover
######################################################################################################################################
####First function used to calculate partitions of interaction turnover for two sites (networks)

beta.sor = function(w1,w2){
    Z<-(w1)-(w2)
    #Presence/absence matrices only
    #w1 and w2 must have equal dimensions
    #Bcc is interaction turnover
    a<-sum(w1[w1==w2])
    b<-length(Z[Z==1])
    c<-length(Z[Z==-1])
    BWN<-((b+c)/(2*a+b+c))
    # For Z, if = 0 it must be either a shared interaction or a lack of an interaction

    #BP: turnover of interactions due to lower trophic level turnover (rows of matrices)
    #BPb: Z=1, rows must be empty in w2, columns must be filled in w2
    #BPc: Z=-1, rows must be empty in w1, columns must be filled in w1

    BPb<-vector("numeric",length=1)
    BPb1<-Z[rowSums(w2)==0,colSums(w2)>0]
    BPb1<-BPb1[BPb1==1]
    BPb<-length(BPb1)

    BPc<-vector("numeric",length=1)
    BPc1<-Z[rowSums(w1)==0,colSums(w1)>0]
    BPc1<-BPc1[BPc1==-1]
    BPc<-length(BPc1)

    #BA: turnover of interactions due to upper trophic level turnover (columns of matrices)
    #BAb: Z=1, rows must be filled in w2, columns must be empty in w2
    #BAc: Z=-1, rows must be filled in w1, columns must be empty in w1

    BAb<-vector("numeric",length=1)
    BAb1<-Z[rowSums(w2)>0,colSums(w2)==0]
    BAb1<-BAb1[BAb1==1]
    BAb<-length(BAb1)

    BAc<-vector("numeric",length=1)
    BAc1<-Z[rowSums(w1)>0,colSums(w1)==0]
    BAc1<-BAc1[BAc1==-1]
    BAc<-length(BAc1)

    #BPA: interaction turnover due to simultaneous turnover of both trophic levels
    #BPAb: Z=1, rows must be empty, columns must be empty
    #BPAc Z=-1, rows must be empty, columns must be empty

    BPAb<-vector("numeric",length=1)
    BPAb1<-Z[rowSums(w2)==0,colSums(w2)==0]
    BPAb1<-BPAb1[BPAb1==1]
    BPAb<-length(BPAb1)

    BPAc<-vector("numeric",length=1)
    BPAc1<-Z[rowSums(w1)==0,colSums(w1)==0]
    BPAc1<-BPAc1[BPAc1==-1]
    BPAc<-length(BPAc1)

    #See Novotny 2009;Trøjelsgaard et al 2015 for more detail
    BPA<-(BPAb+BPAc)/(2*a+b+c)
    BP<-(BPb+BPc)/(2*a+b+c)
    BA<-(BAb+BAc)/(2*a+b+c)

    #Bos is the rewiring, since BWN = BP + BA + BPA + BOS,
    #and we calculated all other partitions, we can subtract to find it
    BOS=(BWN-BP-BA-BPA)
    BOS<-ifelse(BOS<0.00001,0,BOS)
    #Whichever interactions do not change must be conserved,
    #therefore 1-BWN = proportion of links conserved between two networks
    prop.links.conserved=(1-BWN)

    #Additive total for species turnover: simultaneous + Animal+ plant
    BST=(BP+BA+BPA)

    return(list(prop.links.conserved = prop.links.conserved,
                BST = BST, BP = BP, BA = BA, BOS = BOS, BPA = BPA, BWN = BWN))
}
