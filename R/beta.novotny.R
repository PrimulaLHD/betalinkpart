#R function modified from Poisot et al. (2012) to calculate partitions of interaction turnover according to Novotny (2009).
# Simanonok and  Burkle 2014
######################################################################################################################################
####Third function used to calculate partitions of interaction turnover according to Novotny (2009), for two site(matrices)

beta.novotny = function(w1,w2){
    Z<-(w1)-(w2)
    #Presence/absence matrices only
    #w1 and w2 must have equal dimensions
    #Bcc is interaction turnover
    a<-sum(w1[w1==w2])
    b<-length(Z[Z==1])
    c<-length(Z[Z==-1])
    Bcc<-((b+c)/(a+b+c))
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

    #BH: turnover of interactions due to upper trophic level turnover (columns of matrices)
    #BHb: Z=1, rows must be filled in w2, columns must be empty in w2
    #BHc: Z=-1, rows must be filled in w1, columns must be empty in w1

    BHb<-vector("numeric",length=1)
    BHb1<-Z[rowSums(w2)>0,colSums(w2)==0]
    BHb1<-BHb1[BHb1==1]
    BHb<-length(BHb1)

    BHc<-vector("numeric",length=1)
    BHc1<-Z[rowSums(w1)>0,colSums(w1)==0]
    BHc1<-BHc1[BHc1==-1]
    BHc<-length(BHc1)

    #BPH: interaction turnover due to simultaneous turnover of both trophic levels
    #BPHb: Z=1, rows must be empty, columns must be empty
    #BPHc Z=-1, rows must be empty, columns must be empty

    BPHb<-vector("numeric",length=1)
    BPHb1<-Z[rowSums(w2)==0,colSums(w2)==0]
    BPHb1<-BPHb1[BPHb1==1]
    BPHb<-length(BPHb1)

    BPHc<-vector("numeric",length=1)
    BPHc1<-Z[rowSums(w1)==0,colSums(w1)==0]
    BPHc1<-BPHc1[BPHc1==-1]
    BPHc<-length(BPHc1)

    #See Novotny 2009 for more detail
    BPH<-(BPHb+BPHc)/(a+b+c)
    BP<-(BPb+BPc)/(a+b+c)
    BH<-(BHb+BHc)/(a+b+c)

    #B0 is the rewiring, since Bcc = BP + BH + BPH + B0,
    #and we calculated all other partitions, we can subtract to find it
    B0=(Bcc-BP-BH-BPH)
    B0<-ifelse(B0<0.00001,0,B0)
    #Whichever interactions do not change must be conserved,
    #therefore 1-Bcc = proportion of links conserved between two networks
    prop.links.conserved=(1-Bcc)

    #Additive total for species turnover: simultaneous + herbivore + plant
    BST=(BP+BH+BPH)

    return(list(prop.links.conserved = prop.links.conserved, BST = BST, BP = BP, BH = BH, B0 = B0, BPH = BPH, Bcc = Bcc))
}

