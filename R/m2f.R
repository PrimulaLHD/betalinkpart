# interaction matrix to data frame
m2f <- function(m){
    df <- as.data.frame.table(m)
    names(df) <- c("Plant","Animal","Value")
    df<-subset(df,Value>0)
    df
}
