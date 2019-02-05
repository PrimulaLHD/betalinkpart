#multi-site interaction dataframe to interation matices with equal dimensions
multi_equal <- function(df,varnames=c("Plant","Animal","WebID","freq")){
    network_list <- bipartite::frame2webs(df,varnames = varnames,type.out = "list",emptylist = FALSE)
    return(network_list)
}
