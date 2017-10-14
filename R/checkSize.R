checkSize <- function(exp.test, treshold=0.05) {
    t0 <- c()
    t1 <- c()

    for(i in 1:10000) {
        smpl1 <- rexp(20, 1)
        stat <- exp.test(smpl1)

        stat_val <- 0
        if( typeof(stat) == "double" ) {
            stat_val <- stat
        }
        else {
            stat_val <- stat$statistic
        }
        t0 <- c(t0, stat_val)

        smpl10 <- rexp(20, 10)
        stat_val <- exp.test(smpl10)
        stat_val <- 0
        if( typeof(stat) == "double" ) {
            stat_val <- stat
        }
        else {
            stat_val <- stat$statistic
        }
        t1 <- c(t1, stat_val)
    }

    return(1-ecdf(t1)(quantile(t0,1-treshold)))
}
