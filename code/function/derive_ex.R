get.ex <- function(ages, mx){
        
        m <- length(ages)
        n <- c(diff(ages), NA)
        ax <- n/2
        ax[1] <- 0.07 + 1.7*mx[1]
        ax[2] <- 1.5
        ax[m] <- 1 / mx[m]
        qx  <- n*mx / (1 + (n - ax) * mx)
        qx[m] <- 1
        px  <- 1-qx
        lx  <- cumprod(c(1,px))*100000
        dx  <- -diff(lx)
        Lx  <- n*lx[-1] + ax*dx
        lx <- lx[-(m+1)]
        Lx[m] <- lx[m]/mx[m]
        Lx[is.na(Lx)] <- 0 ## in case of NA values
        Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
        Tx  <- rev(cumsum(rev(Lx)))
        ex  <- Tx/lx
        return(ex)
}