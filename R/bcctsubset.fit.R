bcctsubset.fit <-
function(priornum, subset.index, maximal.mod, IP, eta.hat, ini.index, ini.beta, 
    ini.sig, iters, save, name, null.move.prob, a, b, progress) 
{
    if (is.null(name)) {
        name_RJACC <- "RJACC.txt"
        name_MHACC <- "MHACC.txt"
        name_BETA <- "BETA.txt"
        name_MODEL <- "MODEL.txt"
        name_SIG <- "SIG.txt"
    }    else {
        name_RJACC <- paste(name, "RJACC.txt", sep = "")
        name_MHACC <- paste(name, "MHACC.txt", sep = "")
        name_BETA <- paste(name, "BETA.txt", sep = "")
        name_MODEL <- paste(name, "MODEL.txt", sep = "")
        name_SIG <- paste(name, "SIG.txt", sep = "")
    }

	subset.mod<-c()
	for(i in 1:dim(subset.index)[1]){
	subset.mod[i]<-index2model(subset.index[i,])}
	big.X <- maximal.mod$x
    y <- maximal.mod$y
    data <- maximal.mod$data
    curr.index <- ini.index
    curr.mod<-index2model(curr.index)
    curr.X <- big.X[, curr.index == 1]
    curr.ivar <- IP[curr.index == 1, curr.index == 1]
    MODEL <- c()
    BETA <- c()
    curr.beta <- ini.beta
    SIG <- c()
    curr.sig <- ini.sig
    rj_acc <- c()
    mh_acc <- c()
    counter <- 0
    if (progress) {
        pb <- txtProgressBar(min = 0, max = iters, style = 3)
    }
    while (counter < iters) {
        
	uu<-runif(1)
	if(uu<null.move.prob){
	prop_a_model<-list(new.index=curr.index,type="null",total.choices=0,null.move.prob=null.move.prob)} else{
	prop_a_model<-list(new.index=model2index(sample(x=subset.mod[subset.mod!=curr.mod],size=1),dig=length(curr.index)),type="move",total.choices=length(subset.mod)-1,null.move.prob=null.move.prob)}

prop.index <- prop_a_model$new.index
        if (prop_a_model$type == "null") {
            new.beta <- iwls_mh(curr.y = y, curr.X = curr.X, 
                curr.beta = curr.beta, iprior.var = curr.ivar/curr.sig)
            if (all(new.beta == curr.beta)) {
                mh_acc <- c(mh_acc, 0)
            }
            else {
                mh_acc <- c(mh_acc, 1)
            }
            new.index <- curr.index
        }
        if (prop_a_model$type != "null") {
            rho_bot <- (1 - prop_a_model$null.move.prob)/prop_a_model$total.choices
            prop_a_rev <- list(total.choices = length(subset.mod)-1)
            rho_top <- (1 - prop_a_model$null.move.prob)/prop_a_rev$total.choices
            rj <- RJ_update_swap(prop.index = prop.index, curr.index = curr.index, 
                curr.beta = curr.beta, eta.hat = eta.hat, curr.y = y, 
                big.X = big.X, proposal.probs = c(rho_top, rho_bot), 
                i.prop.prior.var = IP[prop.index == 1, prop.index == 
                  1]/curr.sig, i.curr.prior.var = curr.ivar/curr.sig)
            new.beta <- rj$new.beta
            new.index <- rj$new.index
            if (all(curr.index == new.index)) {
                rj_acc <- c(rj_acc, 0)
            }
            else {
                rj_acc <- c(rj_acc, 1)
            }
        }
        if (priornum == 2) {
            iR <- IP[new.index == 1, new.index == 1]
            curr.sig <- 1/rgamma(n = 1, shape = 0.5 * (length(new.beta) - 
                1 + a), rate = 0.5 * (b + as.vector(matrix(new.beta[-1], 
                nrow = 1) %*% iR[-1, -1] %*% matrix(new.beta[-1], 
                ncol = 1))))
        }
        curry <- rep(0, dim(big.X)[2])
        curry[new.index == 1] <- new.beta
        curr.index <- new.index
		curr.mod<-index2model(curr.index)
        curr.beta <- new.beta
        curr.X <- big.X[, curr.index == 1]
        curr.ivar <- IP[curr.index == 1, curr.index == 1]
        BETA <- rbind(BETA, curry)
        SIG <- c(SIG, curr.sig)
        MODEL <- c(MODEL, curr.mod)
        counter <- counter + 1
        if (progress) {
            setTxtProgressBar(pb, counter)
        }
        if (save > 0) {
            if (counter%%save == 0) {
                write.table(file = name_BETA, x = BETA, row.names = FALSE, 
                  col.names = FALSE, append = TRUE)
                write.table(file = name_MODEL, x = MODEL, row.names = FALSE, 
                  col.names = FALSE, append = TRUE)
                write.table(file = name_SIG, x = SIG, row.names = FALSE, 
                  col.names = FALSE, append = TRUE)
                write.table(file = name_RJACC, x = rj_acc, row.names = FALSE, 
                  col.names = FALSE, append = TRUE)
                write.table(file = name_MHACC, x = mh_acc, row.names = FALSE, 
                  col.names = FALSE, append = TRUE)
                rj_acc <- c()
                mh_acc <- c()
                BETA <- c()
                MODEL <- c()
                SIG <- c()
            }
        }
    }
    if (progress) {
        close(pb)
    }
    list(BETA = BETA, SIG = SIG, MODEL = MODEL, rj_acc = rj_acc, 
        mh_acc = mh_acc)
}
