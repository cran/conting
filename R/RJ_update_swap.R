RJ_update_swap <-
function (prop.index, curr.index, curr.beta, eta.hat, curr.y, 
    big.X, proposal.probs, i.prop.prior.var, i.curr.prior.var){

icurrR <- i.curr.prior.var[-1, -1]
ipropR <- i.prop.prior.var[-1, -1]
currR <- chol2inv(chol(icurrR))
propR <- chol2inv(chol(ipropR))
RHO_TOP <- proposal.probs[1]
RHO_BOT <- proposal.probs[2]
curr.X <- big.X[, curr.index == 1]
prop.X <- big.X[, prop.index == 1]
Xcc <- big.X[, curr.index==1 & prop.index==1]
S1 <- matrix(big.X[,curr.index==1 & prop.index==0],nrow = length(curr.y))
S2 <- matrix(big.X[,curr.index==0 & prop.index==1],nrow = length(curr.y))
w <- exp(eta.hat)
curr.LP <- as.vector(curr.X %*% matrix(curr.beta, ncol = 1))

cXW <- t(Xcc) %*% diag(w)
cXWX <- cXW %*% Xcc
icXWX <- chol2inv(chol(cXWX))
icXWX.XW <- icXWX %*% cXW
SWIP1 <- t(S1) %*% diag(w) - t(S1) %*% t(cXW) %*% icXWX.XW
SWIP2 <- t(S2) %*% diag(w) - t(S2) %*% t(cXW) %*% icXWX.XW

## START OF SWAP MOVE

if(dim(S1)[2]>0 & dim(S2)[2]>0){ 

SIG1 <- chol2inv(chol(SWIP1%*%S1))
SIG2 <- chol2inv(chol(SWIP2%*%S2))
MU1<-SIG1%*%SWIP1%*%matrix(eta.hat,ncol=1)
MU2<-SIG2%*%SWIP2%*%matrix(eta.hat,ncol=1)

beta_1 <- curr.beta[prop.index[curr.index == 1] == 1]
beta_2 <- curr.beta[prop.index[curr.index == 1] == 0]
u<-as.vector(rmvnorm(n=1,mean=MU2,sigma=SIG2))
prop.beta <- rep(0, dim(prop.X)[2])

prop.beta[curr.index[prop.index == 1] == 1] <- beta_1 + as.vector(icXWX.XW%*%(S1%*%matrix(beta_2,ncol=1)-S2%*%matrix(u,ncol=1)))
prop.beta[curr.index[prop.index == 1] == 0] <- u

prop.LP <- as.vector(prop.X %*% matrix(prop.beta, ncol = 1))

top <- sum(curr.y * prop.LP) - sum(exp(prop.LP)) + dmvnorm(x = prop.beta[-1],mean = rep(0, length(prop.beta) - 1), sigma = propR, log = TRUE)
bot <- sum(curr.y * curr.LP) - sum(exp(curr.LP)) + dmvnorm(x = curr.beta[-1],mean = rep(0, length(curr.beta) - 1), sigma = currR, log = TRUE)

jac<-dmvnorm(x=beta_2,mean=MU1,sigma=SIG1,log=TRUE)-dmvnorm(x=u,mean=MU2,sigma=SIG2,log=TRUE)

prob <- (RHO_TOP/RHO_BOT) * exp(top - bot + jac)}

## END OF SWAP MOVE

## START OF DEATH MOVE

if(dim(S1)[2]>0 & dim(S2)[2]==0){ 

SIG1 <- chol2inv(chol(SWIP1%*%S1))
MU1<-SIG1%*%SWIP1%*%matrix(eta.hat,ncol=1)

beta_1 <- curr.beta[prop.index[curr.index == 1] == 1]
beta_2 <- curr.beta[prop.index[curr.index == 1] == 0]
prop.beta <- rep(0, dim(prop.X)[2])

prop.beta[curr.index[prop.index == 1] == 1] <- beta_1 + as.vector(icXWX.XW%*%(S1%*%matrix(beta_2,ncol=1)))

prop.LP <- as.vector(prop.X %*% matrix(prop.beta, ncol = 1))

top <- sum(curr.y * prop.LP) - sum(exp(prop.LP)) + dmvnorm(x = prop.beta[-1],mean = rep(0, length(prop.beta) - 1), sigma = propR, log = TRUE)
bot <- sum(curr.y * curr.LP) - sum(exp(curr.LP)) + dmvnorm(x = curr.beta[-1],mean = rep(0, length(curr.beta) - 1), sigma = currR, log = TRUE)

jac<-dmvnorm(x=beta_2,mean=MU1,sigma=SIG1,log=TRUE)

prob <- (RHO_TOP/RHO_BOT) * exp(top - bot + jac)}

## END OF DEATH MOVE

## START OF BIRTH MOVE

if(dim(S1)[2]==0 & dim(S2)[2]>0){ 

SIG2 <- chol2inv(chol(SWIP2%*%S2))
MU2<-SIG2%*%SWIP2%*%matrix(eta.hat,ncol=1)

u<-as.vector(rmvnorm(n=1,mean=MU2,sigma=SIG2))

prop.beta <- rep(0, dim(prop.X)[2])
prop.beta[curr.index[prop.index == 1] == 1] <- curr.beta + as.vector(icXWX.XW%*%(-S2%*%matrix(u,ncol=1)))
prop.beta[curr.index[prop.index == 1] == 0] <- u

prop.LP <- as.vector(prop.X %*% matrix(prop.beta, ncol = 1))

top <- sum(curr.y * prop.LP) - sum(exp(prop.LP)) + dmvnorm(x = prop.beta[-1],mean = rep(0, length(prop.beta) - 1), sigma = propR, log = TRUE)
bot <- sum(curr.y * curr.LP) - sum(exp(curr.LP)) + dmvnorm(x = curr.beta[-1],mean = rep(0, length(curr.beta) - 1), sigma = currR, log = TRUE)

jac<--dmvnorm(x=u,mean=MU2,sigma=SIG2,log=TRUE)

prob <- (RHO_TOP/RHO_BOT) * exp(top - bot + jac)}

if (prob >= runif(1)) {
        new.beta <- prop.beta
        new.index <- prop.index
    }
    else {
        new.beta <- curr.beta
        new.index <- curr.index
    }
list(new.beta = new.beta, new.index = new.index)}
