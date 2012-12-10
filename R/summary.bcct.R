summary.bcct <-
function(object,n.burnin=0,thin=1,cutoff=0.75,statistic="X2",best=NULL,scale=0.1,prob.level=0.95,...){

is1<-inter_stats(object,n.burnin=n.burnin,cutoff=cutoff,thin=thin,prob.level=prob.level)
is2<-mod_probs(object,n.burnin=n.burnin,scale=scale,best=best,thin=thin)
is3<-bayespval(object,n.burnin=n.burnin,thin=thin,statistic=statistic)

est<-list(BETA=object$BETA,MODEL=object$MODEL,SIG=object$SIG,rj_acc=object$rj_acc,mh_acc=object$mh_acc,priornum=object$priornum,maximal.mod=object$maximal.mod,IP=object$IP,eta.hat=object$eta.hat,save=object$save,name=object$name,int_stats=is1,mod_stats=is2,pval_stats=is3)

class(est)<-"sbcct"

est

}
