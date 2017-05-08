#' EcolUtils: Utilities for community ecology analysis.
#'
#' The package \pkg{EcolUtils} provides tools for community ecology analysis not present in up-to-date available packages. Designed with molecular 16S-derived data in mind.
#'
#'@details The package \pkg{EcolUtils} depends on \pkg{vegan} which can be installed from CRAN.
#' 
#' To see the preferable citation of the package, type \code{citation("EcolUtils")}.
#'@docType package
#'@name EcolUtils
#'@author Guillem Salazar <salazar@@icm.csic.es>

NULL


#' Rarefaction of a community matrix with permutations
#'
#' This function generates one randomly rarefied community data frame through \code{n} repeated independent rarefactions.
#' @param x Community data, a matrix-like object.
#' @param sample Subsample size (\code{min(rowSums(x))} as default)
#' @param n Number of independent rarefactions.
#' @param round.out logical; should output be rounded.
#' @details
#' Function \code{rrarefy.perm} generates one randomly rarefied community data frame by computing \code{n} rarefied communities with \code{rrarefy} in \pkg{vegan} and computing the mean.
#' The average value for each cell may (or may not) be rounded by using \code{round.out} parameter.
#' @keywords EcolUtils
#' @return Rarefied community.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' library(vegan)
#' data(varespec)
#' rrarefy.perm(varespec*100)

rrarefy.perm<-function(x,sample=min(rowSums(x)),n=100,round.out=T){
  require(vegan)
  y<-rrarefy(x,sample)
  for (i in 2:n){
    cat("Permutation ",i," out of ",n,"\n")
    y<-y+rrarefy(x,sample)	
  }
  if (round.out==T) y<-round(y/n)
  if (round.out==F) y<-y/n
  y}
  
#' Pairwise comparisons for Permutational Multivariate Analysis of Variance Using Distance Matrices
#'
#' Pairwise comparisons for all pairs of levels of a factor by using for Permutational MANOVA.
#' @param dist.mat Dissimilarity object
#' @param Factor Factor whose levels are to be compared.
#' @param nper Number of permutations.
#' @param corr.method P-value's correction method (from \code{p.adjust}).
#' @details Basically the \code{adonis.pair} function applies the \code{adonis} function from \pkg{vegan} to all pairs of levels of a factor. P-values are then corrected with \code{p.adjust}.
#' @keywords EcolUtils
#' @return Data frame with the R2, p-values and corrected p-values for each pairwise combination.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' library(vegan)
#' data(dune)
#' data(dune.env)
#' adonis.pair(vegdist(dune),dune.env$Management)

adonis.pair<-function(dist.mat,Factor,nper=1000,corr.method="fdr"){
  require(vegan)
  as.factor(Factor)
  comb.fact<-combn(levels(Factor),2)
  pv<-NULL
  R2<-NULL
  for (i in 1:dim(comb.fact)[2]){
    model.temp<-adonis(as.dist(as.matrix(dist.mat)[Factor==comb.fact[1,i] | Factor==comb.fact[2,i],Factor==comb.fact[1,i] | Factor==comb.fact[2,i]])~Factor[Factor==comb.fact[1,i] | Factor==comb.fact[2,i]],permutations=nper)
    pv<-c(pv,model.temp$aov.tab[[6]][1])
    R2<-c(R2,model.temp$aov.tab$R2[1])}
  pv.corr<-p.adjust(pv,method=corr.method)
  data.frame(combination=paste(comb.fact[1,],comb.fact[2,],sep=" <-> "),R2=R2,P.value=pv,P.value.corrected=pv.corr)}

#' Specialist/Generalist classification of OTUs based on niche width and permutation algorithms
#'
#' Classification of OTUs in generalists / specialists / non-significant based on the deviation of niche width indexes (\code{shanon}, \code{levins} or \code{occurrence}) from null values computed with permutation algorithms for community matrices.
#' @param comm.tab Community data, a matrix-like object (samples as rows; OTUs as columns).
#' @param niche.width.method Niche width index (from \code{niche.width} in \pkg{spaa}): \code{levins} (default) or \code{shannon}. Or simply the \code{occurrence}: the number of samples where an OTU occurs.
#' @param n Number of permutations.
#' @param  perm.method Method for null model construction (from \code{permatswap} in \pkg{vegan}). Currently, only \code{quasiswap} (default) has been thoroughly tested.
#' @param  probs Probabilities for confidence interval calculations.
#' @details Basically the \code{spec.gen} function computes a niche width index for each OTU in the \code{comm.tab}. The mean index value and CI for each OTU is computed for \code{n} null matrices created through permutation algorithms. Each OTU is classified as specialist / generalist / non significant if the real value is lower / higher / within the CI.
#' @keywords EcolUtils
#' @return Data frame with the observed niche width value, the mean and CI null values and the classification of each OTU.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' library(RCurl)
#' x<-getURL("https://raw.githubusercontent.com/GuillemSalazar/MolEcol_2015/master/OTUtable_Salazar_etal_2015_Molecol.txt")
#' comm.tab<-read.table(text=x,sep="\t",row.names=1,header=TRUE,comment.char="@@")
#' comm.tab<-t(comm.tab[,1:60])
#' comm.tab<-comm.tab[,which(colSums(comm.tab)>0)]
#' res<-spec.gen(comm.tab,n=100)
#' 
#' comm.tab.bin<-ceiling(comm.tab/max(comm.tab))
#' plot(colSums(comm.tab),colSums(comm.tab.bin)/dim(comm.tab.bin)[1],col=res$sign,pch=19,log="x",xlab="Abundance",ylab="Occurrence")
#' legend("bottomright",levels(res$sign),col=1:3,pch=19,inset=0.01,cex=0.7)

spec.gen<-function(comm.tab,niche.width.method="levins",perm.method="quasiswap",n=1000,probs=c(0.025,0.975)){
  require(spaa)
  require(vegan)
  occurrence<-function(x){apply(ceiling(x/max(x)),2,sum)}
  n<-n
  if (niche.width.method=="occurrence") levin.index.real<-occurrence(comm.tab) else levin.index.real<-as.numeric(niche.width(comm.tab,method=niche.width.method))
  names(levin.index.real)<-colnames(comm.tab)
  
  levin.index.simul<-matrix(NA,ncol=dim(comm.tab)[2],nrow=n)
  for (i in 1:n){
    if (niche.width.method=="occurrence") levin.index.simul[i,]<-occurrence(permatswap(comm.tab,perm.method,times=1)$perm[[1]]) else levin.index.simul[i,]<-as.numeric(niche.width(permatswap(comm.tab,perm.method,times=1)$perm,method=niche.width.method))
  }
  colnames(levin.index.simul)<-colnames(comm.tab)
  levin.index.simul<-as.data.frame(levin.index.simul)
  media<-apply(levin.index.simul,2,mean)
  ci<-apply(levin.index.simul,2,quantile,probs=probs)
  resultats<-data.frame(observed=levin.index.real,mean.simulated=media,lowCI=ci[1,],uppCI=ci[2,],sign=NA)
  for (j in 1:dim(resultats)[1]){
    if (resultats$observed[j]>resultats$uppCI[j]) resultats$sign[j]<-"GENERALIST"
    if (resultats$observed[j]<resultats$lowCI[j]) resultats$sign[j]<-"SPECIALIST"
    if (resultats$observed[j]>=resultats$lowCI[j] & resultats$observed[j]<=resultats$uppCI[j]) resultats$sign[j]<-"NON SIGNIFICANT"
  }
  resultats$sign<-as.factor(resultats$sign)
  resultats}

#' Seasonality test based on autocorrelation and null communities
#'
#' Classification of OTU's seasonality based on the sum of their auto-correaltion function and on null community matrices.
#' @param comm.tab Community data, a matrix-like object (samples as rows; OTUs as columns). Samples should be ordered and representing a time series.
#' @param n Number of permutations.
#' @param  probs Probabilities for confidence interval calculations.
#' @param lag.max Maximum lag at which to calculate the acf. See the \code{acf} function.
#' @details Basically the \code{seasonality.test} function computes the auto-correlation function (acf) for each OTU in the \code{comm.tab} through the \code{acf} function in the \pkg{stats} package. The sum of the absolute values of the acf is computed as the seasonality index for each OTU.
#' This seasonality index and CI for each OTU is also computed for \code{n} null community matrices. The null matrices are created by randomly shuffling the rows in \code{comm.tab}.  Each OTU is classified depending whether the real seasonality index is lower / higher / within the CI.
#' @keywords EcolUtils
#' @return Data frame with the observed seasonality index, the mean and CI null values and the classification of each OTU.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' library(RCurl)
#' # It runs but makes no ecological sense as data does not represent a time-series
#' x<-getURL("https://raw.githubusercontent.com/GuillemSalazar/MolEcol_2015/master/OTUtable_Salazar_etal_2015_Molecol.txt")
#' comm.tab<-read.table(text=x,sep="\t",row.names=1,header=TRUE,comment.char="@@")
#' comm.tab<-t(comm.tab[,1:60])
#' comm.tab<-comm.tab[,which(colSums(comm.tab)>0)]
#' res<-seasonality.test(comm.tab,n=10)

seasonality.test<-function(comm.tab,n=1000,probs=c(0.025, 0.975),lag.max=120,na.action=na.pass) 
{
  require(vegan)
  
  season.index<-function(x){
    acf.all<-apply(as.matrix(x),2,acf,plot=F,lag.max=lag.max,na.action=na.pass)
    acf.all<-sapply(acf.all,"[[",1)
    apply(acf.all,2,function(x) sum(abs(x)))
  }
  
  n<-n
  season.index.real<-season.index(comm.tab)
  
  names(season.index.real) <- colnames(comm.tab)
  season.index.simul<-matrix(NA, ncol = dim(comm.tab)[2],nrow = n)
  for (i in 1:n) {
    season.index.simul[i, ]<-season.index(comm.tab[sample(1:nrow(comm.tab)),])
  }
  colnames(season.index.simul) <- colnames(comm.tab)
  season.index.simul <- as.data.frame(season.index.simul)
  media <- apply(season.index.simul, 2, mean)
  ci <- apply(season.index.simul, 2, quantile, probs = probs)
  resultats <- data.frame(observed = season.index.real, mean.simulated = media,lowCI = ci[1, ], uppCI = ci[2, ], sign = NA)
  for (j in 1:dim(resultats)[1]) {
    if (resultats$observed[j] > resultats$uppCI[j]) 
      resultats$sign[j] <- "SIGNIFICANTLY HIGHER"
    if (resultats$observed[j] < resultats$lowCI[j]) 
      resultats$sign[j] <- "SIGNIFICANTLY LOWER"
    if (resultats$observed[j] >= resultats$lowCI[j] & resultats$observed[j] <= 
        resultats$uppCI[j]) 
      resultats$sign[j] <- "NON SIGNIFICANT"
  }
  resultats$sign <- as.factor(resultats$sign)
  resultats
}

#' Niche value computation for OTUs in a community through abundance-weighted mean and matrix randomization.
#'
#' Computation of the abundance-weighted mean of an environmental variable for all OTUs in a community and statistical comparison to randomized communities.
#' @param comm.tab Community data, a matrix-like object (samples as rows; OTUs as columns).
#' @param env.var Environmental variable as a numeric vector.
#' @param n Number of permutations.
#' @param  probs Probabilities for confidence interval calculations.
#' @details The \code{niche.val} function computes the abundance-weighted mean of an environmental variable for each OTU in the \code{comm.tab}. The mean value and CI for each OTU is computed for \code{n} null matrices. The null matrices are created by randomly shuffling the rows in \code{comm.tab}.  Each OTU is classified depending whether the real niche value is lower / higher / within the CI.
#' @keywords EcolUtils
#' @return Data frame with the observed niche value, the mean and CI null values and the classification of each OTU based on randomizations.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>

niche.val<-function(comm.tab,env.var,n=1000,probs=c(0.025,0.975)){
  require(vegan)
  stat.real<-apply(comm.tab,2,function (x) {weighted.mean(env.var,x,na.rm=T)})
  stat.simul<-matrix(NA,ncol=dim(comm.tab)[2],nrow=n)
  for (i in 1:n){
    print(paste("Rarefaction",i))
    stat.simul[i,]<-apply(comm.tab[sample(1:nrow(comm.tab)),],2,function (x) {weighted.mean(env.var,x,na.rm=T)})
  }
  colnames(stat.simul)<-colnames(comm.tab)
  simul<-as.data.frame(stat.simul)
  media<-apply(stat.simul,2,mean,na.rm=T)
  ci<-apply(stat.simul,2,quantile,probs=c(0.025,0.975),na.rm=T)
  resultats<-data.frame(observed=stat.real,mean.simulated=media,lowCI=ci[1,],uppCI=ci[2,],sign=NA)
  
  classify.sign<-function (x){
    if (is.na(x[1])) NA
    else if (x[1]>x[4]) "HIGHER"
    else if (x[1]<x[3]) "LOWER"
    else if (x[1]>=x[3] & x[1]<=x[4]) "NON SIGNIFICANT"}
  
  resultats$sign<-as.factor(apply(resultats,1,classify.sign))
  resultats
}
