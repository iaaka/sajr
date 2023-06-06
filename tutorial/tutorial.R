#  do not forget to run read counter first. 
# 'java -jar sajr.jar count_reads'
setwd("your working directory") # the one where count_reads was run
#devtools::install_github('iaaka/sajr')
library(SAJR)

pdf('output.pdf')
data = loadSAData(ann.gff='example/a.gff',c(1,2))

#check how many segment we have
length(data)
#see first 10
data[1:10,]
#see annotation only
data$seg[1:10,]

# classify splicing events
data = setSplSiteTypes(data,'example/a.gff')
table(data$seg$sites) # aa denotes alternative acceptor sites; ad - cassete exons; ae - last exons; etc

data.f = data[data$seg$type %in% c('ALT','INT') & data$seg$position %in% c('LAST','INTERNAL','FIRST') & apply(data$i+data$e>=10,1,sum)==2 & apply(data$ir,1,sd,na.rm=TRUE) > 0,]
data.f[1:10,]
length(data.f)
# split genes into alternatives
all.alts = makeAlts(data$seg,'example/a.gff',remove.exn.ext = F)
par(mfrow=c(4,4))
plotAllAlts(all.alts,to.plot = 16)
no.int.ret.alt = filterAlts(all.alts,TRUE)
plotAllAlts(no.int.ret.alt,to.plot = 16)
par(mfrow=c(1,1))
#test for significant differences. Since there are no replicas, using of GLM is a bit redundant, but, just as example
mod = list(f=factor(c('treatment','control')))
data.f.glm = fitSAGLM(data.f,terms(x ~ f),mod)

# calculate p-value. Since there are no replicates, use binomial distribution
# newer use 'overdisp = FALSE' in real work.
data.f.pv = calcSAPvalue(data.f.glm,overdisp = FALSE)
data.f.pv[1:10,]
#plot histogramm of p-values
hist(data.f.pv[,2])
data.f.pv[1:10,]
#make BH correction
data.f.pv[,2] = p.adjust(data.f.pv[,2],method='BH')

#choose significant ones:
data.sign = data.f[data.f.pv[,2] <=0.05,]
length(data.sign)
# print top ten (by amplitude)
data.sign[order(abs(data.sign$ir[,1]-data.sign$ir[,2]),decreasing = TRUE),][1:10,]

##############################################################
###generate an a random data for more sofisticated analysis###
##############################################################
data$seg = data.f$seg[1:1000,]
data$i = matrix(rpois(10000,20),ncol=10)
t = rnorm(1000,20)
t1 = rnorm(1000,t,3)
t2 = rnorm(1000,t,3)

data$e = cbind(matrix(rpois(5000,rnorm(5000,t1)),ncol=5),matrix(rpois(5000,rnorm(5000,t2)),ncol=5))
data$ir = data$i/(data$i+data$e)

colnames(data$i)=colnames(data$e)=colnames(data$ir)=c(paste('c',1:5,sep=''),paste('t',1:5,sep=''))
rownames(data$seg) = rownames(data$i)=rownames(data$e)=rownames(data$ir)=1:1000

# check selfonsistency
plotCorrHM(data$ir,ColSideColors=rep(c('red','blue'),each=5))
plotMDS(data$ir,pch=rep(c(7,19),each=4),col=rep(c('red','blue'),each=5),main='All segments')
plotMDS(data$ir[data$seg$sites == 'ad',],pch=rep(c(7,19),each=4),col=rep(c('red','blue'),each=5),main='Cassette exons')

mod = list(f=factor(c(rep('treatment',times=5),rep('control',times=5))))
data.glm = fitSAGLM(data,terms(x ~ f),mod)
data.pv = calcSAPvalue(data.glm)
data.pv.corr = p.adjust(data.pv[,2],method='BH')
data.sign = data[data.pv.corr <= 0.1,]
clustSegs(data.sign,norm.groups=rep(1,times=10),k=2,rows=1,cols=2,col=c(rep('green',times=5),rep('red',times=5)))
dev.off()
