#' redefine operator for class sajr
#' @export
#' @import MASS
#' @import data.table
#' @import GenomicRanges
"[.sajr" = function(a,i,j){
	r = list()
	for(p in 1:length(names(a))){
		if(p == 1)
			r[[names(a)[p]]] = a[[p]][i,,drop=FALSE]
		else
			r[[names(a)[p]]] = a[[p]][i,j,drop=FALSE]
	}
	class(r)=class(a)
	r
}

#' redefine length fo class sajr
#' @export
length.sajr = function(a){
	length = nrow(a[[1]])
}

getAttr = function(attrs,name){
	attrs = attrs[1]
	attrs = substr(attrs,regexpr(name,attrs,fixed=T)[[1]]+nchar(name),nchar(attrs))
	i = regexpr(";",attrs,fixed=T)[[1]]-1
	if(i == -2)
		i = nchar(attrs)
	substr(attrs,1,i)
}

#' Plots fold change vs expression
#'
#' Calculates mean for each row and then plots for each column (mean - column) vs mean.
#' Adds density and lowess curve. In case of expression, log transformed data is preferable
#'
#' @param e numerical matrix. Expression or inclusion rato.
#' @param xlab graphical parameter
#' @param ylab graphical parameter
#'  
#' @export
plotFC2Exp = function(e,xlab='Mean expression',ylab='Fold change'){
	s = apply(e,1,mean,na.rm=T)
	for(i in 1:dim(e)[2]){
		t = s-e[,i]
		s_ = s[!is.na(t)]
		t = t[!is.na(t)]
		d = kde2d(s_,t,n=100)
		image(d,xlab=xlab,ylab=ylab,main=colnames(e)[i],col=rev(heat.colors(100)))
		contour(d,add=T)
		points(s_,t,type='p',col='#77777706',pch='.')		
		lines(lowess(s_,t),lwd=2,col='blue')
		abline(h=mean(t),col='green')
	}
}

#' Loads gff annotation and expression data from files produces by java SplicingAnalyser
#'
#' This function loads annotaton, read counts for all genes
#'
#' @param ann.gff name of annotation file (relative to base)
#' @param fnames vector of files (SplicingAnalyser output, relative to base, without extensions)
#' @param lib_names vector of names of libraries, fnames used instead if NA
#' @param base base name for all files (optional)
#' 
#' @return object of class sajr and list with 2 dataframes: gene annotation (gene), 
#' and read counts (cnts)
#' finction length and '[' operator redefined for clas sajr: first gives number of segments,
#' second substracts segments and samples
#'  
#' @export
loadGData = function(ann.gff,fnames=c(),lib_names=NA,base=''){
	if(is.na(lib_names[1]))
		lib_names = fnames
	ann.gff = paste(base,ann.gff,sep='')
	fnames = paste(base,fnames,sep='')
	gff = read.table(ann.gff,comment.char='#',sep="\t",as.is=T)
	print("annotation loaded")
	gff = gff[gff[,3] == 'gene',]
	gene = as.data.frame(matrix(ncol=4,nrow = dim(gff)[1]))
	colnames(gene) = c('chr_id','start','stop','strand')
	rownames(gene) = sapply(gff[,9],getAttr,name='gene_id=')
	gene[,1:3] = gff[,c(1,4,5)]
	gene[gff[,7] == '+',4] = 1
	gene[gff[,7] == '-',4] =-1
	cnts = matrix(ncol=length(lib_names),nrow = dim(gene)[1])
	rownames(cnts) = rownames(gene)
	print("annotation parsed")
	if(length(lib_names)>0)
		for(j in 1:length(fnames)){
			cat('\rload',j,fnames[j])
			t = read.table(paste(fnames[j],'.gene',sep=''),header=T,as.is=T,sep='\t',comment.char="#")
			rownames(t) = t[,1]
			cnts[,j] = t[rownames(gene),2]
			gc(verbose = FALSE)
		}
	colnames(cnts) = lib_names
	r = list(gene=gene,cnts=cnts)
	class(r)=c('sajr',class(r))
	r
}


#' Calculates length of constitutive part of gene
#'
#' It could be used to normalize gene read counts to calculate [FR]PKM
#'
#' @param segs output of loadSAData function
#' @param with_first_last should constitutive last of first segment (if any) be used. This 
#' option should be used in agreement with count_only_internal option of SAJR read counter
#' 
#' @return vector of lengths
#'  
#' @export
calcGeneConstLength = function(segs,with_first_last){
	genes = unique(segs$seg[,1])
	filter = segs$seg$type=='EXN'
	if(!with_first_last)
		filter = filter & segs$seg$position %in% c('INTERNAL','ONLY')
	segs = segs$seg[filter,]
	lens = segs[,4]-segs[,3]+1
	tapply(lens,factor(segs[,1]),sum)
}

#' Adds splicing site info to segments annotation
#'
#' Adds columns 'sites' to d$seg. Column contains information about splce sites
#' at segment boundaries: 'a' is acceptor, 'd' is donor, '.' means that both sites coincide.
#' 's' and 'e' denotes TSS or polyA sites respectively
#'
#' @param d sajr object (from loadSAData function)
#' @param gff.fname name of annotation file
#' 
#' @return sajr object with additional info
#'  
#' @export
setSplSiteTypes = function(d,gff.fname){
	ints = loadIntrons(gff.fname)
	sts1 = paste(ints$gene_id,ints$start-1,sep='.')
	sts2 = paste(ints$gene_id,ints$stop,sep='.')
	don = unique(c(sts1[ints[,5]== 1],sts2[ints[,5]==-1]))
	acc = unique(c(sts1[ints[,5]==-1],sts2[ints[,5]== 1]))
	
	seg.site1 = paste(d$seg$gene_id,d$seg$start-1,sep='.')
	seg.site2 = paste(d$seg$gene_id,d$seg$stop   ,sep='.')
	
	f =function(s){
		s = cbind(s %in% don,s %in% acc)
		s = ifelse(s[,1] & s[,2],'.',ifelse(!s[,1] & !s[,2],'-',ifelse(s[,1],'d','a')))
		ifelse(s == '-',ifelse(d$seg$pos=='FIRST','s','e'),s)
	}
	
	seg.site1 = f(seg.site1)
	seg.site2 = f(seg.site2)
	sites = paste(ifelse(d$seg$strand == 1,seg.site1,seg.site2),
								ifelse(d$seg$strand == 1,seg.site2,seg.site1),sep='')
	sites[is.na(d$seg$strand)] = '--'
	sites[!is.na(d$seg$strand) & d$seg$position=='ONLY'] = 'se'
	d$seg$sites = sites
	d
}

#' Calculates RPKM for gene
#'
#' Normolizes read counts to length of constitutive part of gene and sum of gene read counts
#'
#' @param segs output of loadSAData function
#' @param genes output of loadGData function
#' @param with_first_last should constitutive last of first segment (if any) be used. This 
#' option should be used in agreement with count_only_internal option of SAJR read counter
#' 
#' @return genes with $rpkm table
#'  
#' @export
calcRPKM = function(segs,genes,with_first_last){
	len = calcGeneConstLength(segs,with_first_last)
	len = len[rownames(genes$gene)]
	genes$rpkm = sweep(genes$cnts,2,apply(genes$cnts,2,sum)/1e6,'/')
	genes$rpkm = sweep(genes$rpkm,1,len/1e3,'/')
	genes
}

#' Loads intron annotations from gff file
#'
#'
#' @param ann.gff name of annotation file
#' 
#' @return intron_id,gene_id,chromosome,strand,start and stop
#'  
#' @export
loadIntrons= function(ann.gff){
	gff = read.table(ann.gff,comment.char='#',sep="\t",as.is=T)
	print("annotation loaded")
	gff = gff[gff[,3] == 'intron',]
	int = as.data.frame(matrix(ncol=5,nrow = dim(gff)[1]))
	colnames(int) = c("gene_id",'chr_id','start','stop','strand')
	int[,1] = sapply(gff[,9],getAttr,name='gene_id=')
	int[,2:4] = gff[,c(1,4,5)]
	int[gff[,7] == '+',5] = 1
	int[gff[,7] == '-',5] =-1
	int	
}


#' Loads gff annotation and splicing data from files produces by java SplicingAnalyser
#'
#' This function loads annotaton, read counts and inclusion ratio for all segments
#'
#' @param ann.gff name of annotation file (relative to base)
#' @param fnames vector of files (SplicingAnalyser output, relative to base, without extensions)
#' @param lib_names vector of names of libraries, fnames used instead if NA
#' @param base base name for all files (optional)
#' 
#' @return object of class sajr and  list with 4 dataframes: segment annotation (seg), 
#' inclusion ratio (ir), inclusion read count (i) and exclusion read count (e)
#' finction length and '[' operator redefined for clas sajr: first gives number of segments,
#' second substracts segments and samples
#'  
#' @export
loadSAData = function(ann.gff,fnames=c(),lib_names=NA,base=''){
	if(is.na(lib_names[1]))
		lib_names = fnames
	ann.gff = paste(base,ann.gff,sep='')
	fnames = paste(base,fnames,sep='')
	gff = read.table(ann.gff,comment.char='#',sep="\t",as.is=T)
	print("annotation loaded")
	gff = gff[gff[,3] == 'segment',]
	seg = as.data.frame(matrix(ncol=7,nrow = dim(gff)[1]))
	colnames(seg) = c("gene_id",'chr_id','start','stop','strand','type','position')
	rownames(seg) = sapply(gff[,9],getAttr,name='segment_id=')
	seg[,1] = sapply(gff[,9],getAttr,name='gene_id=')
	seg[,2:4] = gff[,c(1,4,5)]
	seg[gff[,7] == '+',5] = 1
	seg[gff[,7] == '-',5] =-1
	seg[,6] = sapply(gff[,9],getAttr,name='type=')
	seg[,7] = sapply(gff[,9],getAttr,name='position=')
	ir = matrix(ncol=length(lib_names),nrow = dim(seg)[1])
	i = matrix(ncol=length(lib_names),nrow = dim(seg)[1])
	e = matrix(ncol=length(lib_names),nrow = dim(seg)[1])
	rownames(ir) = rownames(i) = rownames(e) = rownames(seg)
	print("annotation parsed")
	if(length(lib_names)>0)
		for(j in 1:length(fnames)){
			cat('\rload',j,fnames[j])
			t = read.table(paste(fnames[j],'.seg',sep=''),header=T,as.is=T,sep='\t',comment.char="#")
			rownames(t) = t[,1]
			i[,j] = t[rownames(i),2]
			e[,j] = t[rownames(i),3]
			ir[,j] = as.numeric(t[rownames(i),4])
			ir[is.nan(ir[,j]),j] = NA
			gc(verbose = FALSE)
		}
	colnames(ir) = colnames(i) = colnames(e) = lib_names
	r = list(seg=seg,ir=ir,i=i,e=e)
	class(r)=c('sajr',class(r))
	r
}


#' Plots 2-dimensial MDS 
#' 
#' Function plots 2-dimensial MDS for all supplied features
#'
#' @param d numerical matrix (row*col=features*samples) of feature intensities. rows are segments, cols are samples
#' @param dist.method type of distance to be calculated. If first element is 'cor' than 1-cor distance is used with method specified as second element of dist.method parameter.
#' if first element is not 'cor' it is considered as method for dist function. For example dist.method='euclidean' or dist.method = c('cor','sp')
#' @param points points to plot. Will used instead of d and dist if specified.
#' @param dist distance matrix. Will used instead of d if specified
#' @param pch plot points instead of sample names
#' @param col graphical parameter
#' @param main graphical parameter
#' @param text sample names (colnames(d) used as default)
#' 
#' @export
plotMDS = function(d=NULL,dist.method=c('cor','sp'),points=NULL,dist=NULL,pch=NULL,col='black',main='',text=NULL,...){
	if(!is.null(points))
		mds = list(points=points)
	else{
		if(!is.null(dist))
			dst = as.dist(dist)
		else if(dist.method[1] == 'cor')
			dst = as.dist(1-cor(d,use='p',method = dist.method[2]))
		else
			dst = dist(t(d),method = dist.method)
		mds = isoMDS(dst,k=2)
	}
	
	if(is.null(text))
		text = rownames(mds$points)
	plot(1,type='n',xlim=range(mds$points[,1])*1.1,ylim=range(mds$points[,2])*1.1,xlab="Dimension 1",
			 ylab='Dimension 2',main=main)
	if(is.null(pch))
		text(mds$points[,1],mds$points[,2],labels=text,col=col,...)
	else
		points(mds$points[,1],mds$points[,2],col=col,pch=pch,...)
	invisible(mds)
}


#' Fit binomial GLM
#' 
#' Fit GLM model with binomial distribution for each supplied segments
#'
#' @param data object of sajr class, produced by \code{\link{loadSAData}}
#' @param formula object returned by \code{\link{terms}} function. Should be in form of x ~ sum_of_factors
#' @param terms list with values of factors used in formula
#' @param pseudocount proportion of total (inclusion + exclusion) reads that will be added to both inclusion and exclusion counts to evide overfitting
#' @param .parallel, .progress parameters for plyr::alply
#' @param return.pv logical, should function return p-values instead of models. Makes calculations more memory-efficient.
#' @param overdisp should overdispersion be taked into account. Has effect only if return.pv is true, see \code{\link{calcSAPvalue}} for details
#' @param disp.param  to use external disperison parameters. Has effect only if return.pv is true.
#' 
#' @return list of glm objects
#' 
#' @examples
#' mod = list(brain_reg=factor(c('pfc','pfc','cbc','cbc')),age=factor(c('nb','adult','nb','adult')))
#' alt.glm = fitSAGLM(seg.filtered,terms(x ~ brain_reg + age,keep.order=T),mod,0.05)
#' 
#' @import plyr
#' @export
fitSAGLM = function(data,formula,terms,pseudocount=0,.parallel=FALSE,.progress='none',return.pv=FALSE,overdisp=TRUE,disp.param=NULL){
	formula = terms(formula)
	term.names = attr(formula,'term.labels')
	r = alply(1:length(data),1,function(i){
		t = cbind(data$i[i,],data$e[i,])
		t = t + (t[,1]+t[,2])*pseudocount
		terms$x = t
		t = tryCatch(glm(formula,data=terms,family='quasibinomial'),
											error=function(e){warning(paste(e$message));return(NA)},
											warning=function(w){warning(paste('Segment ',rownames(data$seg)[i],': ',w$message,sep=''));return(NA)})
		if(return.pv)
			t = getPvForGLM(t,overdisp = overdisp,disp.param = disp.param[i],term.names = term.names,sid = rownames(data$seg)[i])
		t
		}
		,.parallel=.parallel,.progress=.progress)
	if(return.pv){
		r = do.call(rbind,r)
		rownames(r) = rownames(data$seg)
		colnames(r) = c('overdispersion',term.names)
	}else{
		names(r) = rownames(data$seg)
		attr(r,'term.labels')=term.names
	}
	r
}

#' Fit binomial GLM for gene expression data
#' Reads from each sample considered as mapped within particular gene (success in binomial trials) and other (failure)
#' 
#' @param data object of sajr class, produced by \code{\link{loadGData}}
#' @param formula object returned by \code{\link{terms}} function. Should be in form of x ~ sum_of_factors
#' @param terms list with values of factors used in formula
#' @param pseudocount proportion of total (inclusion + exclusion) reads that will be added to both inclusion and exclusion counts to evide overfitting
#' @param .parallel, .progress parameters for plyr::alply 
#' 
#' @return list of glm objects
fitDEGLM = function(ge,formula,terms,pseudocount,.parallel=FALSE,.progress='none'){
	t = apply(ge$cnts,2,sum)
	d = list(i = ge$cnts,e=sweep(-ge$cnts,2,t,'+'))
	class(d) = c('sajr','list')
	fitSAGLM(d,formula,terms,pseudocount,.parallel=.parallel,.progress=.progress)
}

#' Calculates p-values for GLM fits
#' 
#' Calculate p-value based on chisq distribution for each fit and factor
#' Could account for overdispersion (see parameter qbinom)
#'
#' @param data.glm list of glm fits from \code{\link{fitSAGLM}}
#' @param overdisp define how to deal with overdispersion:
#' 	if overdisp is TRUE (default), overdispersion will be taken into account (see Details), 
#' 	if overdisp is FALSE binomial distribution will be used. 
#' @param disp.param - to use external disperison parameters
#' @param .parallel, .progress parameters for plyr::alply 
#' 
#' @return matrix with rows for each fit, firts column contains overdispersion, remainings contain p-values for each term
#' 
#' @details if overdisp is TRUE, quasibinomial distribution will be used if overdispersion is more than 1,
#'  binomial distibution will be used instead. 
#'  So, switch it off if you don't like to account for overdispersion or if you don't have replicates.
#' 
#' @examples
#' alt.pv.bq = calcSAPvalue(alt.glm)
#' 
#' @import plyr
#' @export
calcSAPvalue = function(data.glm,overdisp=TRUE,disp.param=NULL,.parallel=FALSE,.progress='none'){
	if(!(overdisp %in% c(TRUE,FALSE)))
		stop(paste('Wrong parameter overdisp Only TRUE or FALSE are allowed'))
	term.names = attr(data.glm,'term.labels')
	term.inx = 2:(1+length(term.names))
	
	res = laply(1:length(data.glm),function(i){
		getPvForGLM(data.glm[[i]],overdisp,disp.param[i],term.names,names(data.glm)[i])
		},.parallel=.parallel,.progress=.progress)
	rownames(res) = names(data.glm)
	colnames(res) = c('overdispersion',term.names)
	res[is.na(res)] = NA #convert NaNs to NA
	res
}

getPvForGLM = function(g,overdisp,disp.param,term.names,sid){
	r = rep(NA,length(term.names)+1)
	if(is.na(g)[1])
		return(r)
	#set disp parameter
	if(is.null(disp.param)){
		r[1] = summary(g)$dispersion #sum(residuals(g, type = 'pearson')^2)/g$df.residual
		if(g$df.residual == 0){
			r[1] = NA
			if(overdisp)
				warning(paste('Segment ',sid,': cannot account for overdispersion in absense of replicates. Either use replicates or set overdisp parameter to FALSE',sep=''))
		}
	}else
		r[1] = disp.param
	disp = max(1,r[1],na.rm=TRUE)
	if(!overdisp) disp = 1
	a = tryCatch(anova(g,test='Chisq',dispersion=disp),error=function(e){return(NULL)})
	if(!is.null(a)){
		r[-1] = a[2:(1+length(term.names)),5]
	}
	r
}

#' Cluster segments by splicing similarity
#' 
#' Performs simple \code{\link{hcluts}} and plots clusters as boxplot. First normalize mean and variance to 0 and 1 (see norm.groups parameter) and then used 1-cor distance.
#'
#' @param data object of sajr class, produced by \code{\link{loadSAData}}. You may want to use only significan ones (see \code{\link{calcSAPvalue}})
#' @param k desired numer of clusters
#' @param rows number of rows in multiplot
#' @param cols number of columns in multiplot
#' @param col graphical parameter
#' @param dist type of distance to be used. Either cor (1-cor will be used) or argument of dist function
#' @param norm.groups factors that denotes in which groups mean should be normolised to 0
#' @param norm.sd whether variance should be normolized to 1
#' @param outline whether outline values shold be plotted
#' 
#' @return clustering for segments. Segments with zero variance (after mean normalization) will omitted
#' 
#' @examples
#' t = clustSegs(age[,c(1,3,2,4)],6,2,3,col=c('orange','gray','orange','gray'),dist='cor',norm.groups=c(1,2,1,2))
#' 
#' @export
clustSegs = function(data,k,rows,cols,col,dist='cor',norm.groups=NULL,norm.sd=F,outline=T,plotPCAandTree=F){
	i = data$ir
	for(g in unique(norm.groups)){
		inx = norm.groups == g
		i[,inx] = sweep(i[,inx],1,apply(i[,inx],1,mean,na.rm=T),'-')
	}
	if(norm.sd)
		i = sweep(i,1,apply(i,1,sd,na.rm=T),'/')
	if(dist=='cor')
		d = as.dist(1-cor(t(i),use='p'))
	else
		d = dist(i,method=dist)
#	d = dist(i)
	tr = hclust(d)
	cl_ = cutree(tr,k=k)
	tab = table(cl_)
	tab.o = order(tab,decreasing=T)
	cl = cl_
	for(c in 1:length(tab))
		cl[cl_ == names(tab)[tab.o[c]]] = c
	#plot tree and PCA
	if(plotPCAandTree){
		mds = cmdscale(d,2)
		plot(mds[,1],mds[,2],pch=16,col=cl,xlab='Coordinate 1',ylab='Coordinate 2',main='MDS')
		hcl = as.dendrogram(tr)
		hcl = dendrapply(hcl, colLab,clusters=cl,1:k)
		plot(hcl,main='Tree',xlab='',ann=F,col=cl,leaflab='none')
	}
	#plot profiles
	par(mfrow=c(rows,cols))
	for(j in 1:k){
		t = i[cl == j,]
		boxplot(t,main=paste('c',j,':',dim(t)[1],':',length(unique(data$seg$gene_id[cl == j])),sep=''),col=col,outline=outline)
		abline(h=0,col='gray',lty=2)
	}
	cl
}

#' Colors three by clusters
#'
#' @param n dendrogram
#' @param clusters clustering of leafs
#' @param cols colors of clusters
colLab = function(n,clusters,cols) {
	leafs = order.dendrogram(n)
	leaf.cl = unique(clusters[leafs])
	if(length(leaf.cl)==1){
		attr(n,'edgePar') = list(col=cols[leaf.cl])
	}
	n
}

#' Removes some kinds of segments from alternatives
#' 
#'
#' @param a output of \code{\link{makeAlts}} method
#' @param filter.introns whether retained introns should be removed  
#' @return alternatives in the same format. Segments specified by parameters are removed. Alternatives without segments are removed from results.
#'  
#' @export 
filterAlts = function(a,filter.introns=F){
	if(!filter.introns && !filter.last.first)
		return(a)
	for(i in 1:dim(a)[1]){
		ss = unlist(strsplit(a[i,'segs'],';'))
		is = unlist(strsplit(a[i,'ints'],';'))
		good = rep(T,times=length(ss))
		if(filter.introns){
			good = !(ss %in% is)
		}
		a[i,'segs'] = paste(ss[good],collapse=';')
	}
	a[a$segs!='',]
}

#' Plot schematic representation of all alternatives
#' 
#' Alternatives are sortered by number of their occurence in genome
#' Segments are shown as rects, junctions are shown as arcs, see \code{\link{plotAlt}} for details.
#'
#' @param a output of \code{\link{makeAlts}} method
#' @param min.cnt plot only alternatives that meats as many times
#' @param to.plot plot not more as to.plot alternatives
#'  
#' @export
plotAllAlts = function(a,min.cnt=0,to.plot=Inf,...){
	a=table(paste(a$ints,a$segs))
	a=sort(a,decreasing=T)
	a = a[a>=min.cnt]
	if(to.plot<length(a))
		a = a[1:to.plot]
	for(i in 1:length(a)){
		t = unlist(strsplit(names(a)[i],' '))
		plotAlt(t[1],t[2],main=a[i],...)
		box(which='figure')
	}	
}

#' Plot schematic representation of alternative
#' 
#' Segments are shown as rects, junctions are shown as arcs
#'
#' @param ints character with list of introns in form of: 0-1;0-3;2-3 (see \code{\link{makeAlts}} for details)
#' @param segs character with list of segments in form of: 0-1;0-3;2-3 (see \code{\link{makeAlts}} for details)
#'  
#' @export
plotAlt = function(ints,segs,main='',col.exn='green',col.alt='yellow',col.int='blue',col.junc='red',lwd.junc=2,int.part=0.5,lf.part=0.3,plot.mid.line=T){
	is.ret.int = unlist(strsplit(segs,';')) %in%  unlist(strsplit(ints,';'))
	ints = parseIntSeg(ints)
	segs = parseIntSeg(segs)
	segs[is.na(segs[,1]),1] = segs[is.na(segs[,1]),2]-lf.part
	segs[is.na(segs[,2]),2] = segs[is.na(segs[,2]),1]+lf.part
	ymax=1/(1-int.part)
	plot(1,t='n',xaxt='n',yaxt='n',bty='n',xlab='',ylab='',main=main,xlim=c(-lf.part,max(ints)+lf.part),ylim=c(0,ymax))
	if(plot.mid.line)
		lines(x=c(0,max(ints)),y=c(0.5,0.5),col='black',lwd=1)
	#plot constitutve exons
	rect(-lf.part,0,0,1,col=col.exn,border=NA)
	rect(max(ints),0,max(ints)+lf.part,1,col=col.exn,border=NA)
	#retained introns
	if(sum(!is.ret.int)>0)
		rect(segs[!is.ret.int,1],0,segs[!is.ret.int,2],1,col=col.alt,border=NA)
	#alternative
	if(sum(is.ret.int)>0)
		rect(segs[is.ret.int,1],0,segs[is.ret.int,2],1,col=col.int,border=NA)
	for(i in 1:dim(ints)[1])
		plotParabola(c(ints[i,1],(ints[i,1]+ints[i,2])/2,ints[i,2]),c(0.5,ymax,0.5),lwd=lwd.junc,col=col.junc)
	sites.pos = unique(c(ints[,1],ints[,2]))
	segments(sites.pos,0,sites.pos,1,col='black',lwd=1)
}

#' Assign segments to splicing alternatives
#'
#' Cuts genes by constitutive segments into regions that 
#' can be spliced alternatively
#' @param a alternatives, output of \code{\link{makeAlts}}
#' @param s segment annotation (first element of output of \code{\link{loadSAData}} function)
#' 
#' @return segments with alt.id field
#'  
#' @export
assignSegs2Alts = function(a,s){
	s$alt.id = NA
	ord = rownames(s)
	s = s[order(s$gene_id,s$start),]
	a = a[order(a$gene_id,a$start),]
	si=1
	ai=1
	while(si <= nrow(s) & ai <= nrow(a)){
		cat('\r',si,ai,nrow(s),nrow(a))
		if(s$gene_id[si] == a$gene_id[ai] & s$start[si] >= a$start[ai] & s$stop[si] <= a$stop[ai]){
			s$alt.id[si] = ai
			si = si + 1
		}else if(s$gene_id[si] < a$gene_id[ai] || (s$gene_id[si] == a$gene_id[ai] & s$start[si] < a$start[ai])){
			si = si + 1
		}else
			ai = ai + 1
	}
	s$alt.id = rownames(a)[s$alt.id]
	s[ord,]
}

#' Splits genes into splicing alternatives
#'
#' Cuts genes by constitutive segments into regions that 
#' can be spliced alternatively
#' @param seg segment annotation (first element of output of \code{\link{loadSAData}} function)
#' @param ann.gff name of gff file
#' @param remove.exn.ext specifies whether alternative TSS/polyA that aren't supplied by alternative splicing (i.e. just an extension of other exons) should be removed
#' 
#' @return data.frame with all found alternatives. 
#' sites is types of splicing sites: [a]cceptor and [d]onor
#' segs is list of alternative segments separated by ';'. Each segment is denoted by ranks of splice sites that form segment's boundaries, rank == NA correspons to alternative TSS or polyA
#' ints the same as segs but for introns
#'  
#' @export
makeAlts = function(seg,ann.gff,remove.exn.ext=F){
	require(data.table)
	gff = read.table(ann.gff,comment.char='#',sep="\t",as.is=T)
	seg = seg[seg$position != 'ONLY',]
	if(remove.exn.ext){
		seg = seg[!(seg$sites %in% c('sa','de')),]
		seg$type[seg$position %in% c('FIRST','LAST')] = 'EXN'
		
		alt.start = table(seg$gene_id[seg$position=='FIRST'])
		alt.start = names(alt.start)[alt.start>1]
		
		alt.stop = table(seg$gene_id[seg$position=='LAST'])
		alt.stop = names(alt.stop)[alt.stop>1]

		seg$type[(seg$gene_id %in% alt.start) & seg$position=='FIRST'] = 'ALT'
		seg$type[(seg$gene_id %in% alt.stop) & seg$position=='LAST'] = 'ALT'
	}
	seg = seg[,c("gene_id",'start','stop','strand','type','position')]
	gff = gff[gff[,3] == 'intron',]
	int = as.data.frame(matrix(ncol=4,nrow = dim(gff)[1]))
	colnames(int) = c("gene_id",'start','stop','strand')
	int[,1] = sapply(gff[,9],getAttr,name='gene_id=')
	int[,2:3] = gff[,c(4,5)]
	int[gff[,7] == '+',4] = 1
	int[gff[,7] == '-',4] =-1
	rm(gff)
	
	int = split(int,factor(int$gene_id))
	seg = split(seg,factor(seg$gene_id))
	
	alts = new.env()
	gid.names = intersect(names(seg),names(int))
	j = 1
	for(gid in gid.names){
		cat('\r',j,'from',length(gid.names),'      ')
		j = j + 1
		t = makeGAlts(seg[[gid]],int[[gid]])
		if(length(t)>0)
			for(i in 1:length(t))
				alts[[paste('t',length(alts),sep='')]] = t[[i]]
	}
	
	r = as.data.frame(rbindlist(as.list(alts)))
	r = lapply(split(r,r$gene_id),function(g){g$alt.id = paste(g$gene_id,1:nrow(g),sep='.a');g})
	r = as.data.frame(rbindlist(r))
	rownames(r) = r$alt.id
	r$alt.id = NULL
	r
}

#' Finds all alternative within gene
#'
#' Cuts gene by constitutive exons
#'
#' @param segs data.frame with all gene segments
#' @param ints data.frame with all gene introns
#' 
#' @return list of dataframes with all gene alternatives
makeGAlts = function(segs,ints){
	strand = segs$strand[1]
	if(strand==-1){
		segs[,2:3] = -segs[,3:2]
		ints[,2:3] = -ints[,3:2]
	}
	segs = segs[order(segs[,2]),]
	ints = ints[order(ints[,2],ints[,3]),]
	sites.coors = unique(c(paste(ints[,2]-1,'d',sep=''),paste(ints[,3],'a',sep='')))
	sites.types = substr(sites.coors,nchar(sites.coors),nchar(sites.coors))
	sites.coors = as.numeric(substr(sites.coors,1,nchar(sites.coors)-1))
	o = order(sites.coors)
	sites.coors = sites.coors[o]
	sites.types = sites.types[o]
	#change coordinates to site ranks
	ints.coor = ints[,2:3]
	segs.coor = segs[,2:3]
	ints[,2] = findInterval(ints.coor[,1]-1,sites.coors)
	ints[,3] = findInterval(ints.coor[,2],sites.coors)
	segs[,2] = findInterval(segs.coor[,1]-1,sites.coors)
	segs[,3] = findInterval(segs.coor[,2],sites.coors)
	segs[segs[,2]==0,2] = 1
	segs[segs[,3]==0,3] = 1
	segs[sites.coors[segs[,2]] != segs.coor[,1]-1,2] = NA
	segs[sites.coors[segs[,3]] != segs.coor[,2],3] = NA
	alts = list()
	in.alt = F
	first.alt = -1
	for(i in 1:dim(segs)[1]){
		if(segs$type[i]!='EXN' && !in.alt){
			first.alt = i
			in.alt = T		
		}
		if(in.alt && (segs$type[i]=='EXN' || i==dim(segs)[1])){
			#make alternative
			finx = max(1,first.alt-1)
			start.coor = ifelse(first.alt==1,segs.coor[1,1],segs.coor[finx,2]+1)
			#if it isn't constitutive it should be last
			stop.coor = ifelse(segs$type[i]!='EXN',segs.coor[dim(segs)[1],2],segs.coor[i,1]-1)
			aranks = c(segs[finx,3],segs[i,2])
			asegs = segs[first.alt:ifelse(segs$type[i]=='EXN',i-1,i),]
			asegs = paste(paste(asegs[,2]-aranks[1],asegs[,3]-aranks[1],sep='-'),collapse=';')
			
			aints = mfi(aranks,ints[,2])
			aints = ints[aints[1]:aints[2],]
			aints = paste(paste(aints[,2]-aranks[1],aints[,3]-aranks[1],sep='-'),collapse=';')
			r = list(gene_id=segs[1,1],start=start.coor,'stop'=stop.coor,
							 sites=paste(sites.types[aranks[1]:aranks[2]],collapse=''),
							 segs=asegs,ints=aints)
			if(strand==-1){
				tmp = r
				r[['start']] = -tmp[['stop']]
				r[['stop']] = -tmp[['start']]
			}
			alts[[length(alts)+1]] = as.data.frame(r,stringsAsFactors=F)
			first.alt = -1
			in.alt = F
		}
	}
	alts
}


#' Finds all elements from within [x[1],x[2]]
#'
#' @param x numeric vector with at least two elements, x[1] have to be not more than x[2]
#' @param y sorted numeric vector
#' 
#' @return first and last indexes of elements in v that are within [x[1],x[2]]
mfi = function(x,v){
	t = findInterval(x,v)
	while(T){
		if(t[1]>1 && v[t[1]-1] == x[1])
			t[1] = t[1] -1
		else
			break
	}
	while(T){
		if(t[2] < length(v) && v[t[2]+1] == x[2])
			t[2] = t[2]+1
		else
			break
	}
	t
}

#' Parse string representation of list of intervals 
#'
#' @param t list of intervals in form of: 0-1;3-6;7-NA
#' 
#' @return numeric matrix with two columns
parseIntSeg=function(t){
	if(nchar(t)==0)
		return(matrix(ncol=2,nrow=0))
	t = unlist(strsplit(t,';'))
	m = matrix(nrow=length(t),ncol=2)
	for(i in 1:length(t)){
		x = unlist(strsplit(t[i],'-'))
		x[x=='NA'] = NA
		m[i,] = as.numeric(x)
	}
	m
}

#' Plot parabola specified by 3 points
#'
#' @param x vector with at leqast 3 values
#' @param y vector with at leqast 3 values
#' @param n number of points to approximate curve
plotParabola = function(x,y,n=30,...){
	x1 = x[2] - x[1]
	x2 = x[3] - x[1]
	y1 = y[2] - y[1]
	y2 = y[3] - y[1]
	c = y[1]
	den = x1*x2*(x2-x1)
	a = (y2*x1-y1*x2)/den
	b = (y1*x2*x2-y2*x1*x1)/den
	c = c - b*x[1]+a*x[1]*x[1]
	b = b - 2*a*x[1]
	p = function(x){a*x*x + b*x +c}
	xi = (0:n)/n*(max(x)-min(x))+min(x)
	yi = p(xi)
	lines(xi,yi,...)
}

#' Makes a summary about common splicing alternnative types
#'
#' @param a output of \code{\link{makeAlts}}
#' @param print whether statistic should be printed
#'  
#' @return vector of counts of simple alternatives
#'  
#' @export
alt.summary = function(a,print=T){
	a = paste(a$ints,a$segs)
	s = function(al,nm){
		cnt = sum(a %in% al)
		if(print)
			cat(nm,': ',cnt,' (',round(cnt/length(a)*100,2),'%)\n',sep='')
		cnt
	}
	if(print)
		cat('Total: ',length(a),'\n',sep='')
	names = c('Retained introns','Cassette exon','Alternative donor','Alternative acceptor','Complex')
	alts = list('0-1 0-1','0-1;0-3;2-3 1-2','0-2;1-2 0-1','0-1;0-2 1-2')
	alts[[length(alts)+1]] = setdiff(unique(a),unlist(alts))
	cnts = numeric(length(names)+1)
	names(cnts) = c('Total',names)
	cnts[1] = length(a)
	for(i in 1:length(names))
		cnts[i+1] = s(alts[[i]],names[i])
	invisible(cnts)
}

#' calculates z-scores
#' 
#' @param d data
#' @param MARGIN
#'  
#' @export
normRows=function(d,MARGIN=1){
	d = sweep(d,MARGIN,apply(d,MARGIN,mean,na.rm=T),'-')
	sweep(d,MARGIN,apply(d,MARGIN,sd,na.rm=T),'/')
}

#' Make sequential palette with n colors 
#' 
#' @param col colors to be used to construct palette
#' @param n desired number of colors
#'  
#' @export
getPal=function(col=c('blue','white','red'),n=200){
	r = character(n)
	if(n == 1){
		i = floor(length(col)/2)
		return(.getPal(col[i],col[i+1],0.2))
	}
	for(i in 0:(n-1)){
		cinx = i/(n-1)*(length(col)-1)+1
		if(cinx == floor(cinx))
			r[i+1] = col[cinx]
		else{
			rate = cinx-floor(cinx)
			cinx = floor(cinx)
			r[i+1] = .getPal(col[cinx],col[cinx+1],1-rate)
		}
	}
	r
}

.getPal = function(c1,c2,r){
	c1 = t(col2rgb(c1,alpha=TRUE)/255)
	c2 = t(col2rgb(c2,alpha=TRUE)/255)
	return(rgb((c1*r+c2*(1-r)),alpha=r*c1[4]+(1-r)*c2[4]))
}


#' plots symmetrical correlation heatmap
#' 
#' @param data data (set of columns) to be used to plot heatmap. 
#' @param cor correlation matrix to be used instead of data
#' @param method correlation method
#' @param main graphical parameter
#' @param norm specifies whether rows should be normolized
#' @param zeroIncenter specifies whether color cols[length(cols)/2] should correspons to zero correlation
#' @param cols colors used to plot heatmap
#' @param reorder heatmap parameter
#' @param hclust.method 
#' @param ... other parameters to be passed to heatmap function
#'  
#' @export
plotCorrHM = function(data=NULL,cor=NULL,method='sp',main='',norm=T,zeroIncenter=norm,cols=NULL,reorder=NULL,hclust.method="complete",...){
	if(is.null(cor)){
		if(norm)
			data = normRows(data)
		c = cor(data,use='pair',method=method)
	}else
		c = cor
	if(is.null(cols))
		cols = getPal(col=c('blue','white','red'),n=200)
	if(zeroIncenter){
		r = round(range(c)*length(cols)/2)+length(cols)/2
		cols=cols[r[1]:r[2]]
	}
	if(is.null(reorder))
		reorder = 1:ncol(c)
	heatmap(c,symm=T,col=cols,main=paste(main,' (',dim(data)[1],')',sep=''),Rowv=reorder,reorderfun=function(d,w){reorder(d,w,agglo.FUN=mean)},
					distfun=function(d){as.dist(1-d)},hclustfun=function(d){hclust(d,method=hclust.method)},...)
	cx = grconvertX(c(0.03,0.13), from='ndc',to='u')
	cy = grconvertY(c(0.95,0.75), from='ndc',to='u')
	lg = round((0:3)/3*(max(c)-min(c))+min(c),2)
	lg[length(lg)] = '1.00'
	plotrix::color.legend(cx[1],cy[1],cx[2],cy[2],legend=lg,rect.col=cols,gradient="y",align='rb')
	invisible(c)
}
