
my.invnorm = function(x) {
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
return(res)
}

fisherp <- function(x) {
    if (any(is.na(x))) {
        warning("Some p-values missing; removing these.")
        x <- stats::na.omit(x)
    }
    if (any(x<=0 | x>1)) stop("P-values must be >0 and <=1.")
    if (length(x)<2) stop("Must have at least two valid p-values.")
    df <- 2*length(x)
    fisherp <- stats::pchisq( -2*sum(log(x)), df, lower.tail=FALSE)
    return(fisherp)
}

run_diff_mod <- function(mm) {
	all_ps = NULL
	ps_dat = list()
	genes = unique(mm$gene)
	pps = vector()
	pps1 = vector()
	fpps = vector()
	for(x in 1:length(genes)) {
		m1 = mm[mm$gene==genes[x],]
		sams = unique(m1$sample_name)
		ps = vector()
		ks = vector()
		mw = vector()
		mw_est = vector()
		poses = vector()
		for(y in 1:length(sams)) {
			m2 = m1[m1$sample_name==sams[y],]
			# maybe do this - restrict to site seen in both haplotypes
			h1 = m2[m2$hap==1,]
			h2 = m2[m2$hap==2,]
			m2 = m2[m2[,2]%in%h1[,2],]
			m2 = m2[m2[,2]%in%h2[,2],]
			mr = m2$mod_ratio
			mr[mr>0.5] = 1
			mr[mr<=0.5] = 0
			t = table(m2$hap, mr)
			poses[y] = max(m2$end)-((max(m2$end) - min(m2$start))/2)
			if(dim(t)[1]==2 && dim(t)[2]==2) {
				ps[y] = fisher.test(t)$p.value
				k = kruskal.test(mod_ratio ~ hap, data = m2)$p.value
				ks[y] = k
				wil = wilcox.test(mod_ratio ~ hap, data = m2, conf.int = TRUE)
				mw[y] = wil$p.value
				mw_est[y] = wil$estimate
				#, conf.int = TRUE
			} else {
				ps[y] = NA
				ks[y] = NA
				mw[y] = NA
			}
			
		}
		ps_dat[[x]] = cbind(ps, ks, mw, mw_est, genes[x], sams, unique(m1$chr), poses)
		all_ps = rbind(all_ps, ps)
	ps[ps>0.9999] = NA
	pps[x] = mean(ps, na.rm=T)
	pps1[x] = min(ps, na.rm=T)
	fpps[x] = fisherp(ps)
	print(x)
	}
	dat = data.frame(do.call(rbind, ps_dat))
	colnames(dat) = c("p", "k", "w", "w_est", "g", "s", "c", "pos")
return(dat)
}

plot_den <- function(snp, g = "LMF1") {
	sg = snp[snp$SNP==g,]
	m1 = mm[mm$gene==g,]
	h1 = m1[m1$hap==1,]
	h2 = m1[m1$hap==2,]
	m1 = m1[m1[,2]%in%h1[,2],]
	m1 = m1[m1[,2]%in%h2[,2],]
	mb = m1[m1$sample_name==sg$SAM[which.min(sg$K)],]
	mh = m1[m1$sample_name==sg$SAM[which.max(sg$K)],]
	mb_h1 = mb[mb$hap==1,]
	mb_h2 = mb[mb$hap==2,]
	mb = mb[mb[,2]%in%mb_h1[,2],]
	mb = mb[mb[,2]%in%mb_h2[,2],]
	mh_h1 = mh[mh$hap==1,]
	mh_h2 = mh[mh$hap==2,]
	mh = mh[mh[,2]%in%mh_h1[,2],]
	mh = mh[mh[,2]%in%mh_h2[,2],]
	par(mfrow=c(1,2))
	d_high1 = density(mb$mod_ratio[mb$hap==1])
	d_high2 = density(mb$mod_ratio[mb$hap==2])
	tit1 = paste("Sample:", sg$SAM[which.min(sg$K)], " | -log10(P):", round(-log10(sg$K[which.min(sg$K)]), 3))
	plot(d_high1$x, d_high1$y, type="l", col="blue", ylim=c(0, max(c(max(d_high1$y), max(d_high2$y)))), main=tit1, xlab="Mod Ratio", ylab="Density")
	matplot(d_high2$x, d_high2$y, type="l", col="red", add=T)
	legend("topleft", c("Hap 1", "Hap 2"), col=c("blue", "red"), lty="solid", lwd=2, bty="n")
	d_high1 = density(mh$mod_ratio[mh$hap==1])
	d_high2 = density(mh$mod_ratio[mh$hap==2])
	tit2 = paste("Sample:", sg$SAM[which.max(sg$K)], " | -log10(P):", round(-log10(sg$K[which.max(sg$K)]), 3))
	plot(d_high1$x, d_high1$y, type="l", col="blue", ylim=c(0, max(c(max(d_high1$y), max(d_high2$y)))), main=tit2, xlab="Mod Ratio", ylab="Density")
	matplot(d_high2$x, d_high2$y, type="l", col="red", add=T)
	legend("topleft", c("Hap 1", "Hap 2"), col=c("blue", "red"), lty="solid", lwd=2, bty="n")
}

plot_den_samples <- function(snp, g = "LMF1", sam1, sam2) {
	sg = snp[snp$SNP==g,]
	m1 = mm[mm$gene==g,]
	h1 = m1[m1$hap==1,]
	h2 = m1[m1$hap==2,]
	m1 = m1[m1[,2]%in%h1[,2],]
	m1 = m1[m1[,2]%in%h2[,2],]
	mb = m1[m1$sample_name==sam1,]
	mh = m1[m1$sample_name==sam2,]
	mb_h1 = mb[mb$hap==1,]
	mb_h2 = mb[mb$hap==2,]
	mb = mb[mb[,2]%in%mb_h1[,2],]
	mb = mb[mb[,2]%in%mb_h2[,2],]
	mh_h1 = mh[mh$hap==1,]
	mh_h2 = mh[mh$hap==2,]
	mh = mh[mh[,2]%in%mh_h1[,2],]
	mh = mh[mh[,2]%in%mh_h2[,2],]
	par(mfrow=c(1,2))
	d_high1 = density(mb$mod_ratio[mb$hap==1])
	d_high2 = density(mb$mod_ratio[mb$hap==2])
	tit1 = paste("Sample:", sam1, " | -log10(P):", round(-log10(sg$W[sg$SAM==sam1]), 3))
	plot(d_high1$x, d_high1$y, type="l", col="blue", ylim=c(0, max(c(max(d_high1$y), max(d_high2$y)))), main=tit1, xlab="Mod Ratio", ylab="Density")
	matplot(d_high2$x, d_high2$y, type="l", col="red", add=T)
	legend("topleft", c("Hap 1", "Hap 2"), col=c("blue", "red"), lty="solid", lwd=2, bty="n")
	d_high1 = density(mh$mod_ratio[mh$hap==1])
	d_high2 = density(mh$mod_ratio[mh$hap==2])
	tit2 = paste("Sample:", sam2, " | -log10(P):", round(-log10(sg$W[sg$SAM==sam2]), 3))
	plot(d_high1$x, d_high1$y, type="l", col="blue", ylim=c(0, max(c(max(d_high1$y), max(d_high2$y)))), main=tit2, xlab="Mod Ratio", ylab="Density")
	matplot(d_high2$x, d_high2$y, type="l", col="red", add=T)
	legend("topleft", c("Hap 1", "Hap 2"), col=c("blue", "red"), lty="solid", lwd=2, bty="n")
}

main_manplot <- function(snp, cut, pin = 10) {
	t = table(snp$SNP[snp$W<cut])
	g1 = ggplot(snp[snp$SNP%in%names(t[t>=pin]),], aes(x= fct_reorder(SNP, as.numeric(CHR)), y=-log10(W), color=factor(CHR))) +    
	geom_point()+ theme_bw()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	geom_hline(yintercept=-log10(cut))+xlab("Gene | Sample")+ylab("-log10(P)")+
	guides(color = guide_legend(title = "Chromosome"))
}

multi_plot <- function(snp, outfile) {
	gs_list = list()
	for(x in 1:22) {
		gs_list[[x]] = ggplot(snp[snp$CHR==x,], aes(x= fct_reorder(SNP, as.numeric(CHR)), y=-log10(W), color=factor(SAM))) +    
		geom_point()+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=-log10(cut))+xlab("Gene")+
		guides(color = guide_legend(title = "Sample")) +ggtitle(paste("Chromosome:", x))+ theme(legend.position = "none")

	}
	pdf(outfile, height=20, width=40)
	ggarrange(plotlist =gs_list, ncol=6, nrow=4)
	dev.off()
}

plot_est_p <- function(dat, gen="TERT") {
	plot(my.invnorm(dat$w_est[dat$g==gen]), -log10(dat$w[dat$g==gen]), col=(-log10(dat$w[dat$g==gen]) > -log10(cut) )+1, pch=20, ylab="-log10(p)", xlab="MW Estimate", cex=1.2, main=gen)
	abline(h=-log10(cut), lty='dashed')
}

get_data <- function(cpg_file) {
	dat = read.table(, header=T)
	dat$w_est_n = my.invnorm(dat$w_est)
	dat = dat[complete.cases(dat),]
	snp = data.frame(CHR=as.numeric(as.character(dat$c)), BP=as.numeric(as.character(dat$pos)), P=as.numeric(as.character(dat$p)), K=as.numeric(as.character(dat$k)), W=as.numeric(as.character(dat$w)), W_est=as.numeric(as.character(dat$w_est)), SNP=dat$g, SAM=dat$s)
	snp =snp[complete.cases(snp),]
return(snp)
}


library(ggplot2)
library(forcats)
library(ggpubr)

#### MAIN

cpg_file="samples_cpg_pileup_66_5mC_gene_assoc.txt"
snp = get_data(cpg_file)
cut = 0.01 / (length(unique(snp$SNP))*length(unique(snp$SAM))) 
cut = 0.01 / nrow(snp)


#### PLOTS

#main_manplot(snp, cut)
#multi_plot(snp, "hap_diff_meth_294_clinically_important_genes_66.pdf")

#plot_den(snp)
#plot_den(snp, "HYAL1")
#plot_den(snp, "DEAF1")
#plot_den(snp, "KCNV2")

#plot_est_p(dat, "NLRP2")
#plot_den_samples(snp, "NLRP2", "HG00323", "HG00242")


