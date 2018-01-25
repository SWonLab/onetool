onetool.plot.pedigree <- function(v, fn) {
	ok <- TRUE
	if (class(try(library(kinship2), T)) == 'try-error') {
		if (!file.exists(Sys.getenv("R_LIBS_USER"))) dir.create(Sys.getenv("R_LIBS_USER"), recursive=T)
		install.packages("kinship2", repos="http://healthstat.snu.ac.kr/CRAN", lib = Sys.getenv("R_LIBS_USER"))
		.libPaths(c( .libPaths(), Sys.getenv("R_LIBS_USER")))
		if (class(try(library(kinship2), T)) == 'try-error') {
			cat("Failed to load [kinship2] package, pedigree will not be drawn\n")
			ok <- FALSE
		}
	}
	if (ok == TRUE) {
		pdf(fn, width=10)
		affected <- as.numeric(v[,6])
		ped <- pedigree(id=v[,2], dadid=v[,3], momid=v[,4], sex=as.numeric(v[,5]), famid=v[,1], affected=rep(0,length(affected)))
		fid <- unique(ped$famid)
		plot(ped[fid[1]])
		mtext(paste("Family [", fid[1], "], # = ", sum(ped$famid==fid[1]), sep=""), line=1, cex=1.5)
		for (i in 2:length(fid)) {
			plot(ped[fid[i]])
			mtext(paste("Family [", fid[i], "], # = ", sum(ped$famid==fid[i]), sep=""), line=1, cex=1.5)
		}
		dev.off()
	}
}
