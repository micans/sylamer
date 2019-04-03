
fname <- "out.parse-syloscope.pl"

library(R.utils)
library(evd)

args <- R.utils::commandArgs(asValues=TRUE)

fname          <- args$sylsum
qtile          <- as.numeric(args$quantile)
filter         <- as.numeric(args$filter)
ftable         <- args$table

if (is.null(fname) || is.null(qtile) || is.null(ftable)) {
   stop("need --sylsum=fname --quantile=fraction --table=eightmer-table")
}

tb <- read.table(fname, as.is=T, row.names=1, header=T)
mirna <- read.table(ftable, as.is=T, row.names=NULL, header=T)
rownames(mirna) <- mirna$eight1A

# tbe <- fgev(tb$maxabs[tb$maxabs > quantile(tb$maxabs, qtile) & tb$maxabs < quantile(tb$maxabs, 1-qtile)])
tbe <- fgev(tb$maxabs[tb$maxabs < quantile(tb$maxabs, 1-qtile)])
eee <- as.list(tbe$estimate[c("loc", "scale", "shape")])

##  options(warn=-1)
##  ks <- ks.test(tb$maxabs[tb$maxabs > quantile(tb$maxabs, qtile) & tb$maxabs < quantile(tb$maxabs, 1-qtile)],
##  "pgev", loc=eee$loc, scale=eee$scale, shape=eee$shape)
##  options(warn=0)

tb$pval <- apply(tb, 1, function (x) { pgev(as.numeric(x["maxabs"]), eee$loc, eee$scale, eee$shape, lower=F) })
tb$handle <- mirna[rownames(tb), "handle"]
tb$lsomirna <- mirna[rownames(tb), "lsomirna"]

if (filter > 0) {
   tb <- tb[tb$pval <= filter,]
}
tb <- tb[order(tb$pval),]
write.table(tb, file="", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

