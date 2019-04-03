
# library(R.utils)
args <- R.utils::commandArgs(asValues=TRUE)

# args <- commandArgs(asValues=TRUE)

pdfout <- "t.pdf"
fnsyl <- ""

title    <- "Sylamer%landscape"
selectedwords <- "Selected words"
secondlegend <- ""

topOligos <- 3
lowOligos <- 3
oligosChosen <- c()
n_max <- 500
ymin <- 0
ymax <- 0
dyad <- ""
prefix_max <- 0
fnannot <- ""
fntable <- ""
fngenes <- ""
species <- ""
annotfilter <- "dddd,six3"
batch <- TRUE
hili_cons_mirna <- FALSE
bonferroni <- FALSE
# annotfilter=""


if (!is.null(args[["data"]])) {
   fnsyl <- args$data
}
if (!is.null(args[["pdf"]])) {
   pdfout <- args$pdf
}
if (!is.null(args[["dyad"]])) {
   dyad <- args$dyad
}
if (!is.null(args$top)) {
   topOligos <- as.integer(args$top)
}
if (!is.null(args$ymin)) {
   ymin <- as.integer(args$ymin)
}
if (!is.null(args$species)) {
   species <- args$species
}
if (!is.null(args$interactive)) {
   batch <- FALSE
}
if (!is.null(args$table)) {
   fntable <- args$table
}
if (!is.null(args$ymax)) {
   ymax <- as.integer(args$ymax)
}
if (!is.null(args$bot)) {
   lowOligos <- args$bot
}
if (!is.null(args$bonferroni)) {
   bonferroni <- TRUE
}
if (!is.null(args$add)) {
   oligosChosen <- unlist(strsplit(args$add, ","))
   print(oligosChosen)
}
if (!is.null(args$"annot")) {
   fnannot <- args$"annot"
}
if (!is.null(args$"hili-all-cons")) {
   hili_cons_mirna <- TRUE
}
if (!is.null(args$"add-colours")) {
   colorsChosen <- unlist(strsplit(args$"add-colours", ","))
}
if (!is.null(args[["title"]])) {
   title <- args$title
}
               # weird name for this option.
if (!is.null(args[["add-legend"]])) {
   selectedwords <- args$"add-legend"
}
               # e.g. for the top words of the EVD list.
if (!is.null(args[["second-legend"]])) {
   secondlegend <- args$"second-legend"
}
if (!is.null(args[["prefix-max"]])) {
   prefix_max <- args$"prefix-max"
}
if (!is.null(args[["genes"]])) {
   fngenes <- args$genes
}

thegenes <- c()

if (nchar(fngenes) > 0) {
   thegenes <- read.table(fngenes, header=F, row.names=1)
}

rt <- function (s, mode) {
   have_colnames <- regexpr("c", mode) >= 0
   have_rownames <- regexpr("r", mode) >= 0
   check_names <- regexpr("n", mode) >= 0
   if (have_rownames) read.table(s, check.names=check_names, sep="\t", header=have_colnames, row.names=1, as.is=T, quote="", comment.char="") else
   read.table(s, check.names=check_names, sep="\t", header=have_colnames, as.is=T, quote="", comment.char="")
}


allstrings <- function (nmer, k)
{  ret <- unlist(lapply(1:(nchar(nmer)-k+1), function(x) { substring(nmer, x, k+x-1) }))
   if (k == 6 && nchar(nmer) == 8 && regexpr("six3", annotfilter) > 0) {
      ret <- ret[2:3]
   }
   ret
}


   #  y is a list of strings; each string a comma-separated list of microRNAs.

normalise_mirna_names <- function (y) {

         # y1 is simply all microRNA names
   y1 <- unlist(strsplit(y, ","))

   # if (regexpr("dddd", annotfilter) > 0) {                    # get rid of high-name mirs.
   #    y2 <- grep("\\d\\d\\d\\d", y1, perl=TRUE, invert=TRUE, value=TRUE)
   # }
   # y1 <- gsub("-miR-", "-", y1)

         # species is either already specified, or we choose it to be the most common
         # occurring tag.
   if (nchar(species) == 0) {
         # get "bta", "mmu" et cetera.
      z <- unlist(lapply(y1, function(x) { strsplit(x, "-")[[1]][[1]] }))
         # count them.
      z2 <- aggregate(z, by=list(z), FUN=length)
         # and get the most frequent one.
      species <- z2$Group.1[order(z2$x)[nrow(z2)]]
   }

   s1 <- grep(sprintf("^%s\\b", species), y1, perl=TRUE)

   # hierverder. document logic. make sure names like hsap-miR-087111 go through.
   # I guess, even if the species name is mmu ...
   # Document Document Document Document Document Document
   # what about mmus or hsap (some funky long name)?

   if (length(s1) > 0) {
      x <- y1[s1]
      sel_family <- (regexpr("-\\d+[a-z]", x) > 0) & ! (regexpr("-\\d+-[35]p$", x) > 0)
      mir_family <- x[sel_family]
      mir_loners <- x[!sel_family]

      pat_mirs <- unique(gsub(".*?(-let-7|-\\d+)[a-z]?(-[35]p)?", "\\1([a-z])\\2$", mir_family, perl=TRUE, fixed=FALSE))
      pat_resubstitute <- gsub("(.*?let-7|\\d+)[a-z]?(-[35]p)?", "\\1#\\2", mir_family, perl=TRUE, fixed=FALSE)

      seen_mirs <- rep(FALSE, length(mir_family))

      display_result <- c(mir_loners)

      for (m in pat_mirs) {
         re_result <- regexpr(m, mir_family)
         if (sum((re_result > 0) & seen_mirs)) {
            print(paste("Error 1 in simplifying these names: ", mir_family))
         } else if (!sum(re_result)) {
            print(paste("Error 2 in simplifying these names: ", mir_family))
         }
         indices <- re_result > 0
         nhits <- length(unique(re_result[indices]))
         if (nhits > 1) {
            print(paste("Error 3 in simplifying these names: ", mir_family[indices]))
         } else if (nhits == 1) {
            family <- paste(sub(sprintf(".*?%s.*", m), "\\1", mir_family[indices]), collapse="")
            if (nchar(family) > 1) {
               family <- sprintf("[%s]", family)
            }
            display_result <- c(display_result, unique(sub("#", family, pat_resubstitute[indices])))
         }
      }
      display_result <- display_result[order(as.integer(gsub(".*?(\\d+).*", "\\1", display_result)))]
# print(display_result)
   } else {
      if (length(y1) > 0) {
         display_result <- c(y1[[1]])
      } else {
         display_result <- c()
      }
   }
   display_result <- gsub(sprintf("%s-", species), "", display_result)
   display_result <- gsub("miR-", "", display_result)
   display_result <- sub("^ ", "", display_result, perl=TRUE)
   paste(display_result, collapse=", ")
}


annot <- list()
t1 <- NULL

ouall <- c()         # bad organisation: will be used later.

# stopifnot(FALSE)
# what's x1 x2 .... why no pretty names?

if (nchar(fnannot) > 0) {
   t1 <- rt(fnannot, "r")
} else if (nchar(fntable) > 0) {
   x1 <- rt(fntable, "cr")
   rownames(x1) <- x1$eight1A
   if (hili_cons_mirna) {
      ouall <- rownames(x1[grep("NATIVE_CONS", x1$type),])
   }
   x2 <- x1[grep("ALIEN", x1$type, invert=TRUE),]
   x3 <- lapply(x2$lsomirna, normalise_mirna_names)
   names(x3) <- rownames(x2)
   t1 <- as.data.frame(cbind(V2=x3))

   stopifnot(batch)
}

if (!is.null(t1)) {
                        # fixme: 1:2 hardcoded to implement six3 filter.
   t6 <- unique(unlist(lapply(1:2, function(y) { unlist(lapply(rownames(t1), function(x) { allstrings(x, 6)[[y]] })) })))
   t7 <- unique(unlist(lapply(1:2, function(y) { unlist(lapply(rownames(t1), function(x) { allstrings(x, 7)[[y]] })) })))
   t8 <- unique(unlist(lapply(1:1, function(y) { unlist(lapply(rownames(t1), function(x) { allstrings(x, 8)[[y]] })) })))
   annot[t6] <- ""
   annot[t7] <- ""
   annot[t8] <- ""
                     # fixme: six3 filter hardcoded in allstrings and below 1:2 range.
                     # why does six use split and seven does not?
                     # hierverder.
                     # Also, try to sort seed/1 and seed/2.
   for (i in 1:2) {
      sixes <- unlist(lapply(rownames(t1), function(x) { allstrings(x, 6)[[i]] }))
      tmp <- t1[,"V2"]
      # tmp2 <- lapply(split(tmp, sixes), function(x) { paste(x, collapse=",") })
      tmp2 <- lapply(split(tmp, sixes), function(x) { sprintf("%s/%d", paste(x, collapse=","), 3-i) })
      annot[names(tmp2)] <- paste(annot[names(tmp2)], tmp2, sep = " - ")
   }
   for (i in 1:2) {
      # stopifnot(FALSE)
      sevenses <- unlist(lapply(rownames(t1), function(x) { allstrings(x, 7)[[i]] }))
      tmp <- t1[,"V2"]
      tmp2 <- lapply(split(tmp, sevenses), function(x) { sprintf("%s/%d", paste(x, collapse=","), 3-i) })
      annot[names(tmp2)] <- paste(annot[names(tmp2)], tmp2, sep = " - ")
   }
   eightses <- unlist(lapply(rownames(t1), function(x) { allstrings(x, 8)[[1]] }))
   # annot[eightses] <- paste(annot[eightses], paste(t1[,"V2"], rep("1", length(eightses)), sep="/"))
   annot[eightses] <- t1[,"V2"]
   annot <- lapply(annot, function(x) { sub(" - ", "", x) } )

   annot <- lapply(annot, function(x) { ifelse(nchar(x) > 45, sub(",[^,]*", "", substr(x, 1, 45)), x) })
   annot[nchar(annot) == 0] <- "na"
}



rich.colors <- function (n, palette = "temperature", rgb = FALSE, plot = FALSE)
{
    if (n <= 0)
        return(character(0))
    palette <- match.arg(palette, c("temperature", "blues"))
    x <- seq(0, 1, length = n)
    if (palette == "temperature") {
        r <- 1/(1 + exp(20 - 35 * x))
        g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
        b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
    }
    else {
        r <- 0.6 * x + 0.4 * x^2
        g <- 1.5 * x - 0.5 * x^2
        b <- 0.36 + 2.4 * x - 2 * x^2
        b[x > 0.4] <- 1
    }
    rgb.m <- matrix(c(r, g, b), ncol = 3, dimnames = list(as.character(seq(length = n)),
        c("red", "green", "blue")))
    rich.vector <- apply(rgb.m, 1, function(v) rgb(v[1], v[2],
        v[3]))
    if (rgb)
        attr(rich.vector, "rgb") <- rgb.m
    if (plot) {
        opar <- par("fig", "plt")
        par(fig = c(0, 1, 0, 0.7), plt = c(0.15, 0.9, 0.2, 0.95))
        plot(NA, xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01),
            xlab = "Spectrum", ylab = "", xaxs = "i", yaxs = "i",
            axes = FALSE)
        title(ylab = "Value", mgp = c(3.5, 0, 0))
        matlines(x, rgb.m, col = colnames(rgb.m), lty = 1, lwd = 3)
        matpoints(x, rgb.m, col = colnames(rgb.m), pch = 16)
        axis(1, at = 0:1)
        axis(2, at = 0:1, las = 1)
        par(fig = c(0, 1, 0.75, 0.9), plt = c(0.08, 0.97, 0,
            1), new = TRUE)
        midpoints <- barplot(rep(1, n), col = rich.vector, border = FALSE,
            space = FALSE, axes = FALSE)
        axis(1, at = midpoints, labels = 1:n, lty = 0, cex.axis = 0.6)
        par(opar)
    }
    return(rich.vector)
}


title <- gsub("%", " ", title)
sylamer <- "sylamer"

# Check to see if the file exists
if (!file.exists(fnsyl)) {
   stop(paste("Can't find the expected sylamer data file: ",fnsyl,sep=""), call.=FALSE)
}

# Read in and process the result
syltab <- read.table(fnsyl, sep="\t", row.names=1, header=T, check.names=F)
syltab <- cbind("0"=0,syltab)   # To add an initial column of 0s

pdf(file=pdfout, width=11, height=8)

pval_levels <- list()
p <- which(colnames(syltab) %in% c("pval_levels"))
if (length(p) > 0) {
   pval_levels <- syltab[,"pval_levels"]
   names(pval_levels) <- rownames(syltab)
   syltab <- syltab[,-c(p)]
}

xVals <- as.numeric(colnames(syltab))

u <- which(rownames(syltab) %in% c("utrlen"))
g <- which(rownames(syltab) %in% c("gccontent"))
n <- which(rownames(syltab) %in% c("nmasked"))

utrlen <- c()
if (length(u) > 0) {
   utrlen <- syltab["utrlen",]
}

gccontent <- c()
if (length(g) > 0) {
   gccontent <- syltab["gccontent",]
}

if (length(c(u,n,g)) > 0) {
   syltab <- syltab[-c(u,n,g),]
}


# Generate the plot area
if (ymin == 0) {
   ymin   <- min(syltab)
}
if (ymax == 0) {
   ymax   <- max(syltab)
}

yrange <- ymax - ymin
par(bg="white", plt=c(0.1,0.2, 0.9,0.8), mar=c(2,2,2,2), oma=c(2,2,2,2))
plot(NULL, xlab="Sorted sequences", ylab="log10(enrichment P-value)", axes=T, 
     main=title,
     ylim=c(round(ymin-yrange/10),round(ymax+yrange/10)), xlim=range(xVals))

# It can save time to plot no more than ~1,000 lines (particularly for all words of length 7 or 8)

if (n_max > 0) {
   syltab <- syltab[order(apply(syltab,1,function(x) max(abs(x))),decreasing=TRUE),]
   if (n_max > nrow(syltab)) {
      n_max <- nrow(syltab)
   }
} else {
   n_max <- nrow(syltab)
}

# Plot the background lines
for (i in 1:n_max) {
   lines(xVals,syltab[i,], col='grey')
}

# Draw a reference line at 0
abline(h=0)

# Up/Down best words
oligosUp <- c()
oligosDown <- c()

thek <- nchar(rownames(syltab)[1])

interval <- 1:ncol(syltab)
if (prefix_max > 0) {
   interval <- 1:prefix_max
}

if (topOligos > 0) { # Only if I really want these plots
   # oligosUp <- names((sort(apply(syltab,1, function(x) {max(x[is.finite(x)])}),decreasing=TRUE))[1:topOligos])
   if (topOligos > nrow(syltab)) {
      topOligos <- nrow(syltab)
   }
   oligosUp <- names((sort(apply(syltab,1, function(x) {max(x[interval])}),decreasing=TRUE))[1:topOligos])
}
if (lowOligos > 0) {
   # oligosDown <- names((sort(apply(syltab,1,function(x) {min(x[is.finite(x)])}),decreasing=FALSE))[1:lowOligos])
   oligosDown <- names((sort(apply(syltab,1,function(x) {min(x[interval])}),decreasing=FALSE))[1:lowOligos])
}
if (length(oligosChosen) > 0) {
   oligosChosen <- unique(unlist(lapply(oligosChosen, function(x) { allstrings(x, thek) })))
   oligosChosen <- oligosChosen[oligosChosen %in% rownames(syltab)]
   print(paste("have", length(oligosChosen), "oligos"))
   print(dim(syltab))
   oup <- oligosUp[!oligosUp %in% oligosChosen]
   oud <- oligosDown[!oligosDown %in% oligosChosen]
} else {
   oup <- oligosUp
   oud <- oligosDown
}
if (length(ouall) > 0) {
   paste("New code, take out this obstacle if you want to use it")
   stop()
   if (thek == 6) {
      ouall <- unique(unlist(lapply(ouall, function(x) { substring(x, 2, 7) })))
   } else if (thek == 7) {
      ouall <- unique(unlist(lapply(ouall, function(x) { substring(x, 1, 7) })))
   }
}

oligosAll <- c(oligosChosen,oup,oud, ouall)
                   # hierverder

print(ouall)



# oligosAll <- c(oligosChosen,oup,oud)
# 
# oligosAll <- list()
# length(oligosAll) <- 14
# 
# oligosAll[1:4] <- oligosChosen[1:4]
# oligosAll[11:14] <- oligosChosen[5:8]
# 
# ouall <- c(oup,oud)
# for (i in seq_along(ouall)) {
#    oligosAll[i+4] <- ouall[i]
# }



# if (length(oligosAll) != length(oligosChosen)) {
#    oligosAll <- names(sort(apply(syltab,1, function(x) {max(abs(x[is.finite(x)]))})[oligosAll],decreasing=TRUE))
# }


dyadise <- function(x) {
   if (nchar(dyad) > 0) {
      len <- nchar(x)
      halflen <- len / 2
      return(paste(substr(x, 1, halflen), dyad, substr(x, halflen+1, len), sep=""))
   }
   return(x)
}



if (length(oligosAll) > 0) {
   colors   <- rich.colors(length(oligosAll))
   names(colors) <- oligosAll
   par(family="mono")
   for (i in rev(seq_along(ouall))) {
      lines(xVals, syltab[ouall[i],], col=colors[ouall[i]], lwd=2)
   }
   for (i in rev(seq_along(oligosDown))) {
      lines(xVals, syltab[oligosDown[i],], col=colors[oligosDown[i]], lwd=2)
   }
   for (i in rev(seq_along(oligosUp))) {
      lines(xVals, syltab[oligosUp[i],], col=colors[oligosUp[i]], lwd=2)
   }
   for (i in rev(seq_along(oligosChosen))) {
      lines(xVals, syltab[oligosChosen[i],], col=colors[oligosChosen[i]], lwd=2)
   }

   if (bonferroni) {
      abline(h=log10(0.05 / nrow(syltab)), col="red")
      abline(h=-log10(0.05 / nrow(syltab)), col="red")
   }
   if (topOligos >0) {
      # lg <- oligosUp
      # len <- nchar(oligosUp[[1]])
      # halflen <- len / 2
      # if (length(dyad) > 0) {
      #    lg <- sapply(oligosUp, function(x) { paste(substr(x, 1, halflen), dyad, substr(x, halflen+1, len)) })
      # }
      # strings <- oligosUp
      # if (length(pval_levels) > 0) {
      #    strings <- paste(oligosUp, sprintf("%.2f", pval_levels[oligosUp]))
      # }
      if (length(annot) > 0) {
         a <- annot[oligosUp]
         a[is.na(names(a))] <- "na"
         thetext <- sprintf("%s %s", oligosUp, a)
      } else {
         thetext <- oligosUp
      }
      legend('topright', inset=c(0.0,0.0), legend=thetext, lwd=2, lty=1, horiz=FALSE, col=colors[unique(c(oligosUp,oligosChosen))],
             cex=0.8, bg='white')   # , title="Words with highest peak")
   }
   if (lowOligos >0) {
      if (length(annot) > 0) {
         a <- annot[oligosDown]
         a[is.na(names(a))] <- "na"
         thetext <- sprintf("%s %s", oligosDown, a)
      } else {
         thetext <- oligosDown
      }
      legend('bottomright', inset=c(0.0,0.0), legend=thetext, lwd=2, lty=1, horiz=FALSE, col=colors[oligosDown],
             cex=0.8, bg='white') # , title="Words with lowest peak")
   }
   if (length(oligosChosen) >0) {
      if (length(annot) > 0) {
         a <- annot[oligosChosen]
         a[is.na(names(a))] <- "na"
         thetext <- sprintf("%s %s", oligosChosen, a)
      } else {
         thetext <- oligosChosen
      }
      legend('bottomleft', inset=c(0.0,0.0), legend=thetext, lwd=2, lty=1, horiz=F, ncol=1, col=colors[oligosChosen],
             bg='white', title=selectedwords, cex=0.8)
   }
   if (nchar(secondlegend) > 0) {
      legend('topright', inset=c(0.0,0.0), legend=unlist(strsplit(secondlegend, " ## ")), ncol=1, cex=0.8, bg='white')
   }
# print(length(oligosChosen))
}

if (length(utrlen)) {
   points(xVals, utrlen[1,] / 100, col="black", pch=".")
   print("plotting utr lengths")
}

if (length(gccontent)) {
   points(xVals, gccontent[1,] * 100, col="red", pch=".")
   print("plotting gccontent")
}

dev.off()
