library(ape)
#opening tree
phy <- read.tree("tree.tree")
phy
length(phy$tip.label)
is.ultrametric(phy)
plot(phy, show.tip.label=F)
is.rooted(phy)
#Loading the data
mydata3 <- read.csv("mydata3.csv")
mydata3


#Setting the char pca1
charpca1.p <- mydata3[,2] 
names(charpca1.p) = as.character(mydata3[,1])
charpca1.p
#Setting the char bio18 variation
sdpca1 <- rep(0.5, 140)
names(sdpca1) = as.character(mydata3[,1])
sdpca1

library(diversitree)
p <- starting.point.quasse(phy, charpca1.p)
p
xr <- range(charpca1.p) + c(-1, 1) * 20 * p["diffusion"]
xr
linear.x <- make.linear.x(xr[1], xr[2])
make <- function(lambda, mu) make.quasse(phy, charpca1.p, sdpca1, lambda, mu, sampling.f=0.50)
nodrift <- function(f) constrain(f, drift ~ 0)
f.c.c <- make(constant.x, constant.x)
f.l.c <- make(linear.x, constant.x)
f.s.c <- make(sigmoid.x, constant.x) 
f.h.c <- make(noroptimal.x, constant.x) 
 
control <- list(parscale = 0.1, reltol = 0.001)
biopca1.mle.c.c <- find.mle(nodrift(f.c.c), p, lower = 0, control = control, verbose = 0)
biopca1.p.c <- biopca1.mle.c.c$par
biopca1.p.c
biopca1.p.l.c <- c(biopca1.p.c[1], l.m = 0, biopca1.p.c[2:3])
biopca1.p.l.c
biopca1.p.s.c <- c(biopca1.p.c[1], biopca1.p.c[1], mean(charpca1.p), 1, biopca1.p.c[2:3])
biopca1.p.s.c
biopca1.p.h.c <- c(biopca1.p.c[1], biopca1.p.c[1], mean(charpca1.p), 1, biopca1.p.c[2:3])
biopca1.p.h.c
biopca1.mle.d.c.c <- find.mle(f.c.c, coef(biopca1.mle.c.c, TRUE), control = control, verbose = 0)
biopca1.mle.l.c <- find.mle(nodrift(f.l.c), biopca1.p.l.c, control = control, verbose = 0)
biopca1.mle.d.l.c <- find.mle(f.l.c, coef(biopca1.mle.l.c, TRUE), control = control, verbose = 0)
biopca1.mle.s.c <- find.mle(nodrift(f.s.c), biopca1.p.s.c, control = control, verbose = 0)
biopca1.mle.d.s.c <- find.mle(f.s.c, coef(biopca1.mle.s.c, TRUE), control = control, verbose = 0)
biopca1.mle.h.c <- find.mle(nodrift(f.h.c), biopca1.p.h.c, control = control, verbose = 0)
biopca1.mle.d.h.c <- find.mle(f.h.c, coef(biopca1.mle.h.c, TRUE), control = control, verbose = 0)
biopca1.anova <- anova(biopca1.mle.c.c, OU.constant.constant = biopca1.mle.d.c.c, BM.linear.constant = biopca1.mle.l.c, OU.linear.constant = biopca1.mle.d.l.c, BM.sigmoid.constant = biopca1.mle.s.c, OU.sigmoid.constant = biopca1.mle.d.s.c, BM.hump.constant = biopca1.mle.h.c, OU.hump.constant = biopca1.mle.d.h.c)
biopca1.anova
save.image(file="biopca1myEnvironment.RData")
