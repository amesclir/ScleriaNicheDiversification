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


#Setting the char pca3
charpca3.p <- mydata3[,4] 
names(charpca3.p) = as.character(mydata3[,1])
charpca3.p
#Setting the char bio18 variation
sdpca3 <- rep(0.5, 140)
names(sdpca3) = as.character(mydata3[,1])
sdpca3

library(diversitree)
p <- starting.point.quasse(phy, charpca3.p)
p
xr <- range(charpca3.p) + c(-1, 1) * 20 * p["diffusion"]
xr
linear.x <- make.linear.x(xr[1], xr[2])
make <- function(lambda, mu) make.quasse(phy, charpca3.p, sdpca3, lambda, mu, sampling.f=0.50)
nodrift <- function(f) constrain(f, drift ~ 0)
f.c.c <- make(constant.x, constant.x)
f.l.c <- make(linear.x, constant.x)
f.s.c <- make(sigmoid.x, constant.x) 
f.h.c <- make(noroptimal.x, constant.x) 
 
control <- list(parscale = 0.1, reltol = 0.001)
biopca3.mle.c.c <- find.mle(nodrift(f.c.c), p, lower = 0, control = control, verbose = 0)
biopca3.p.c <- biopca3.mle.c.c$par
biopca3.p.c
biopca3.p.l.c <- c(biopca3.p.c[1], l.m = 0, biopca3.p.c[2:3])
biopca3.p.l.c
biopca3.p.s.c <- c(biopca3.p.c[1], biopca3.p.c[1], mean(charpca3.p), 1, biopca3.p.c[2:3])
biopca3.p.s.c
biopca3.p.h.c <- c(biopca3.p.c[1], biopca3.p.c[1], mean(charpca3.p), 1, biopca3.p.c[2:3])
biopca3.p.h.c
biopca3.mle.d.c.c <- find.mle(f.c.c, coef(biopca3.mle.c.c, TRUE), control = control, verbose = 0)
biopca3.mle.l.c <- find.mle(nodrift(f.l.c), biopca3.p.l.c, control = control, verbose = 0)
biopca3.mle.d.l.c <- find.mle(f.l.c, coef(biopca3.mle.l.c, TRUE), control = control, verbose = 0)
biopca3.mle.s.c <- find.mle(nodrift(f.s.c), biopca3.p.s.c, control = control, verbose = 0)
biopca3.mle.d.s.c <- find.mle(f.s.c, coef(biopca3.mle.s.c, TRUE), control = control, verbose = 0)
biopca3.mle.h.c <- find.mle(nodrift(f.h.c), biopca3.p.h.c, control = control, verbose = 0)
biopca3.mle.d.h.c <- find.mle(f.h.c, coef(biopca3.mle.h.c, TRUE), control = control, verbose = 0)
biopca3.anova <- anova(biopca3.mle.c.c, OU.constant.constant = biopca3.mle.d.c.c, BM.linear.constant = biopca3.mle.l.c, OU.linear.constant = biopca3.mle.d.l.c, BM.sigmoid.constant = biopca3.mle.s.c, OU.sigmoid.constant = biopca3.mle.d.s.c, BM.hump.constant = biopca3.mle.h.c, OU.hump.constant = biopca3.mle.d.h.c)
biopca3.anova
save.image(file="biopca3myEnvironment.RData")
