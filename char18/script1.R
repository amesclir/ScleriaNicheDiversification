library(ape)
#opening tree
phy <- read.tree("tree.tree")
phy
length(phy$tip.label)
is.ultrametric(phy)
plot(phy, show.tip.label=F)
is.rooted(phy)
#Loading the data
mydata <- read.csv("mydata.csv")
mydata
mydata2 <- read.csv("mydata2.csv")
mydata2


setdiff(phy$tip.label, mydata[,1]) 
setdiff(mydata[,1],phy$tip.label) 
#The tree is ready for the analyses!!!!!!


#Setting the char bio18
char18.p <- mydata[,19] 
names(char18.p) = mydata[,1]
char18.p
#Setting the char bio18 variation
sd18 <- sqrt(mydata2[,19])
names(sd18) = mydata2[,1]
sd18
library(diversitree)
p <- starting.point.quasse(phy, char18.p)
p
xr <- range(char18.p) + c(-1, 1) * 20 * p["diffusion"]
xr
linear.x <- make.linear.x(xr[1], xr[2])
make <- function(lambda, mu) make.quasse(phy, char18.p, sd18, lambda, mu, sampling.f=0.50)
nodrift <- function(f) constrain(f, drift ~ 0)
f.c.c <- make(constant.x, constant.x)
f.l.c <- make(linear.x, constant.x)
f.s.c <- make(sigmoid.x, constant.x) 
f.h.c <- make(noroptimal.x, constant.x) 
 
control <- list(parscale = 0.1, reltol = 0.001)
bio18.mle.c.c <- find.mle(nodrift(f.c.c), p, lower = 0, control = control, verbose = 0)
bio18.p.c <- bio18.mle.c.c$par
bio18.p.c
bio18.p.l.c <- c(bio18.p.c[1], l.m = 0, bio18.p.c[2:3])
bio18.p.l.c
bio18.p.s.c <- c(bio18.p.c[1], bio18.p.c[1], mean(char18.p), 1, bio18.p.c[2:3])
bio18.p.s.c
bio18.p.h.c <- c(bio18.p.c[1], bio18.p.c[1], mean(char18.p), 1, bio18.p.c[2:3])
bio18.p.h.c
bio18.mle.d.c.c <- find.mle(f.c.c, coef(bio18.mle.c.c, TRUE), control = control, verbose = 0)
bio18.mle.l.c <- find.mle(nodrift(f.l.c), bio18.p.l.c, control = control, verbose = 0)
bio18.mle.d.l.c <- find.mle(f.l.c, coef(bio18.mle.l.c, TRUE), control = control, verbose = 0)
bio18.mle.s.c <- find.mle(nodrift(f.s.c), bio18.p.s.c, control = control, verbose = 0)
bio18.mle.d.s.c <- find.mle(f.s.c, coef(bio18.mle.s.c, TRUE), control = control, verbose = 0)
bio18.mle.h.c <- find.mle(nodrift(f.h.c), bio18.p.h.c, control = control, verbose = 0)
bio18.mle.d.h.c <- find.mle(f.h.c, coef(bio18.mle.h.c, TRUE), control = control, verbose = 0)
bio18.anova <- anova(bio18.mle.c.c, OU.constant.constant = bio18.mle.d.c.c, BM.linear.constant = bio18.mle.l.c, OU.linear.constant = bio18.mle.d.l.c, BM.sigmoid.constant = bio18.mle.s.c, OU.sigmoid.constant = bio18.mle.d.s.c, BM.hump.constant = bio18.mle.h.c, OU.hump.constant = bio18.mle.d.h.c)
bio18.anova
save.image(file="bio18myEnvironment.RData")
