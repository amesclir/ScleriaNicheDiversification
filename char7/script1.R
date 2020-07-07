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

#Setting the char bio7
char7.p <- mydata[,8] 
names(char7.p) = mydata[,1]
char7.p
#Setting the char bio7 variation
sd7 <- sqrt(mydata2[,8])
names(sd7) = mydata2[,1]
sd7
library(diversitree)
p <- starting.point.quasse(phy, char7.p)
p
xr <- range(char7.p) + c(-1, 1) * 20 * p["diffusion"]
xr
linear.x <- make.linear.x(xr[1], xr[2])
make <- function(lambda, mu) make.quasse(phy, char7.p, sd7, lambda, mu, sampling.f=0.50)
nodrift <- function(f) constrain(f, drift ~ 0)
f.c.c <- make(constant.x, constant.x)
f.l.c <- make(linear.x, constant.x)
f.s.c <- make(sigmoid.x, constant.x) 
f.h.c <- make(noroptimal.x, constant.x) 
 
control <- list(parscale = 0.1, reltol = 0.001)
bio7.mle.c.c <- find.mle(nodrift(f.c.c), p, lower = 0, control = control, verbose = 0)
bio7.p.c <- bio7.mle.c.c$par
bio7.p.c
bio7.p.l.c <- c(bio7.p.c[1], l.m = 0, bio7.p.c[2:3])
bio7.p.l.c
bio7.p.s.c <- c(bio7.p.c[1], bio7.p.c[1], mean(char7.p), 1, bio7.p.c[2:3])
bio7.p.s.c
bio7.p.h.c <- c(bio7.p.c[1], bio7.p.c[1], mean(char7.p), 1, bio7.p.c[2:3])
bio7.p.h.c
bio7.mle.d.c.c <- find.mle(f.c.c, coef(bio7.mle.c.c, TRUE), control = control, verbose = 0)
bio7.mle.l.c <- find.mle(nodrift(f.l.c), bio7.p.l.c, control = control, verbose = 0)
bio7.mle.d.l.c <- find.mle(f.l.c, coef(bio7.mle.l.c, TRUE), control = control, verbose = 0)
bio7.mle.s.c <- find.mle(nodrift(f.s.c), bio7.p.s.c, control = control, verbose = 0)
bio7.mle.d.s.c <- find.mle(f.s.c, coef(bio7.mle.s.c, TRUE), control = control, verbose = 0)
bio7.mle.h.c <- find.mle(nodrift(f.h.c), bio7.p.h.c, control = control, verbose = 0)
bio7.mle.d.h.c <- find.mle(f.h.c, coef(bio7.mle.h.c, TRUE), control = control, verbose = 0)
bio7.anova <- anova(bio7.mle.c.c, OU.constant.constant = bio7.mle.d.c.c, BM.linear.constant = bio7.mle.l.c, OU.linear.constant = bio7.mle.d.l.c, BM.sigmoid.constant = bio7.mle.s.c, OU.sigmoid.constant = bio7.mle.d.s.c, BM.hump.constant = bio7.mle.h.c, OU.hump.constant = bio7.mle.d.h.c)
bio7.anova
save.image(file="bio7myEnvironment.RData")
