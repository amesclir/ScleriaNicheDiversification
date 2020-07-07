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


#Setting the char bio15
char15.p <- mydata[,16] 
names(char15.p) = mydata[,1]
char15.p
#Setting the char bio15 variation
sd15 <- sqrt(mydata2[,16])
names(sd15) = mydata2[,1]
sd15
library(diversitree)
p <- starting.point.quasse(phy, char15.p)
p
xr <- range(char15.p) + c(-1, 1) * 20 * p["diffusion"]
xr
linear.x <- make.linear.x(xr[1], xr[2])
make <- function(lambda, mu) make.quasse(phy, char15.p, sd15, lambda, mu, sampling.f=0.50)
nodrift <- function(f) constrain(f, drift ~ 0)
f.c.c <- make(constant.x, constant.x)
f.l.c <- make(linear.x, constant.x)
f.s.c <- make(sigmoid.x, constant.x) 
f.h.c <- make(noroptimal.x, constant.x) 
 
control <- list(parscale = 0.1, reltol = 0.001)
bio15.mle.c.c <- find.mle(nodrift(f.c.c), p, lower = 0, control = control, verbose = 0)
bio15.p.c <- bio15.mle.c.c$par
bio15.p.c
bio15.p.l.c <- c(bio15.p.c[1], l.m = 0, bio15.p.c[2:3])
bio15.p.l.c
bio15.p.s.c <- c(bio15.p.c[1], bio15.p.c[1], mean(char15.p), 1, bio15.p.c[2:3])
bio15.p.s.c
bio15.p.h.c <- c(bio15.p.c[1], bio15.p.c[1], mean(char15.p), 1, bio15.p.c[2:3])
bio15.p.h.c
bio15.mle.d.c.c <- find.mle(f.c.c, coef(bio15.mle.c.c, TRUE), control = control, verbose = 0)
bio15.mle.l.c <- find.mle(nodrift(f.l.c), bio15.p.l.c, control = control, verbose = 0)
bio15.mle.d.l.c <- find.mle(f.l.c, coef(bio15.mle.l.c, TRUE), control = control, verbose = 0)
bio15.mle.s.c <- find.mle(nodrift(f.s.c), bio15.p.s.c, control = control, verbose = 0)
bio15.mle.d.s.c <- find.mle(f.s.c, coef(bio15.mle.s.c, TRUE), control = control, verbose = 0)
bio15.mle.h.c <- find.mle(nodrift(f.h.c), bio15.p.h.c, control = control, verbose = 0)
bio15.mle.d.h.c <- find.mle(f.h.c, coef(bio15.mle.h.c, TRUE), control = control, verbose = 0)
bio15.anova <- anova(bio15.mle.c.c, OU.constant.constant = bio15.mle.d.c.c, BM.linear.constant = bio15.mle.l.c, OU.linear.constant = bio15.mle.d.l.c, BM.sigmoid.constant = bio15.mle.s.c, OU.sigmoid.constant = bio15.mle.d.s.c, BM.hump.constant = bio15.mle.h.c, OU.hump.constant = bio15.mle.d.h.c)
bio15.anova
save.image(file="bio15myEnvironment.RData")
