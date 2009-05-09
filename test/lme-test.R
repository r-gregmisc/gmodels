library(gmodels)
library(lme4)

sleepstudy$dayGroup <- cut(sleepstudy$Days, seq(-1,9,by=2), include=T)

# ci example
fm2 <- lmer(Reaction ~ dayGroup + (1|Subject) + (0+Days|Subject), sleepstudy)
ci(fm2)


# estimable examples
estimable(fm2, c( 0, -1, 1, 0,  0 )  ) # list all terms
estimable(fm2, c("dayGroup(1,3]"=-1, "dayGroup(3,5]"=1)) # just the nonzero terms
estimable(fm2, c("dayGroup(1,3]"=-1, "dayGroup(3,5]"=1), n.sim=5000 ) # more simulations...


# fit.contrast example
fit.contrast( fm2, "dayGroup",
              rbind("0-1 vs 3-4"=c(-1,0,1,0,0),
                    "3-4 vs 5-6"=c(0,0,-1,1,0)
                  ),
            conf=0.95 )

