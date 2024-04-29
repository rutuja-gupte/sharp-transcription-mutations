library(MASS)

lm1 <- lm(mpg ~ 1, mtcars)
lm2 <- lm(mpg ~ (.)^2, mtcars)
lm.forward = stepAIC(lm1, direction="forward", scope=list(upper=lm2,lower=lm1))
summary(lm.forward)

lm.backward = stepAIC(lm1, direction="backward", scope=list(upper=lm2,lower=lm1))
summary(lm.backward)

lm.both = stepAIC(lm1, direction="both", scope=list(upper=lm2,lower=lm1))
summary(lm.both)

rm(mtcars)
rm(lm.backward)
rm(lm.forward)
rm(lm.both)
rm(lm1)
rm(lm2)

data <- read.csv("processed_data/weighted.csv")
data <- data[, c("SNM", "ref.type", "CDS", "collision", "rep.ratio", "net_tot")]

glm1 <- glm(SNM ~ 1, data, family="binomial")
glm2 <- glm(SNM ~ (.)^2, data, family="binomial")

lm.both = stepAIC(lm1, direction="both", scope=list(upper=lm2,lower=lm1))
summary(lm.both)