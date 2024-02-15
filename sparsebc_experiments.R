print(R.home())
install.packages("sparseBC")
library(sparseBC)


gammas <- 2 ^ seq(3, 20, 1)

print("lung500")

lung500 = read.csv('data/lung500.csv', head = FALSE)
lung500 = scale(data.matrix(lung500))

lung_kr <- sparseBC.choosekr(lung500, k=1:8, r=1:25, lambda=0, percent=0.125, TRUE)
lung_lambda <- sparseBC.BIC(lung500, lung_kr$estimated_kr[1],
                            lung_kr$estimated_kr[2], gammas)

lung_fit <- sparseBC(lung500, lung_kr$estimated_kr[1],
                     lung_kr$estimated_kr[2], lung_lambda$lambda)

saveRDS(lung_fit, "sparsebc_lung500.RDS")

print("lym")
lym <- t(read.csv("./data/lymphoma.x.txt", header=FALSE, sep=" "))
lym <- scale(lym)

lym_kr <- sparseBC.choosekr(lym, k=1:8, r=1:25, lambda=0, percent=0.125, TRUE)
lym_lambda <- sparseBC.BIC(lym, lym_kr$estimated_kr[1],
                           lym_kr$estimated_kr[2], gammas)

lym_fit <- sparseBC(lym, lym_kr$estimated_kr[1],
                     lym_kr$estimated_kr[2], lym_lambda$lambda)

saveRDS(lym_fit, "sparsebc_lym.RDS")


print("srcbt")
srbct <- t(read.csv("./data/srbct.x.txt", header=FALSE, sep=" "))
srbct <- scale(srbct)

srbct_kr <- sparseBC.choosekr(srbct, k=1:8, r=1:25, lambda=0, percent=0.125, TRUE)
srbct_lambda <- sparseBC.BIC(srbct, srbct_kr$estimated_kr[1],
                             srbct_kr$estimated_kr[2], gammas)

srbct_fit <- sparseBC(srbct, srbct_kr$estimated_kr[1],
                    srbct_kr$estimated_kr[2], srbct_lambda$lambda)

saveRDS(srbct_fit, "sparsebc_srbct.RDS")


print("brain")
brain <- t(read.csv("./data/brain.x.txt", header=FALSE, sep=" "))
brain <- scale(brain)
brain_kr <- sparseBC.choosekr(brain, k=1:8, r=1:25, lambda=0, percent=0.125, TRUE)
brain_lambda <- sparseBC.BIC(brain, brain_kr$estimated_kr[1],
                             brain_kr$estimated_kr[2], gammas)

brain_fit <- sparseBC(brain, brain_kr$estimated_kr[1],
                      brain_kr$estimated_kr[2], brain_lambda$lambda)

saveRDS(brain_fit, "sparsebc_brain.RDS")


print("colon")
colon <- t(read.csv("./data/colon.x.txt", header=FALSE, sep=" "))
colon <- scale(colon)
colon_kr <- sparseBC.choosekr(colon, k=1:8, r=1:25, lambda=0, percent=0.125, TRUE)
colon_lambda <- sparseBC.BIC(colon, colon_kr$estimated_kr[1],
                             colon_kr$estimated_kr[2], gammas)

colon_fit <- sparseBC(colon, colon_kr$estimated_kr[1],
                      colon_kr$estimated_kr[2], colon_lambda$lambda)

saveRDS(colon_fit, "sparsebc_colon.RDS")

print("leuk")
leuk <- t(read.csv("./data/leukemia.x.txt", header=FALSE, sep=" "))
leuk <- scale(leuk)
leuk_kr <- sparseBC.choosekr(leuk, k=1:8, r=1:25, lambda=0, percent=0.125, TRUE)
leuk_lambda <- sparseBC.BIC(leuk, leuk_kr$estimated_kr[1],
                             leuk_kr$estimated_kr[2], gammas)

leuk_fit <- sparseBC(leuk, leuk_kr$estimated_kr[1],
                      leuk_kr$estimated_kr[2], leuk_lambda$lambda)

saveRDS(leuk_fit, "sparsebc_leuk.RDS")
