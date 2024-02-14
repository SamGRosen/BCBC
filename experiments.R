print(R.home())
install.packages("./cvxbiclustr", repos=NULL, type="source")
install.packages("dbscan")
install.packages("future.apply")

source("./RSCobra.R")

plan(sequential)
# plan(multisession, workers = 12)

print("lung500")

lung500 = read.csv('data/lung500.csv', head = FALSE)
lung500 = scale(data.matrix(lung500))

cv_500 <- tune_rscobra(lung500,
                       lambdas=2 ^ seq(-14, -6, 1),
                       nus=c(0.15),
                       gammas=2 ^ seq(3, 14, 1),
                       recalculate_weights = c(TRUE),
                       k_cols = c(10),
                       k_rows = c(4),
                       percent_noise = 0.25)

saveRDS(cv_500, "./cv_lung500.RDS")
print(warnings())

print("lym")

lym <- t(read.csv("./data/lymphoma.x.txt", header=FALSE, sep=" "))
lym <- scale(lym)

cv_lym <- tune_rscobra(lym,
                       lambdas=2 ^ seq(-14, -6, 1),
                       nus=c(0.15),
                       gammas=2 ^ seq(3, 14, 1),
                       recalculate_weights = c(TRUE),
                       k_cols = c(10),
                       k_rows = c(4),
                       percent_noise = 0.25)

saveRDS(cv_lym, "./cv_lym.RDS")
print(warnings())

print("srbct")

srbct <- t(read.csv("./data/srbct.x.txt", header=FALSE, sep=" "))
srbct <- scale(srbct)

cv_srbct <- tune_rscobra(srbct,
                         lambdas=2 ^ seq(-14, -6, 1),
                         nus=c(0.15),
                         gammas=2 ^ seq(3, 14, 1),
                         recalculate_weights = c(TRUE),
                         k_cols = c(10),
                         k_rows = c(4),
                         percent_noise = 0.25)

saveRDS(cv_srbct, "./cv_srbct.RDS")
print(warnings())

print("brain")

brain <- t(read.csv("./data/brain.x.txt", header=FALSE, sep=" "))
brain <- scale(brain)

cv_brain <- tune_rscobra(brain,
                         lambdas=2 ^ seq(-14, -6, 1),
                         nus=c(0.15),
                         gammas=2 ^ seq(3, 14, 1),
                         recalculate_weights = c(TRUE),
                         k_cols = c(10),
                         k_rows = c(4),
                         percent_noise = 0.25)

saveRDS(cv_brain, "./cv_brain.RDS")
print(warnings())

print("colon")
colon <- t(read.csv("./data/colon.x.txt", header=FALSE, sep=" "))
colon <- scale(colon)

cv_colon <- tune_rscobra(colon,
                         lambdas=2 ^ seq(-14, -6, 1),
                         nus=c(0.15),
                         gammas=2 ^ seq(3, 14, 1),
                         recalculate_weights = c(TRUE),
                         k_cols = c(10),
                         k_rows = c(4),
                         percent_noise = 0.25)

saveRDS(cv_colon, "./cv_colon.RDS")
print(warnings())

print("leuk")
leuk <- t(read.csv("./data/leukemia.x.txt", header=FALSE, sep=" "))
leuk <- scale(leuk)

cv_leuk <- tune_rscobra(leuk,
                        lambdas=2 ^ seq(-14, -6, 1),
                        nus=c(0.15),
                        gammas=2 ^ seq(3, 14, 1),
                        recalculate_weights = c(TRUE),
                        k_cols = c(10),
                        k_rows = c(4),
                        percent_noise = 0.25)

saveRDS(cv_leuk, "./cv_leuk.RDS")
print(warnings())

