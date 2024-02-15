print(R.home())
install.packages("./cvxbiclustr", repos=NULL, type="source")
install.packages("dbscan")

source("./RSCobra.R")


gammas <- 2 ^ seq(3, 20, 1)

print("lung500")

lung500 = read.csv('data/lung500.csv', head = FALSE)
lung500 = scale(data.matrix(lung500))

lung500_wts <- fast_gkn_weights(
  t(lung500),
  k_row = 4,
  k_col = 10,
  phi = 0.5,
  approx = 0
)

lung_result <- cobra_validate(t(lung500),
        lung500_wts$E_row,
        lung500_wts$E_col,
        lung500_wts$w_row,
        lung500_wts$w_col,
        gamma = gammas,
        max_iter = 300)

saveRDS(lung_result, "./cobra_lung500.RDS")
print(warnings())

print("lym")

lym <- t(read.csv("./data/lymphoma.x.txt", header=FALSE, sep=" "))
lym <- scale(lym)

lym_wts <- fast_gkn_weights(
  t(lym),
  k_row = 4,
  k_col = 10,
  phi = 0.5,
  approx = 0
)

lym_result <- cobra_validate(t(lym),
                             lym_wts$E_row,
                             lym_wts$E_col,
                             lym_wts$w_row,
                             lym_wts$w_col,
                             gamma = gammas,
                             max_iter = 300)

saveRDS(lym_result, "./cobra_lym.RDS")
print(warnings())

print("srbct")

srbct <- t(read.csv("./data/srbct.x.txt", header=FALSE, sep=" "))
srbct <- scale(srbct)

srbct_wts <- fast_gkn_weights(
  t(srbct),
  k_row = 4,
  k_col = 10,
  phi = 0.5,
  approx = 0
)

srbct_result <- cobra_validate(t(srbct),
                               srbct_wts$E_row,
                               srbct_wts$E_col,
                               srbct_wts$w_row,
                               srbct_wts$w_col,
                               gamma = gammas,
                               max_iter = 300)

saveRDS(srbct_result, "./cobra_srbct.RDS")
print(warnings())

print("brain")

brain <- t(read.csv("./data/brain.x.txt", header=FALSE, sep=" "))
brain <- scale(brain)

brain_wts <- fast_gkn_weights(
  t(brain),
  k_row = 4,
  k_col = 10,
  phi = 0.5,
  approx = 0
)

brain_result <- cobra_validate(t(brain),
                               brain_wts$E_row,
                               brain_wts$E_col,
                               brain_wts$w_row,
                               brain_wts$w_col,
                               gamma = gammas,
                               max_iter = 300)

saveRDS(brain_result, "./cobra_brain.RDS")
print(warnings())

print("colon")

colon <- t(read.csv("./data/colon.x.txt", header=FALSE, sep=" "))
colon <- scale(colon)

colon_wts <- fast_gkn_weights(
  t(colon),
  k_row = 4,
  k_col = 10,
  phi = 0.5,
  approx = 0
)

colon_result <- cobra_validate(t(colon),
                               colon_wts$E_row,
                               colon_wts$E_col,
                               colon_wts$w_row,
                               colon_wts$w_col,
                               gamma = gammas,
                               max_iter = 300)

saveRDS(colon_result, "./cobra_colon.RDS")
print(warnings())

print("leuk")
leuk <- t(read.csv("./data/leukemia.x.txt", header=FALSE, sep=" "))
leuk <- scale(leuk)


leuk_wts <- fast_gkn_weights(
  t(leuk),
  k_row = 4,
  k_col = 10,
  phi = 0.5,
  approx = 0
)

leuk_result <- cobra_validate(t(leuk),
                               leuk_wts$E_row,
                               leuk_wts$E_col,
                               leuk_wts$w_row,
                               leuk_wts$w_col,
                               gamma = gammas,
                               max_iter = 300)

saveRDS(leuk_result, "./cobra_leuk.RDS")
print(warnings())
