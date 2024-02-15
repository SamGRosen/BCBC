install.packages("sparseBC")
library(sparseBC)


gammas <- 2 ^ seq(3, 20, 1)


run_sparsebc <- function(X, name) {
  print(paste(name, "dim", dim(X)))
  
  tuning_start <- Sys.time()
  kr <-
    sparseBC.choosekr(
      X,
      k = 1:8,
      r = 1:25,
      lambda = 0,
      percent = 0.125,
      TRUE
    )
  lambda <- sparseBC.BIC(X,
                         kr$bestK,
                         kr$bestR,
                         gammas)
  print(paste(
    "Tuning time:",
    as.numeric(Sys.time() - tuning_start, units = "secs")
  ))
  
  fit_start <- Sys.time()
  sparsebc_fit <- sparseBC(lung500, kr$bestK,
                           kr$bestR, lambda$lambda)
  print(paste("Fit time:", as.numeric(Sys.time() - fit_start, units = "secs")))
  print(warnings())
  
  saveRDS(sparsebc_fit, paste0("./sparsebc_", name, ".RDS"))
}


lung500 = read.csv('data/lung500.csv', head = FALSE)
lung500 = scale(data.matrix(lung500))
run_sparsebc(lung500, "lung500")


lym <-
  t(read.csv("./data/lymphoma.x.txt", header = FALSE, sep = " "))
lym <- scale(lym)
run_sparsebc(lym, "lym")


srbct <-
  t(read.csv("./data/srbct.x.txt", header = FALSE, sep = " "))
srbct <- scale(srbct)
run_sparsebc(srbct, "srbct")


brain <-
  t(read.csv("./data/brain.x.txt", header = FALSE, sep = " "))
brain <- scale(brain)
run_sparsebc(brain, "brain")


colon <-
  t(read.csv("./data/colon.x.txt", header = FALSE, sep = " "))
colon <- scale(colon)
run_sparsebc(colon, "colon")


leuk <-
  t(read.csv("./data/leukemia.x.txt", header = FALSE, sep = " "))
leuk <- scale(leuk)
run_sparsebc(leuk, "leuk")


prostate <-
  t(read.csv("./data/prostate.x.txt", header = FALSE, sep = " "))
prostate <- scale(prostate)
run_sparsebc(prostate, "prostate")