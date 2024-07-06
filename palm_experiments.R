install.packages("./cvxbiclustr", repos = NULL, type = "source")
install.packages("dbscan")
install.packages("future.apply")

source("./RSCobra.R")

plan(sequential, split=TRUE)
# plan(multicore, workers = 4)


run_palm <- function(X, name, gammas = 2 ^ seq(10, 18, 1)) {
  print(paste(name, "dim", dim(X)))

  cv_start <- Sys.time()
  cv_results <- cv.BCBC(
    X,
    lambdas = 2 ^ seq(-5,-2, 0.25),
    nus = c(2),
    gammas = gammas,
    recalculate_weights = c(TRUE),
    k_cols = c(10),
    k_rows = c(4),
    tols=c(10^(-4.5)),
    percent_noise = 0.25,
    progress = TRUE
  )

  print(paste(
    "CV time:",
    as.numeric(Sys.time() - cv_start, units = "secs"),
    "seconds"
  ))
  print(warnings())

  saveRDS(cv_results, paste0("./palm_", name, ".RDS"))
}

#
# lung500 = read.csv('data/lung500.csv', head = FALSE)
# lung500 = scale(data.matrix(lung500))
# run_palm(lung500, "lung500")

lym <-
  t(read.csv("./data/lymphoma.x.txt", header = FALSE, sep = " "))
lym <- scale(lym)
run_palm(lym, "lym_tuned3", 2^seq(14, 17, 0.25))

#
# srbct <-
#   t(read.csv("./data/srbct.x.txt", header = FALSE, sep = " "))
# srbct <- scale(srbct)
# run_palm(srbct, "srbct_static", 2^seq(3, 13, 0.5))
#
#
# brain <-
#   t(read.csv("./data/brain.x.txt", header = FALSE, sep = " "))
# brain <- scale(brain)
# run_palm(brain, "brain_static", 2^seq(3, 13, 0.5))

#
# colon <-
#   t(read.csv("./data/colon.x.txt", header = FALSE, sep = " "))
# colon <- scale(colon)
# run_palm(colon, "colon")
#

leuk <-
  t(read.csv("./data/leukemia.x.txt", header = FALSE, sep = " "))
leuk <- scale(leuk)
run_palm(leuk, "leuk3", gammas = 2^seq(13, 16, 0.25))

#
# prostate <-
#   t(read.csv("./data/prostate.x.txt", header = FALSE, sep = " "))
# prostate <- scale(prostate)
# run_palm(prostate, "prostate")
