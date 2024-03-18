install.packages("./cvxbiclustr", repos = NULL, type = "source")
install.packages("dbscan")
install.packages("future.apply")

source("./RSCobra.R")

plan(sequential)
# plan(multisession, workers = 12)


run_rscobra <- function(X, name) {
  print(paste(name, "dim", dim(X)))
  
  cv_start <- Sys.time()
  cv_results <- tune_bcbc(
    X,
    model=rscobra,
    lambdas = 2 ^ seq(-14,-6, 1),
    nus = c(0.15),
    gammas = 2 ^ seq(3, 14, 1),
    recalculate_weights = c(TRUE),
    k_cols = c(10),
    k_rows = c(4),
    percent_noise = 0.25,
    tols=c(1e-5)
  )
  
  print(paste(
    "CV time:",
    as.numeric(Sys.time() - cv_start, units = "secs"),
    "seconds"
  ))
  print(warnings())
  
  saveRDS(cv_results, paste0("./rscobra_", name, ".RDS"))
}


lung500 = read.csv('data/lung500.csv', head = FALSE)
lung500 = scale(data.matrix(lung500))
run_rscobra(lung500, "lung500")


lym <-
  t(read.csv("./data/lymphoma.x.txt", header = FALSE, sep = " "))
lym <- scale(lym)
run_rscobra(lym, "lym")


srbct <-
  t(read.csv("./data/srbct.x.txt", header = FALSE, sep = " "))
srbct <- scale(srbct)
run_rscobra(srbct, "srbct")


brain <-
  t(read.csv("./data/brain.x.txt", header = FALSE, sep = " "))
brain <- scale(brain)
run_rscobra(brain, "brain")


colon <-
  t(read.csv("./data/colon.x.txt", header = FALSE, sep = " "))
colon <- scale(colon)
run_rscobra(colon, "colon")


leuk <-
  t(read.csv("./data/leukemia.x.txt", header = FALSE, sep = " "))
leuk <- scale(leuk)
run_rscobra(leuk, "leuk")


prostate <-
  t(read.csv("./data/prostate.x.txt", header = FALSE, sep = " "))
prostate <- scale(prostate)
run_rscobra(prostate, "prostate")
