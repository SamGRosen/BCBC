install.packages("./cvxbiclustr", repos = NULL, type = "source")
install.packages("dbscan")

source("./RSCobra.R")

gammas <- 2 ^ seq(3, 20, 1)

run_cobra <- function(X, name) {
  print(paste(name, "dim", dim(X)))
  
  weight_start <- Sys.time()
  wts <- fast_gkn_weights(
    t(X),
    k_row = 4,
    k_col = 10,
    phi = 0.5,
    approx = 0
  )
  print(paste(
    "Weight time:",
    as.numeric(Sys.time() - weight_start, units = "secs"),
    "seconds"
  ))
  
  cv_start <- Sys.time()
  result <- cobra_validate(
    t(X),
    wts$E_row,
    wts$E_col,
    wts$w_row,
    wts$w_col,
    gamma = gammas,
    max_iter = 300
  )
  print(paste(
    "CV time:",
    as.numeric(Sys.time() - cv_start, units = "secs"),
    "seconds"
  ))
  print(warnings())
  
  saveRDS(result, paste0("./cobra_", name, ".RDS"))
}


lung500 = read.csv('data/lung500.csv', head = FALSE)
lung500 = scale(data.matrix(lung500))
run_cobra(lung500, "lung500")


lym <-
  t(read.csv("./data/lymphoma.x.txt", header = FALSE, sep = " "))
lym <- scale(lym)
run_cobra(lym, "lym")


srbct <-
  t(read.csv("./data/srbct.x.txt", header = FALSE, sep = " "))
srbct <- scale(srbct)
run_cobra(srbct, "srbct")


brain <-
  t(read.csv("./data/brain.x.txt", header = FALSE, sep = " "))
brain <- scale(brain)
run_cobra(brain, "brain")


colon <-
  t(read.csv("./data/colon.x.txt", header = FALSE, sep = " "))
colon <- scale(colon)
run_cobra(colon, "colon")


leuk <-
  t(read.csv("./data/leukemia.x.txt", header = FALSE, sep = " "))
leuk <- scale(leuk)
run_cobra(leuk, "leuk")


prostate <-
  t(read.csv("./data/prostate.x.txt", header = FALSE, sep = " "))
prostate <- scale(prostate)
run_cobra(prostate, "prostate")
