dimension <- 1
signal_dim <- 2
sim_path <- "Simulations_test"
maxbatchsize <- 1e6

# Much coarser discretisation and delta grid
Ndiscretisation <- 500
deltagrid <- seq(50, 1, -1) / 50   # delta = 1, 0.98, ..., 0.02

# Fewer trajectories, but batches large enough to test the whole machinery
Ntraj_batch <- 100
alphagrid <- c(0.2, 0.1, 0.05, 0.01)

Ncores <- 4
Ntrajmin <- 10^4

Namedeltagrid <- "1to.02by.02"

library(posir)
library(future)
library(doFuture)
logger::log_threshold(logger::DEBUG, namespace = "posir")

# plan(multisession, workers = Ncores)


# future::plan(future::sequential)
future::plan(multisession, workers = 4)
unlink("Simulations_test", recursive = TRUE)
Q <- compute_quantiles(
  Ntot = Ntrajmin,
  Batchsize = Ntraj_batch,
  Ndis = Ndiscretisation,
  deltagrid = deltagrid,
  gridname = Namedeltagrid,
  alphagrid = alphagrid,
  maxmatrixsize = maxbatchsize,
  sim_root_dir = sim_path,
  d = dimension,
  signal_dim = signal_dim
)

plan(sequential)

print(Q)
