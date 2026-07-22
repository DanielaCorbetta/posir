## general parameters

dimension <- 1
signal_dim <- 2
sim_path <- "/scratch/dcorbett/posird/Simulations"
maxbatchsize <- 5e7

## parameters for trajectories simulation

Ndiscretisation <- 20000
deltagrid <- seq(200, 1, -1) / 200

## Parameters for simulation quantiles estimation

Ntraj_batch <- 500
alphagrid <- seq(500, 1, -1) / 1000

## Parameters of Monte-Carlo estimations

Ncores <- 32
Ntrajmin <- 5e5

## File names parameters

Namedeltagrid <- "1to.005by.005"

## library loading

library(posir)
library(future)
library(doFuture)

logger::log_threshold(logger::DEBUG, namespace = "posir")

## Exec

plan(multisession, workers = Ncores)

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
