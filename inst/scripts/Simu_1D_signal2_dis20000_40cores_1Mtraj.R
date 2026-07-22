
## general parameters

dimension = 1
sim_path="Simulations"
maxbatchsize = 10**7
signal_dim = 2

## parameters for trajectories simulation

Ndiscretisation = 20000
deltagrid = seq(200,1,-1)/200 # en ordre décroissant, valeur 1 permise

## Parameters for simulation quantiles estimation

Ntraj_batch = 500
alphagrid = seq(500,1,-1)/1000

## Parameters of Monte-carlo estimations

Ncores = 32
Ntrajmin = 5e5

## File names parameters

Namedeltagrid = "1to.005by.005"

# library loading

library(posir)
library(future)
logger::log_threshold(logger::DEBUG, namespace = "posir")

# Exec

plan(multisession, workers = Ncores)
compute_quantiles(Ntrajmin, Ntraj_batch, Ndiscretisation, deltagrid,
                  Namedeltagrid, alphagrid, maxbatchsize, sim_path,
                  d=dimension, signal_dim = signal_dim)
plan(sequential)
