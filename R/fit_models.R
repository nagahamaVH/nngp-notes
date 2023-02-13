library(rstan)
# library(cmdstanr)
source("R/NNMatrix.R")

# ------------------- Setup ---------------------------------------------------
options(mc.cores = parallel::detectCores())

data_board <- pins::board_folder("./data", versioned = T)
model_board <- pins::board_folder("models", versioned = T)

seed <- 852
models <- c("poisson_nngp.stan", "poisson_gp.stan")
neighbor_size <- c(5, 10, 30)

# ------------------------ Stan parameters ------------------------------------
n_chain <- 4
n_iter <- 1000
# -----------------------------------------------------------------------------

datasets <- pins::pin_list(data_board)

pb <- progress::progress_bar$new(
  total = length(datasets) * (1 + length(neighbor_size))
)

pb$tick(0)
for (j in 1:length(datasets)) {
  data <- data_board |>
    pins::pin_read(datasets[j])
  
  data_version <- data_board |>
    pins::pin_versions(datasets[j]) |>
    tail(1) |>
    dplyr::pull(version)
  
  for (i in 1:length(models)){
    model_file <- file.path("stan", models[i])
    # mod_cl <- cmdstan_model(model_file, cpp_options = list(stan_opencl = T))
    mod_cl <- stan_model(model_file)
    
    # Adding extra data for stan data block
    if (models[i] == "poisson_nngp.stan") {
      for (k in 1:length(neighbor_size)) {
        pb$tick()
        m <- neighbor_size[k]
        nn <- NNMatrix(data$coords, m)
        
        data$NN_ind <- nn$NN_ind
        data$NN_dist <- nn$NN_dist
        data$NN_distM <- nn$NN_distM
        data$M <- m
        
        t_init <- proc.time()
        # fit_cl <- mod_cl$sample(
        #   data = data, 
        #   chains = n_chain,
        #   iter_warmup = floor(n_iter / 2),
        #   iter_sampling = n_iter,
        #   show_messages = F,
        #   seed = seed)
        fit_cl <- sampling(
          mod_cl,
          data = data,
          chains = n_chain,
          iter = n_iter,
          seed = seed,
          refresh = 0)
        t_total <- proc.time() - t_init
        
        model_meta <- list(
          file = model_file,
          n_neighbors = m,
          data = data_version,
          n_it = n_iter,
          n_chain = n_chain,
          time = unname(t_total[3]),
          # model_code = mod_cl$code() |>
          #   paste(collapse = "\n")
          model_code = paste(mod_cl@model_code, collapse = "\n")
        )
        
        tryCatch({
          pins::pin_write(
            model_board, fit_cl, 
            name = paste0(models[i], "_n", data$N, "_m", m), 
            type = "rds", metadata = model_meta)
          rm(fit_cl)
        }, error = function(e){
          cat("ERROR:", conditionMessage(e), "\n")
        })
      }
    } else{
      pb$tick()
      
      m <- 0

      t_init <- proc.time()
      # fit_cl <- mod_cl$sample(
      #   data = data, 
      #   chains = n_chain,
      #   iter_warmup = floor(n_iter / 2),
      #   iter_sampling = n_iter,
      #   seed = seed,
      #   refresh = 0)
      fit_cl <- sampling(
        mod_cl,
        data = data,
        chains = n_chain,
        iter = n_iter,
        seed = seed,
        refresh = 0)
      t_total <- proc.time() - t_init
      
      model_meta <- list(
        file = model_file,
        n_neighbors = m,
        data = data_version,
        n_it = n_iter,
        n_chain = n_chain,
        time = t_total[3],
        # model_code = mod_cl$code() |>
        #   paste(collapse = "\n")
        model_code = paste(mod_cl@model_code, collapse = "\n")
      )
      
      tryCatch({
        pins::pin_write(
          model_board, fit_cl, 
          name = paste0(models[i], "_n", data$N, "_m", m), 
          type = "rds", metadata = model_meta)
        rm(fit_cl)
      }, error = function(e){
        cat("ERROR:", conditionMessage(e), "\n")
      })
    }
  }
}
