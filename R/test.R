source("R/packages.R")
source("R/functions.R")


china_data <- readd("data")

hubei <- china_data %>%
  filter(Province == 'Hubei')

data <- hubei$data[[1]]
N <- hubei$N

R0 <- 6.2
beta <- 1/8
I0_factor <- 10
new_times = 1:70

fun <- 'SIRX'

method <- 'quarantine'
parms <- switch(method,
  'confirmed'  = c(beta = beta, alpha = R0 * beta, kappa = beta, kappa0 = beta / 2, I0_factor = I0_factor, N = 1e7),
  'shutdown'   = c(beta = beta, alpha = R0 * beta, kappa =    0, kappa0 =     beta, I0_factor = I0_factor, N = 1e7),
  'quarantine' = c(beta = beta, alpha = R0 * beta, kappa = beta, kappa0 =        0, I0_factor = I0_factor, N = 1000)
)

lower <- switch(method,
  'confirmed'  = c(beta = beta, alpha = R0 * beta, kappa = 0, kappa0 = 0, I0_factor = 0.001, N = 1000),
  'shutdown'   = c(beta = beta, alpha = R0 * beta, kappa = 0, kappa0 = 0, I0_factor = 1    , N = 100000),
  'quarantine' = c(beta = beta, alpha = R0 * beta, kappa = 0.000001, kappa0 = 0, I0_factor = 0.001, N = 10)
)

upper <-  switch(method,
  'confirmed'  = c(beta = beta, alpha = R0 * beta, kappa = Inf, kappa0 = Inf, I0_factor = Inf, N = {if(exists("N")) N else 115000000}),
  'shutdown'   = c(beta = beta, alpha = R0 * beta, kappa = 0,   kappa0 = Inf, I0_factor = Inf, N = {if(exists("N")) N else 115000000}),
  'quarantine' = c(beta = beta, alpha = R0 * beta, kappa = Inf, kappa0 = 0,   I0_factor = Inf, N = {if(exists("N")) N else 1000000000})
)

opt <- nls.lm(
  par = parms, 
  fn = RSS,
  # method = "L-BFGS-B",
  lower = lower,
  upper = upper,
  data = data,
  fun = fun,
  N = N,
  control = nls.lm.control(
    ftol = .Machine$double.eps^2
    # ptol = .Machine$double.eps
  )
) 

init <- data %>%
  slice(1) %>%
  transmute(
    X = Infected / N,
    I = X * I0_factor,
    S = 1 - I - X,
    R = 0
  ) %>% 
  unlist() %>%
  set_names(c("X", "I", "S", "R"))


new <- ode(y = init, times = new_times, func = eval(as.symbol(fun)), parms = opt$par, method = 'ode45') %>%
  as.data.frame() %>%
  map_dfc(as.numeric) %>%
  mutate_at(vars(-time), ~ .x * opt$par["N"])
