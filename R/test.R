source("R/packages.R")   # load packages (packages in R are similar to modules in python, but without as many namespace issues)
source("R/functions.R")  # this is where I've defined the functions - the relevant functions are at the top of the file.



# get the data
china_data <- readd("data")
print(china_data)

hubei <- china_data %>%
  filter(Province == 'Hubei')

# total infections/date data
data <- hubei$data[[1]]
print(data)

N <- hubei$N

# set parameters for the model
R0 <- 6

beta <- 1/8
I0_factor <- 10
new_times = 1:70

fun <- 'SIRX'

method <- 'shutdown'
parms <- switch(method,
                'confirmed'  = c(beta = beta, alpha = R0 * beta, kappa = beta, kappa0 = beta / 2, I0_factor = I0_factor, N = 1e7),
                'shutdown'   = c(beta = beta, alpha = R0 * beta, kappa =    0, kappa0 =     beta, I0_factor = I0_factor, N = 1e7),
                'quarantine' = c(beta = beta, alpha = R0 * beta, kappa = beta, kappa0 =        0, I0_factor = I0_factor, N = 10000)
)

lower <- switch(method,
                'confirmed'  = c(beta = beta, alpha = R0 * beta, kappa = 0, kappa0 = 0, I0_factor = 0.001, N = 1000),
                'shutdown'   = c(beta = beta, alpha = R0 * beta, kappa = 0, kappa0 = 0, I0_factor = 1    , N = 100000),
                'quarantine' = c(beta = beta, alpha = R0 * beta, kappa = 0, kappa0 = 0, I0_factor = 0.001, N = 1000)
)

upper <-  switch(method,
                 'confirmed'  = c(beta = beta, alpha = R0 * beta, kappa = Inf, kappa0 = Inf, I0_factor = Inf, N = {if(exists("N")) N else 115000000}),
                 'shutdown'   = c(beta = beta, alpha = R0 * beta, kappa = 0,   kappa0 = Inf, I0_factor = Inf, N = {if(exists("N")) N else 115000000}),
                 'quarantine' = c(beta = beta, alpha = R0 * beta, kappa = Inf, kappa0 = 0,   I0_factor = Inf, N = {if(exists("N")) N else 1000000000})
)

opt <- nls.lm(
  par = parms, 
  fn = RSS,
  lower = lower,
  upper = upper,
  control = nls.lm.control(maxiter = 1000, maxfev = 5000),
  data = data,
  fun = fun,
  N = N)  %>%
  print()

# the print methods automatically gives pretty output when printing.  It automatically kicks in when you just type "opt", so to dig into the contents, you can either use "str(opt)"
# or just look at the elements

# access just the parameters
opt$par
# If you type "opt$", Rstudio brings up a list of all the elements within opt (such as par).  
# ?nls.lm or help("nls.lm") brings up help on the function.  help quality varies a lot, but sometimes it tells you what the different elements of the output are.

init <- data %>%
  slice(1) %>%
  transmute(
    X = Infected / opt$par["N"],
    I = X * opt$par["I0_factor"],
    S = 1 - I - X,
    R = 0
  ) %>% 
  select(S, I, X, R) %>%
  unlist() %>%
  set_names(c("S", "I", "X", "R"))

new <- ode(y = init, times = new_times, func = eval(as.symbol(fun)), parms = opt$par, method = 'ode45') %>%
  as.data.frame() %>%
  map_dfc(as.numeric) %>%
  mutate_at(vars(-time), ~ .x * opt$par["N"])
