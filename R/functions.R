
# from https://www.r-bloggers.com/epidemiology-how-contagious-is-novel-coronavirus-2019-ncov/
SIR <- function(time, N, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

# equivalent to dxdt equations in SIRXconfirmedModel code, with parameters changed to match paper.
# the shutdown/quarantine models are simplifications of this, with kappa0 or kappa set to 0.
SIRX <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- - alpha * S * I - kappa0 * S 
    dI <- + alpha * S * I - beta * I - kappa0 * I - kappa * I
    dX <- + kappa * I + kappa0 * I
    dR <- + kappa0 * S
    list(c(dS, dI, dX, dR))
  })
}


# combines "residual" and "SIRX" functions in SIRXconfirmedModel code.
RSS <- function(parameters, data, init, fun, N) {
  if(fun == 'SIR') {                      # SIR model (directly parallels the SIRX code, but can be ignored.
    init <- data %>%
      slice(1) %>%
      transmute(
        S = N - Infected,
        I = Infected,
        R = 0
      ) %>%
      unlist()
    
    fit <- ode(
      y = init,
      times = data$day,
      func = eval(as.symbol(fun)),
      parms = parameters,
      method = 'ode45',
      N = N
    ) %>%
      as_tibble() %>%
      map_dfc(as.numeric)
    
    res <- fit$I - data$Infected
  } else {
    # set initial values of S, I, R, and X
    init <- data %>%
      filter(Infected > 0) %>%
      slice(1) %>%                         # keep only the first row of the data set
      transmute(                           # create the following columns, and keep only these columns
        X = Infected / parameters["N"],
        I = X * parameters["I0_factor"],
        S = 1 - I - X,
        R = 0
      ) %>% 
      select(S, I, X, R) %>%               # reorders the columns to match (S, I, X, R) format in SIRX function
      unlist() %>%                         # removes the dataframe structure from the result to produce a vector with the 4 terms
      set_names(c("S", "I", "X", "R"))     # fix names of the 4 elements since these did not ocome through correctly.
    
    fit <- ode(                            # wrapper around the various ordinary differntial equation solvers in R
      y = init,                            # the initial state values
      times = data$day,                    # time sequence for which output is wanted
      func = eval(as.symbol(fun)),         # either SIR or SIRX - eval(as.symbol(.)) translates the string to a function call
      parms = parameters,                  # the parameters being optimized
      method = 'ode45'                     # this is the same implementation as in the original code
    ) %>%
      as_tibble() %>%                      # these 2 lines format the output in a convienent way
      map_dfc(as.numeric)
    
    res <- parameters["N"] * fit$X - data$Infected  # calculate the residuals
  }

  res
}


# the functions below are used as part of a pipeline-management toolkit.
# You can ignore them - they use complex data structures in R and will probably be difficult to interpret unless you know R.


load_data <- function(set){
  # for now, let's grab the initial China data
  if(set == 'simple')
    data = tibble(
      country = 'China',
      Infected = c(45, 62, 121, 198, 291, 440, 571, 830, 1287, 1975, 2744, 4515, 5974, 7711, 9692, 11791, 14380, 17205, 20440),
      day = 1:length(Infected),
      N = 1400000000
    ) %>%
      nest(data = c(day, Infected))
  if(set == 'china') {
    data <- readLines("data/all_confirmed_cases_with_population.json") %>%
      paste(collapse = "") %>%
      fromJSON() %>%
      map_df(~ .x %>%
               as_tibble() %>%
               rename(
                 N = population,
                 day = times,
                 Infected = cases
               ) %>%
               nest(data = -N),
             .id = 'Province'
      )
  }
  data
  
}
plot_infections <- function(data, group) {
  group = enquo(group)
  p1 <- data %>%
    unnest(data) %>%
    ggplot(aes(day, Infected, group = !!group)) + 
    geom_point() + 
    theme_minimal() + 
    labs(
      x = "Day",
      y = "Infected"
    )
  p1 + 
    geom_line() + 
    p1 + 
    geom_smooth(method = 'lm', se = FALSE) + 
    scale_y_log10() + 
    plot_annotation(title = "Confirmed Cases 2019-nCoV China")
}

predict_infections <- function(data, fun, new_times = 1:70, I0_factor = 10, R0 = 6, beta = 1/8, method = 'confirmed') {
  # browser()
  data %>%
    mutate(
      # 
      opt = map2(data, N,
        function(data, N, fun = fun, R0 = 6.2, beta = 1/8, I0_factor = 10, method = c("confirmed", "shutdown", "quarantine")) {
          if(fun == 'SIR') {
            # set parameters/ranges for SIR model.
            # These follow the code from the website above.  Note that this beta is equivalent to alpha in the SIRX model.
              lower <- c(0, 0)
              parms  <- c(beta = .5, gamma = .5)
              upper <- c(1, 1)

          } else if(fun == 'SIRX') {
            # set initial parameters/ranges for the SIRX model.
            # these should match the initial values defined in the fit functions.
            # alpha == eta (rate from S -> I)
            # beta  == rho (rate from I -> R)
            # kappa0       (rate from I+S -> X; general quarantine rate)
            # kappa        (rate from I -> X; higher quarantine rate for infected individuals)
            method = match.arg(method)
            
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
          }

          # minimizer using the Levenberg-Marquardt algorithm (same as in python code)
          nls.lm(
            par = parms, 
            fn = RSS,
            lower = lower,
            upper = upper,
            control = nls.lm.control(maxiter = 1000, maxfev = 5000),
            data = data,
            fun = fun,
            N = N)
        }, 
        fun = fun, R0 = R0, beta = beta, I0_factor = I0_factor, method = method
      ),
      fit = pmap(list(data, N, opt),
        function(data, N, opt, new_times, I0_factor, fun) {
          if(fun == 'SIR') {
            init <- data %>%
              slice(1) %>%
              transmute(
                S = N - Infected,
                I = Infected,
                R = 0
              ) %>%
              unlist()
            pred <- ode(y = init, times = new_times, func = eval(as.symbol(fun)), parms = opt$par, N = N, method = 'ode45') %>%
              as.data.frame() %>%
              map_dfc(as.numeric)
          } else {
            init <- data %>%
              filter(Infected > 0) %>%
              slice(1) %>%
              transmute(
                X = Infected / N,
                I = X * I0_factor,
                S = 1 - I - X,
                R = 0
              ) %>% 
              select(S, I, X, R) %>%
              unlist() %>%
              set_names(c("S", "I", "X", "R"))
            
            pred <- ode(y = init, times = new_times, func = eval(as.symbol(fun)), parms = opt$par, method = 'ode45') %>%
              as.data.frame() %>%
              map_dfc(as.numeric) %>%
              mutate_at(vars(S, I, R, X), ~ . * opt$par["N"])
          }
          
          pred

        }, new_times = new_times, fun = fun, I0_factor = I0_factor)
    ) %>%
    select(matches("country"), matches("Province"), opt, fit)
  
}

plot_predictions <- function(country, data, predictions) {
  p1 <- predictions %>%
    select(
      time,
      Infected = I,
      # Susceptible = S,
      Recovered = R,
      # "Symp. Quarantined" = matches("X")
      ) %>%
    pivot_longer(cols = -time, names_to = "status", values_to = 'count') %>%
    ggplot(aes(time, count, color = status)) + 
    labs(title = country, x = 'Day') + 
    geom_line(size = 1.5) +
    theme_minimal() + 
    geom_point(data = data, 
      aes(x = day, y = Infected),
      inherit.aes = FALSE,
      size = 2,
      shape = 1
    )
  p1
}
