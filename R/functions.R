
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

# equivalent to dxdt equations in SIRXconfirmedModel code, with parameters changed to match standard notation.
# why oh why do they use alpha/beta in their paper when standard is beta/gamma?????
# (beta = eta/alpha, gamma = rho/beta)
# the shutdown/quarantine models are simplifications of this, with kappa0 or kappa set to 0.
SIRX <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- - beta * S * I - kappa0 * S                           # S = susceptible population *able to get infected*
    dI <- + beta * S * I - gamma * I - kappa0 * I - kappa * I   # I = infected individuals *able to infect other people* 
    dX <- + kappa * I + kappa0 * I                              # X = infected individuals *removed from general population/quarantined*
    dH <- + kappa0 * S                                          # H = healthy, quarantined individuals (protected from getting infected)
    dR <- + gamma * I                                           # R = recovered/removed population (cannot get reinfected)
    list(c(dS, dI, dX, dH, dR))
  })
}
SIRX.quarantine <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- - beta * S * I 
    dI <- + beta * S * I - gamma * I - kappa * I
    dX <- + kappa * I
    dH <- 0
    dR <- gamma * I
    list(c(dS, dI, dX, dH, dR))
  })
}

# combines "residual" and "SIRX" functions in SIRXconfirmedModel code.
RSS <- function(parameters, data, init, fun, N, SIRX.method = 'confirmed') {
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
    if(SIRX.method == 'quarantine') fun <- 'SIRX.quarantine'
    # set initial values of S, I, R, and X
    init <- data %>%
      filter(Infected > 0) %>%
      slice(1) %>%                         # keep only the first row of the data set
      transmute(                           # create the following columns, and keep only these columns
        X = Infected / parameters["N"],
        I = X * parameters["I0_factor"],
        S = 1 - I - X,
        H = 0,
        R = 0
      ) %>% 
      select(S, I, X, H, R) %>%               # reorders the columns to match (S, I, X, R) format in SIRX function
      unlist() %>%                         # removes the dataframe structure from the result to produce a vector with the 4 terms
      set_names(c("S", "I", "X", "H", "R"))     # fix names of the 4 elements since these did not ocome through correctly.
    
    # browser()
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
  # data sets
  data <- switch(set,
    "simple" = tibble(
      country = 'China',
      Infected = c(45, 62, 121, 198, 291, 440, 571, 830, 1287, 1975, 2744, 4515, 5974, 7711, 9692, 11791, 14380, 17205, 20440),
      day = 1:length(Infected),
      N = 1400000000
    ) %>%
      nest(data = c(day, Infected)),
    "china" = readLines("data/all_confirmed_cases_with_population.json") %>%
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
      ),
    "USA" = {
      pop_data <- read.xlsx("https://www2.census.gov/programs-surveys/popest/tables/2010-2019/state/totals/nst-est2019-01.xlsx", startRow = 4) %>%
        as_tibble() %>%
        filter(
          !is.na(Census)
        ) %>%
        select(Region = X1, N = `2019`) %>%
        mutate(
          state_name = str_remove(Region, "\\."),
          state_abb = datasets::state.abb[match(state_name, datasets::state.name)] %>%
            replace(state_name == 'Puerto Rico', 'PR') %>%
            replace(state_name == 'District of Columbia', "DC"),
          state_region = datasets::state.region[match(state_name, datasets::state.name)] %>%
            as.character() %>%
            replace(state_name == 'District of Columbia', "South")
        ) %>% 
        filter(!is.na(state_abb))
      
      covid <- get_states_daily() %>%
        full_join(pop_data, ., by = c("state_abb" = 'state')) %>%
        mutate(
          state_name = state_name %>%
            replace(state_abb == 'AS', "American Samoa") %>%
            replace(state_abb == 'GU', "Guam") %>%
            replace(state_abb == 'MP', "Mariana Islands") %>%
            replace(state_abb == 'VI', 'Virgin Islands'),
          state_region = state_region %>%
            replace(state_abb %in% c('AS', 'GU', 'MP', 'VI', 'PR'), "Territories"),
          N = N %>%
            replace(state_abb == 'GU', 167772) %>%
            replace(state_abb == 'VI', 104578) %>%
            replace(state_abb == 'MP', 51994) %>%
            replace(state_abb == 'AS', 55689)
        )
      covid_small <- covid %>%
        select(state_region, N, state_name, state_abb, date, Infected = positive) %>%
        group_by(state_name) %>%
        arrange(date) %>%
        fill(Infected) %>%
        mutate(day = 1:n()) %>%
        nest(data = c(date, day, Infected))
      
      data = covid %>%
        select(-Region) %>%
        nest(raw_data = -c(state_region, N, state = state_name, state_abb)) %>%
        left_join(covid_small) %>%
        mutate(
          max_cases = map_int(data, ~ max(.x$Infected))
        ) %>%
        filter(max_cases > 20)
    }
    # 'world' = {
    #   world_pop <- read.xlsx("http://gapm.io/dl_pop", sheet = 7) %>%
    #     as_tibble() %>%
    #     select(geo, name, indicator, contains("2020"))
    #   
    #   dat <- read_csv("https://covid.ourworldindata.org/data/ecdc/full_data.csv")
    # }
  )
  data
}
plot_infections <- function(data, x = day, group) {
  group = enquo(group)
  x = enquo(x)
  
  p1 <- data %>%
    unnest(data) %>%
    group_by(!!group) %>%
    filter(max(Infected) > 20) %>%
    ggplot(aes(!!x, Infected, group = !!group, color = !!group)) + 
    geom_point() + 
    geom_line() + 
    theme_minimal() + 
    labs(
      x = "Day",
      y = "Infected"
    )
  p1 + 
    geom_line() + 
    p1 + 
    geom_line() + 
    # geom_smooth(method = 'lm', se = FALSE) + 
    scale_y_log10()
}

predict_infections <- function(data, fun, new_times = 1:70, I0_factor = 10, R0 = 6.2, gamma = 1/8, method = 'confirmed') {
  data %>%
    mutate(
      # 
      opt = pmap(list(data, N),
        function(data, N, fun = fun, R0 = 6.2, gamma = 1/8, I0_factor = 10, method = c("confirmed", "shutdown", "quarantine")) {
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
            # beta  == rho (rate from S -> "R"; presumably "healthy quarantined?"  not clear!!) 
            # kappa0       (rate from I+S -> X; general quarantine rate)
            # kappa        (rate from I -> X; higher quarantine rate for infected individuals)
            method = match.arg(method)
            
            parms <- switch(method,
              'confirmed'  = c(gamma = gamma, beta = R0 * gamma, kappa = gamma, kappa0 = gamma / 2, I0_factor = I0_factor, N = N),
              'shutdown'   = c(gamma = gamma, beta = R0 * gamma, kappa =    0, kappa0 =     gamma, I0_factor = I0_factor, N = N),
              'quarantine' = c(gamma = gamma, beta = R0 * gamma, kappa = gamma, kappa0 =        0, I0_factor = I0_factor, N = N)
            )
            
            lower <- switch(method,
              'confirmed'  = c(gamma = gamma, beta = R0 * gamma, kappa = 0, kappa0 = 0, I0_factor = 0.001, N = N),
              'shutdown'   = c(gamma = gamma, beta = R0 * gamma, kappa = 0, kappa0 = 0, I0_factor = 1    , N = N),
              'quarantine' = c(gamma = gamma, beta = R0 * gamma, kappa = 0, kappa0 = 0, I0_factor = 0.001, N = N)
            )
            
            upper <-  switch(method,
              'confirmed'  = c(gamma = gamma, beta = R0 * gamma, kappa = Inf, kappa0 = Inf, I0_factor = Inf, N = {if(exists("N")) N else 115000000}),
              'shutdown'   = c(gamma = gamma, beta = R0 * gamma, kappa = 0,   kappa0 = Inf, I0_factor = 30, N = {if(exists("N")) N else 115000000}),
              'quarantine' = c(gamma = gamma, beta = R0 * gamma, kappa = Inf, kappa0 = 0,   I0_factor = Inf, N = {if(exists("N")) N else 1000000000})
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
            N = N,
            SIRX.method = method)
        }, 
        fun = fun, R0 = R0, gamma = gamma, I0_factor = I0_factor, method = method
      ),
      fit = pmap(list( data, N, opt),
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
                X = Infected / opt$par["N"],
                I = X * opt$par["I0_factor"],
                S = 1 - I - X,
                H = 0,
                R = 0
              ) %>% 
              select(S, I, X, H, R) %>%
              unlist() %>%
              set_names(c("S", "I", "X", "H", "R"))
            
            pred <- ode(y = init, times = new_times, func = eval(as.symbol(fun)), parms = opt$par, method = 'ode45') %>%
              as.data.frame() %>%
              map_dfc(as.numeric) %>%
              mutate_at(vars(S, I, R, X), ~ . * opt$par["N"])
          }
          
          pred

        }, new_times = new_times, fun = fun, I0_factor = I0_factor)
    ) %>%
    select(matches("country", ignore.case = TRUE), matches("Province", ignore.case = TRUE), contains("state", ignore.case = TRUE), opt, fit)
  
}

plot_predictions <- function(country, data, predictions, group = NULL) {
  
  group <- enquo(group)
  p1 <- predictions %>%
    select(
      time,
      Infected = I,
      # Susceptible = S,
      # Recovered = R,
      "Symp. Quarantined" = matches("X")
      ) %>%
    pivot_longer(cols = -time, names_to = "status", values_to = 'count') %>%
    ggplot(aes(time, count, color = status)) + 
    labs(title = country, x = 'Day') + 
    geom_line(size = 1.5) +
    theme_minimal() + 
    geom_point(data = data, 
      aes(x = day, y = Infected, shape = !!group),
      inherit.aes = FALSE,
      size = 2
    ) + 
    scale_shape_manual(values = c(16, 1))
  p1
}
