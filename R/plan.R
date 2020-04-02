plan <- drake_plan(
  china_data = load_data(set = "china"),
  usa_data = load_data(set = 'USA'),
  `usa_data_2020-04-01` = load_data(set = 'USA') %>%
    mutate(
      data = map(data, ~ .x %>%
                   mutate(new = date > as.Date("2020-03-29")))
    ),
  
  # Initial plots
  initial_plot_china = plot_infections(china_data, Province),
  initial_plot_usa = plot_infections(usa_data, x = date, state_name),
  
  
  # SIR_predictions = predict_infections(data, new_times = 1:70, fun = "SIR"),
  SIRX_confirmed_china = predict_infections(china_data, new_times = 1:70, fun = "SIRX"),
  SIRX_confirmed_usa = predict_infections(usa_data, new_times = 1:70, fun = "SIRX"),
  # SIRX_quarantined = predict_infections(data, new_times = 1:70, fun = "SIRX", method = 'quarantine'),
  
  
  # SIR_prediction_plot = pmap(list(data$Province, data$data, SIR_predictions$fit), plot_predictions),
  SIRX_plot_china = pmap(list(china_data$Province, china_data$data, SIRX_confirmed_china$fit), plot_predictions),
  SIRX_plot_usa = pmap(list(`usa_data_2020-04-01`$state_name, `usa_data_2020-04-01`$data, SIRX_confirmed_usa$fit), plot_predictions, group = new) %>%
    set_names(usa_data$state_abb),
  # 
  # SIRX_stats = SIRX_predictions %$%
  #   map2_df(Province, opt, ~ .y$par %>%
  #     as.matrix() %>%
  #     t() %>%
  #     as_tibble() %>%
  #     transmute(
  #       Province = .x,
  #       P = kappa0 / (kappa0 + kappa),
  #       Q = (kappa0 + kappa) / (beta + kappa0 + kappa),
  #       R0_free = alpha / beta,
  #       R0_eff = alpha / (beta + kappa0 + kappa)
  #     )
  # )
)
