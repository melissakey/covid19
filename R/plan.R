plan <- drake_plan(
  data = load_data(set = "china"),
  initial_plot = plot_infections(data),
  
  # SIR_predictions = predict_infections(data, new_times = 1:70, fun = "SIR"),
  SIRX_confirmed = predict_infections(data, new_times = 1:70, fun = "SIRX"),
  
  # SIR_prediction_plot = pmap(list(data$Province, data$data, SIR_predictions$fit), plot_predictions),
  SIRX_prediction_plot = pmap(list(data$Province, data$data, SIRX_confirmed$fit), plot_predictions),
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
