divide_and_conquor <- function(df, status = status, pool_size, depth = 0) {
  # browser()
  if(depth == 0) {
    return(nrow(df))
  } else {
    qs <- enquo(status)
    N <- nrow(df)
    
    n_pools = N %/% pool_size + (N %% pool_size > 0)
    
    round_results <- df %>%
      mutate(
        pool  = sample(rep(1:n_pools, pool_size), N) 
      ) %>%
      nest(data = -pool) %>%
      mutate(
        test_result = map_lgl(data, ~ 
                                .x %>%
                                filter(!!qs > 0) %>%
                                nrow() %>%
                                magrittr::is_greater_than(0)
        ),
        n = map_int(data, ~ nrow(.x))
      )
    
    n_round <- nrow(round_results)
    
    sum_results <- round_results %>%
      filter(test_result) %$%
      map2_int(data, n, ~ divide_and_conquor(.x, status = status, pool_size = .y %/% 2 + (.y %% 2 > 0), depth = depth - 1)) %>%
      sum() %>%
      add(n_round)
  }
  return(sum_results)
}
