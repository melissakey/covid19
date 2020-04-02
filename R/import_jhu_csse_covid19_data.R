clean_jhd_to_long <- function(df) {
  df_str <- deparse(substitute(df))
  var_str <- substr(df_str, 1, str_length(df_str) - 4)
  
  df %>% 
    select(-`Province/State`, -Lat, -Long) %>%
    rename(country = `Country/Region`) %>%
    mutate(iso3c = countrycode(country,
                               origin = "country.name",
                               destination = "iso3c")) %>%
    select(-country) %>%
    filter(!is.na(iso3c)) %>%
    group_by(iso3c) %>%
    summarise_at(vars(-group_cols()), sum) %>% 
    pivot_longer(
      -iso3c, 
      names_to = "date_str", 
      values_to = var_str
    ) %>%
    ungroup() %>%
    mutate(date = mdy(date_str)) %>%
    select(iso3c, date, !! sym(var_str)) 
}