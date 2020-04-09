source(here::here("R/packages.R"))
source(here::here("R/divide_and_conquor.R"))

n_per_pool <- 12

prob_pos <- 0.06

reps <- 100
N <- 1200

parms_df <- list(
  n_per_pool = c(8, 12, 16, 20),
  N = c(2000, 5000, 10000),
  prob_positive = c(.02, .05, .08),
  depth = 1:3
) %>%
  cross_df()

results <- parms_df %>%
  # slice(1:20) %>%
  mutate(.,
         result = pmap(., function(n_per_pool, N, prob_positive, depth = depth) {
           map_dbl(1:reps, ~ {
             test_population <- tibble(
               ind = 1:N,
               status =  rbinom(N, 1, prob_pos),                   # generate N samples, with prob_pos of a positive tests
             )
             
             res <- divide_and_conquor(test_population, status, pool_size = n_per_pool, depth = depth)
             res
           })
         }),
         mean = map_dbl(result, ~ mean(.x)),
         sd = map_dbl(result, ~ sd(.x))
  )
# identity()

# generate N samples, with prob_pos of a positive tests

# pool at random
p1 <- results %>% 
  group_by(prob_positive) %>%
  ggplot(
    aes(
      n_per_pool,
      mean /N,
      color = factor(depth),
      shape = factor(prob_positive)
    )
  ) +
  geom_point(position = position_dodge(width = .5)) +
  geom_linerange(aes(ymin = (mean - 2 * sd) / N, ymax = (mean + 2*sd) / N), position = position_dodge(width = .5)) + 
  labs(
    shape = "Expected proportion of positive results",
    color = "Maximum number of batch tests per sample",
    x = "Expected # of Tests / Total # Samples",
    y = "# of Samples per Pool"
  ) + 
  theme_bw() +
  theme(legend.position = 'bottom') + 
  facet_wrap(~ N) + 
  scale_color_few() + 
  NULL

p2 <- results %>% 
  group_by(prob_positive) %>%
  ggplot(
    aes(
      n_per_pool,
      mean,
      color = factor(depth),
      shape = factor(prob_positive),
      # group = interaction(prob_positive, depth)
    )
  ) +
  geom_point(position = position_dodge(width = .5)) +
  geom_linerange(aes(ymin = (mean - 2 * sd), ymax = (mean + 2*sd)), position = position_dodge(width = .5)) + 
  theme_bw() +  labs(
    shape = "Expected proportion of positive results",
    color = "Maximum number of batch tests per sample",
    x = "Expected # of Tests",
    y = "# of Samples per Pool"
  ) + 
  theme(legend.position = 'bottom') + 
  facet_wrap(~ N, scales = "free_y") + 
  scale_color_few() + 
  NULL
