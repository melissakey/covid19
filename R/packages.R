# Standard 
library(drake)
library(tidyverse)
library(conflicted)
library(here)
library(tidyselect)
library(magrittr)
library(glue)
library(rjson)

# Plotting
library(scales)
library(ggthemes)
library(patchwork)

# Analysis
library(deSolve)
library(minpack.lm)


## conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("which", "base")
conflict_prefer("ggsave", 'ggplot2')
conflict_prefer("discard", "purrr")
conflict_prefer("summarize", "dplyr")