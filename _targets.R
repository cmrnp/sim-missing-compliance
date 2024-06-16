library(targets)
library(tarchetypes)
library(crew)
library(tidyverse)
library(glue)
library(marginaleffects)
library(insight)
library(estimatr)
library(broom)
library(mice)
library(withr)
library(sandwich)
library(cowplot)

source("R/data-generation.R")
source("R/estimators.R")
source("R/missingness-methods.R")
source("R/scenarios.R")
source("R/summaries.R")

# Set overall targets options
tar_option_set(
  seed = 202406,
  # error = "null"
)

# Set parallel options depending on hostname
if (Sys.info()["nodename"] == "dev2") {
  tar_option_set(
    controller = crew_controller_pbs(
      workers = 144,
      pbs_walltime_hours = 168,
      script_lines = "
#PBS -q batch
#PBS -A CEBU1
module load R/4.3.2
"
    )
  )
} else {
  tar_option_set(
    controller = crew_controller_local(workers = 4)
  )
}

# Define targets
list(
  tar_target(dgp_params, get_dgp_params(scenario_list)),
  tar_target(scenario_params, get_scenario_params(scenario_list, dgp_params)),
  tar_map_rep(
    sim_reps,
    command = 
      run_scenario(filter(scenario_params, scenario_name == name)),
    values = scenario_list %>%
#      filter(outcome_missingness == "no", sample_size == "small",
#             missingness_mechanism %in% c("none", "mar_strong")) %>%
      select(name = scenario_name),
    names = any_of("name"),
    batches = 20,
    reps = 100,
    combine = TRUE,
  ),
  tar_target(sim_reps_summary, summarise_sim_reps(sim_reps)),
  tar_target(
    save_results_no_missing, 
    plot_results_no_missing(sim_reps_summary), 
    format = "file"
  ),
  tar_target(
    save_results_null, 
    plot_results_null(sim_reps_summary), 
    format = "file"
  ),
  tar_target(
    save_results_main, 
    plot_results_main(sim_reps_summary), 
    format = "file"
  )
)
