library(here)
source(here("packages.R"))
source(here("functions.R"))

# Set overall targets options
tar_option_set(
  seed = 202411,
  format = "parquet",
  error = "null",
)

# Set parallel options depending on hostname
if (Sys.info()["nodename"] == "login001.meerkat.mcri.edu.au") {
  tar_option_set(
    controller = crew_controller_slurm(
      name = "misscompl",
      workers = 160,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = 2,
        time_minutes = 24 * 60,
        partition = "prod_med",
        script_lines = "
#SBATCH --account=cebu
module load r/4.4.1
"
      )
    )
  )
} else {
  tar_option_set(
    controller = crew_controller_local(workers = 4)
  )
}

# for a subset of scenarios:
# scenario_list <- scenario_list %>%
#   filter(
#     missingness_mechanism == "mar_strong",
#     sample_size == "large",
#     treatment_effect == "power80"
#   )

# Define targets
list(
  tar_target(
    dgp_params,
    get_dgp_params(scenario_list),
    deployment = "main"
  ),
  tar_target(
    scenario_params,
    get_scenario_params(scenario_list, dgp_params),
    deployment = "main"
  ),

  tar_map_rep(
    sim_reps,
    command =
      run_scenario(filter(scenario_params, scenario_name == name)),
    values = scenario_list %>%
      select(name = scenario_name),
    names = any_of("name"),
    batches = 100,
    reps = 100,
    combine = TRUE,
  ),

  # Generate summary data frame of indiviudal replicate results
  tar_target(
    sim_reps_summary,
    summarise_sim_reps(sim_reps),
    deployment = "main"
  )

)
