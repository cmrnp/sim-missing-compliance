library(here)
source(here("packages.R"))
source(here("functions.R"))

# Set overall targets options
tar_option_set(
  seed = 202407,
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
#PBS -A cebu1
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
      filter(missingness_mechanism == "none") %>%
      #filter(outcome_missingness == "no", sample_size == "small",
      #       missingness_mechanism %in% c("none", "mar_strong")) %>%
      select(name = scenario_name),
    names = any_of("name"),
    batches = 20,
    reps = 100,
    combine = TRUE,
  ),
  
  # Generate summary data frame of indiviudal replicate results
  tar_target(sim_reps_summary, summarise_sim_reps(sim_reps)),
  
  # Plot results for no-missing-data scenarios
  tar_target(
    plot_results_no_missing, 
    make_plot_results_no_missing(sim_reps_summary)
  ),
  tar_target(
    save_results_no_missing, 
    save_plot_results_no_missing(plot_results_no_missing), 
    format = "file"
  ),

  # null treatment effect
  tar_target(
    plot_results_null, 
    make_plot_results_null(sim_reps_summary), 
  ),
  tar_target(
    save_results_null, 
    save_plot_results_null(plot_results_null), 
    format = "file"
  ),
  
  # Plot results for all scenarios with missing data
  tar_map(
    expand_grid(
      sample_size = c("small", "large"),
      outcome_missingness = c("no", "yes"),
    ),
    
    # non-null treatment effect
    tar_target(
      plot_results_main, 
      make_plot_results_main(sim_reps_summary, outcome_missingness, sample_size),
    ),
    tar_target(
      save_results_main, 
      save_plot_results_main(plot_results_main, outcome_missingness, sample_size), 
      format = "file"
    )
  )
)
