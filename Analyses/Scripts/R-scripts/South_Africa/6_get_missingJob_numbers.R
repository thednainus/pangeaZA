# I am using array jobs to run the search for the maximum likelihood trees
# However, several jobs passed the walltime and got terminated
# This script is to get the jobs that failed to complete and send them
# to the cluster using new PBS script that will allow saving the
# checkpoint files created by RAxML-NG

# directory in which ML files are located
ML_dir <- "~/Box Sync/my_R_packages/pangea/Analyses/Trees/South_Africa/env/ML/results"
ML_filenames <- Sys.glob(file.path(ML_dir, "*.bestTree"))

# gets only filenames
filename <- basename(ML_filenames)

# separated each file names by using the symbol "_"
completed_jobs <- strsplit(filename, "_")

# numbers oj jobs that completed
compl_job_num <- sort(as.numeric(unlist(lapply(completed_jobs, function(x) x[2]))))

# total number of jobs that I run in the cluster
total_jobs <- seq(1,50)

# get number of jobs to re-run
setdiff(total_jobs, compl_job_num)
