# Create output folder and subfolders
using Dates
# To output folder
cd("output")

# Save run id for use by other scripts
run_date_time = string(now())
write("run_date_time.txt", run_date_time)

# Make folder for this batch of output
out_str = "out_"*run_date_time
mkdir(out_str)
cd(out_str)

# Make log and parameter folders
mkdir("logs")
mkdir("params")
mkdir("decompositions")

# Message about batch?
# write("run_date_time.txt", run_date_time)

# Return to root
cd("../..")
