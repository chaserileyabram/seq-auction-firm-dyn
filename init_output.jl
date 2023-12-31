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
# batch_message = "No message"
batch_message = "Run with manufacturing industry moments from ABS. 
These moments are a rough first pass, and should be used with caution.
Trying with small nI and nT to check stuff works at scale with speed. Will move to properly large numbers after.
"
write("batch_message.txt", batch_message)

# Return to root
cd("../..")
