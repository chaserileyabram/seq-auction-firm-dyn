# Combine outputs into a single spreadsheet

using XLSX
# Go into output folder
cd("output")

# Go into out_XXX folder
run_date_time = read("run_date_time.txt", String)
out_name = "out_"*run_date_time
cd(out_name)


# Initialize combined output sheet
XLSX.openxlsx("params_combined.xlsx", mode="w") do xfo
    sheet = xfo[1]
    XLSX.rename!(sheet, "params")
    sheet["A1"] = "anzsic"
    sheet["B1"] = "period"
    sheet["C1"] = "lambda"
    sheet["D1"] = "delta"
    sheet["E1"] = "chi"
    sheet["F1"] = "xi"
end

cd("params")
# Get names of successful output
# files = sort(readdir())
files = filter(f -> f[1:6] == "params", readdir())
println("files: ", files)
cd("..")

for i in 1:length(files)
    file = files[i]
    row = XLSX.readdata("params/"*file, "params!A2:F2")

    XLSX.openxlsx("params_combined.xlsx", mode="rw") do xfo
        sheet = xfo[1]
        sheet["A"*string(i+1)] = row
    end
end


# Message about output?
# write("run_date_time.txt", run_date_time)

# Return to root
cd("../..")



