# Combine outputs into a single spreadsheet

using XLSX
# Go into output folder
cd("output")

# Go into out_XXX folder
run_date_time = read("run_date_time.txt", String)
out_name = "out_"*run_date_time
cd(out_name)


# Initialize combined output sheet
XLSX.openxlsx("params/params_combined.xlsx", mode="w") do xfo
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
filter!(f -> f != "params_combined.xlsx", files)
sort!(files)
println("files: ", files)
cd("..")

for i in 1:length(files)

    file = files[i]

    ind = parse(Int, file[8:end-5])

    row = XLSX.readdata("params/"*file, "params!A2:F2")

    XLSX.openxlsx("params/params_combined.xlsx", mode="rw") do xfo
        sheet = xfo[1]
        sheet["A"*string(ind)] = row
    end
end


# Initialize combined decomps sheet
XLSX.openxlsx("decompositions/decompositions_combined.xlsx", mode="w") do xfo
    sheet = xfo[1]
    XLSX.rename!(sheet, "decomps")
    sheet["A1"] = "anzsic"
    sheet["B1"] = "period"
    sheet["C1"] = "stay_dw_0"
    sheet["D1"] = "stay_dw_1"
    sheet["E1"] = "stay_dw_lambda"
    sheet["F1"] = "stay_dw_chi"
    sheet["G1"] = "lambda_share"
    sheet["H1"] = "chi_residual"
    sheet["I1"] = "chi_share"
    sheet["J1"] = "lambda_residual"
end

cd("decompositions")
# Get names of successful output
# files = sort(readdir())
files = filter(f -> f[1:14] == "decompositions", readdir())
filter!(f -> f != "decompositions_combined.xlsx", files)
sort!(files)
println("files: ", files)
cd("..")

for i in 1:length(files)

    file = files[i]

    ind = parse(Int, file[16:end-5])

    row = XLSX.readdata("decompositions/"*file, "decomps!A2:J2")

    XLSX.openxlsx("decompositions/decompositions_combined.xlsx", mode="rw") do xfo
        sheet = xfo[1]
        sheet["A"*string(ind)] = row
    end
end


# Message about output?
# write("run_date_time.txt", run_date_time)

# Return to root
cd("../..")



