# Chase Abram

# Start log
row = ARGS[1]
row_str = string(row)

using XLSX
include("seq_auction_firm_dyn.jl")
# using Logging

# io = open("log_"*row_str*".txt", "w+")

run_date_time = read("output/run_date_time.txt", String)
out_name = "out_"*run_date_time
redirect_stdio(stdout="output/"*out_name*"/logs/log_"*row_str) do

    # Go to proper directory here

    ##
    # Load packages needed, as well as model setup code
    
    ##
    
    ##

    # Use row number of iteration to grab moments
    # asdf


    xf = XLSX.readxlsx("example_moments.xlsx")

    sh = xf["anzsic_period"]


    anzsic = sh["A"*row_str]
    period = sh["B"*row_str]
    e_to_u = 0.2 #sh["F"*row_str]
    switch_dwq = sh["H"*row_str]
    lp_pareto = 1.1 # -sh["M"*row_str]
    lp_drift = -sh["O"*row_str]

    # mcal = seq_auction_firm_dyn()
    mcal = calibrate_model([lp_drift, lp_pareto, e_to_u, switch_dwq])

    ##
    # Write parameters to spreadsheet

    # "output/"*out_name*"/logs/log_"*row_str
    XLSX.openxlsx("output/"*out_name*"/params/params_"*row_str*".xlsx", mode="w") do xfo
        sheet = xfo[1]
        XLSX.rename!(sheet, "params")
        sheet["A1"] = "anzsic"
        sheet["A2"] = anzsic

        sheet["B1"] = "period"
        sheet["B2"] = period

        sheet["C1"] = "lambda"
        sheet["C2"] = mcal.lambda

        sheet["D1"] = "delta"
        sheet["D2"] = mcal.delta

        sheet["E1"] = "chi"
        sheet["E2"] = mcal.chi

        sheet["F1"] = "xi"
        sheet["F2"] = mcal.xi


        # will add a row from "A5" to "E5"
        # sheet["A5"] = collect(1:5) # equivalent to `sheet["A5", dim=2] = collect(1:4)`

        # # will add a column from "B1" to "B4"
        # sheet["B1", dim=1] = collect(1:4)

        # # will add a matrix from "A7" to "C9"
        # sheet["A7:C9"] = [ 1 2 3 ; 4 5 6 ; 7 8 9 ]
    end


    # Other stuff to save, maybe in a batch-labelled folder?

end







