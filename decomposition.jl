# Decompose average earnings changes into search frictions and LP churn effects 

# Get rows of parameters
row0 = parse(Int, ARGS[1])
row1 = row0 + 1
row_str0 = string(row0)
row_str1 = string(row1)

using XLSX
include("seq_auction_firm_dyn.jl")

run_date_time = read("output/run_date_time.txt", String)
out_name = "out_"*run_date_time
redirect_stdio(stdout="output/"*out_name*"/logs/log_"*row_str0*"_decomposition") do

    # Go to proper directory here

    ##
    # Load packages needed, as well as model setup code
    
    ##
    
    ##

    # Use row number of iteration to grab moments
    # asdf


    xf = XLSX.readxlsx("output/"*out_name*"/params/params_"*row_str0*".xlsx")
    sh = xf["params"]

    anzsic = sh["A2"]
    period = sh["B2"]
    
    lambda0 = sh["C2"]
    delta0 = sh["D2"]
    chi0 = sh["E2"]
    xi0 = sh["F2"]

    xf = XLSX.readxlsx("output/"*out_name*"/params/params_"*row_str1*".xlsx")
    sh = xf["params"]

    lambda1 = sh["C2"]
    delta1 = sh["D2"]
    chi1 = sh["E2"]
    xi1 = sh["F2"]

    m0 = seq_auction_firm_dyn(lambda = lambda0, delta = delta0, chi = chi0, xi = xi0)
    mcmc!(m0)
    println("m0 done")
    flush(stdout)

    mlambda = seq_auction_firm_dyn(lambda = lambda1, delta = delta0, chi = chi0, xi = xi0)
    mcmc!(mlambda)
    println("mlambda done")
    flush(stdout)

    mchi = seq_auction_firm_dyn(lambda = lambda0, delta = delta0, chi = chi1, xi = xi0)
    mcmc!(mchi)
    println("mchi done")
    flush(stdout)

    m1 = seq_auction_firm_dyn(lambda = lambda1, delta = delta0, chi = chi1, xi = xi0)
    mcmc!(m1)
    println("m1 done")
    flush(stdout)
    
    ##
    # Write decomp to spreadsheet

    XLSX.openxlsx("output/"*out_name*"/decompositions/decompositions_"*row_str0*".xlsx", mode="w") do xfo
        sheet = xfo[1]
        XLSX.rename!(sheet, "decomps")
        sheet["A1"] = "anzsic"
        sheet["A2"] = anzsic

        sheet["B1"] = "period"
        sheet["B2"] = period

        sheet["C1"] = "stay_dw_0"
        dw0 = stay_mean_wage_change(m0)
        println("dw0: ", dw0)
        sheet["C2"] = dw0

        sheet["D1"] = "stay_dw_1"
        dw1 = stay_mean_wage_change(m1)
        println("dw1: ", dw1)
        sheet["D2"] = dw1

        println()

        sheet["E1"] = "stay_dw_lambda"
        dwlambda = stay_mean_wage_change(mlambda)
        println("dwlambda: ", dwlambda)
        sheet["E2"] = dwlambda

        sheet["F1"] = "stay_dw_chi"
        dwchi = stay_mean_wage_change(mchi)
        println("dwchi: ", dwchi)
        sheet["F2"] = dwchi

        println()
        println("lambda0 = ", lambda0, ", lambda1 = ", lambda1)

        sheet["G1"] = "lambda_share"
        lambda_share = (dwlambda - dw0)/(dw1 - dw0)
        println("lambda_share: ", lambda_share)
        sheet["G2"] = lambda_share

        sheet["H1"] = "chi_residual"
        chi_residual = (dw1 - dwlambda)/(dw1 - dw0)
        println("chi_residual: ", chi_residual)
        sheet["H2"] = chi_residual

        println()
        println("chi0 = ", chi0, ", chi1 = ", chi1)

        sheet["I1"] = "chi_share"
        chi_share = (dwchi - dw0)/(dw1 - dw0)
        println("chi_share: ", chi_share)
        sheet["I2"] = chi_share

        sheet["J1"] = "lambda_residual"
        lambda_residual = (dw1 - dwchi)/(dw1 - dw0)
        println("lambda_residual: ", lambda_residual)
        sheet["J2"] = lambda_residual
    end


    # Other stuff to save, maybe in a batch-labelled folder?

end




# Baseline

# Full change

# Frictions only

# LP churn only


# Output 

# Raw calcs of each of 4 above

# Share just frictions (residual is LP churn)

# Share just LP churn (residual is frictions)

# Do also for J2J?














