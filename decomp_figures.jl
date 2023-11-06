# Chase Abram
# Analyze decompositions across industries or sectors

# Output of interest
cd("output/out_2023-10-31T12:24:03.598/")

##
# Read params
xf = XLSX.readxlsx("params/params_combined.xlsx")
sh = xf["params"]

lambda0 = fill(NaN, 1000)
lambda1 = fill(NaN, 1000)
chi0 = fill(NaN, 1000)
chi1 = fill(NaN, 1000)
for i in 2:2:1000
    # println(i)
    ind = Int(i/2)
    l0 = sh["C"*string(i)]
    if !ismissing(l0)
        lambda0[ind] = l0
        chi0[ind] = sh["E"*string(i)]
    end

    l1 = sh["C"*string(i+1)]
    if !ismissing(l1)
        lambda1[ind] = l1
        chi1[ind] = sh["E"*string(i+1)]
    end
end

good_inds = .!isnan.(lambda0) .&& .!isnan.(lambda1)

lambda0 = lambda0[good_inds]
lambda1 = lambda1[good_inds]
chi0 = chi0[good_inds]
chi1 = chi1[good_inds]

bw = 0.05
plt = density(lambda0, label = L"\lambda_0", xlabel = L"\lambda", ylabel = "density", bandwidth=bw)
density!(lambda1, label = L"\lambda_1", bandwidth = bw)
display(plt)
savefig("/Users/chaseabram/UChiGit/seq-auction-firm-dyn/figures/lambda_dist.svg")

bw = 0.02
plt = density(chi0, label = L"\chi_0", xlabel = L"\chi", ylabel = "density", bandwidth=bw)
density!(chi1, label = L"\chi_1", bandwidth = bw)
display(plt)
savefig("/Users/chaseabram/UChiGit/seq-auction-firm-dyn/figures/chi_dist.svg")

##
bw = 0.07
plt = density(lambda1 - lambda0, label = false, xlabel = L"\Delta \lambda", ylabel = "density", bandwidth=bw)
display(plt)
savefig("/Users/chaseabram/UChiGit/seq-auction-firm-dyn/figures/dlambda_dist.svg")

bw = 0.07
plt = density(chi1 - chi0, label = false, xlabel = L"\Delta \chi", ylabel = "density", bandwidth=bw)
display(plt)
savefig("/Users/chaseabram/UChiGit/seq-auction-firm-dyn/figures/dchi_dist.svg")

scatter(lambda1 - lambda0, chi1 - chi0)
scatter(lambda0, chi0)
scatter(lambda1, chi1)

##
plt = scatter(lambda0, lambda1, 
xlabel = L"\lambda_0", ylabel = L"\lambda_1", label = L"(\lambda_0, \lambda_1)",
xlim = [0, 0.6], ylim = [0, 0.6])
line45 = LinRange(0,0.6,1000)
plot!(line45, line45, label = "45-degree")
display(plt)
savefig("/Users/chaseabram/UChiGit/seq-auction-firm-dyn/figures/lambda_scatter.svg")

plt = scatter(chi0, chi1, 
xlabel = L"\chi_0", ylabel = L"\chi_1", label = L"(\chi_0, \chi_1)",
xlim = [0, 0.25], ylim = [0, 0.25])
line45 = LinRange(0,0.25,1000)
plot!(line45, line45, label = "45-degree")
display(plt)
savefig("/Users/chaseabram/UChiGit/seq-auction-firm-dyn/figures/chi_scatter.svg")



##
# Read decomps
xf = XLSX.readxlsx("decompositions/decompositions_combined.xlsx")
sh = xf["decomps"]

# Get the lambda shares
lambda_shares = []
chi_shares = []
dw0 = []
dw1 = []
for i in 2:1000
    l_sh = sh["G"*string(i)]
    if !ismissing(l_sh)
        push!(lambda_shares, l_sh)
        push!(chi_shares, sh["I"*string(i)])

        push!(dw0, sh["C"*string(i)])
        push!(dw1, sh["D"*string(i)])
    end
end

inds = findall(abs.(lambda_shares) .< 3 .&& abs.(chi_shares) .< 3)
lambda_shares = lambda_shares[inds]
chi_shares = chi_shares[inds]

# filter!(x -> x > -3 && x < 3, lambda_shares)
# filter!(x -> x > -3 && x < 3, chi_shares)

# plt = histogram(lambda_shares, normalize=true)
# display(plt)
bw = 0.2
plt = density(chi_shares, normalize=true, label = L"\chi"*" share", bandwidth = bw, 
xlim = [-0.1,1.1], xlabel = "Share of change explained", ylabel = "density")
density!(1 .- lambda_shares, normalize=true, label = L"\lambda"*" residual share", bandwidth = bw)
display(plt)
savefig("/Users/chaseabram/UChiGit/seq-auction-firm-dyn/figures/chi_share.svg")


##

plt = scatter(lambda_shares, 1 .- chi_shares, label = "lambda share vs residual from chi")
display(plt)

##

plt = density(dw1 - dw0, normalize=true)
display(plt)

plt = scatter(dw0, lambda_shares)
display(plt)

plt = scatter(dw1, lambda_shares)
display(plt)

plt = scatter(dw1 - dw0, lambda_shares)
display(plt)





