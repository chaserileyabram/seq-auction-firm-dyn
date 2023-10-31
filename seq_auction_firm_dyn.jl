# Chase Abram

# Sequential auction labor search with firm producutivity dynamics
# cd("/Users/chaseabram/Library/CloudStorage/Dropbox/Income Inequality Origins")
##
# Packages 
using LinearAlgebra
using Distributions
using Parameters
using QuadGK
using NLsolve
using Roots

using Plots
default(linewidth = 5.0, 
legend = true, 
size = (800,600))
using LaTeXStrings

using StatsPlots

using Random

using KernelDensity

using Statistics
using StatsBase

using DataFrames
using FixedEffectModels

##
@with_kw mutable struct seq_auction_firm_dyn
    
    # Preferences
    rho = 0.05
    # u = (x -> log(x))
    # up = (x -> 1/x)
    # uinv = (x -> exp(x))
    u = (x -> x)
    up = (x -> 1)
    uinv = (x -> x)

    # offer rate
    lambda = 0.2
    # separation rate
    delta = 0.05
    
    # Firm productivity
    # arrival rate
    chi = 0.1
    # z grid points
    n_z = 15
    
    # # mean productivity
    # mu_z = 10
    # # standard dev of shock
    # sigma_z = 1.0
    # # autocorrelation coef 
    # rho_z = 0.8
    # # transition matrix
    # zgrid = Array(rouwenhorst(n_z, mu_z, sigma_z, rho_z)[1])
    # ztrans = rouwenhorst(n_z, mu_z, sigma_z, rho_z)[2]

    # drift
    mu_z = -0.01
    # noise
    sigma_z = 0.05
    xi = 2.0
    # transition matrix
    # zgrid = Array(kesten(n_z, mu_z, sigma_z)[1])
    zgrid = Array(kesten(n_z, xi)[1])
    # ztrans = kesten(n_z, mu_z, sigma_z)[2]
    ztrans = kesten(n_z, xi)[2]
    
    # stationary dist
    f = find_statdist(ztrans)
    # CDF
    F = cumsum(f)

    # outside option
    b = zgrid[1] * 0.9999


    # pshape = 2.1
    # pscale = 1.0
    # pdist = Pareto(pshape, pscale)
    # F = (x -> cdf(pdist,x))
    # f = (x -> pdf(pdist, x))

    V0 = u(b)/rho

    q = (x -> x)

    phi = (x -> x)

    V = (x -> x)

    unemp = delta / (delta + lambda)
    # L = (x -> delta * F(x)/(delta + lambda*(1 - F(x))))
    # l = (x -> delta * (delta + lambda) / ((delta + lambda * (1 - F(x)))^2) * f(x))

    # G = (x -> x)
    # Gtilde = ((q,p) -> (delta + lambda*(1-F(p)))^2 / (delta + lambda*(1 - F(q)))^2)

    # Surplus
    S = ((rho + delta)*I + chi*(I - ztrans)) \ (zgrid .- b)

    # MCMC parameters
    # Seed
    seed = 1234
    # workers
    # nI = 1000
    nI = 10000
    # firms
    nJ = 100
    # time periods
    nT = 1000

    # Firm TFP
    zjt_ind = zeros(Int64, nJ,nT)
    zjt_val = zeros(nJ,nT)

    # Worker-firm pairs
    jit = zeros(Int64, nI,nT)
    zit_val = zeros(nI,nT)
    # worker share of Surplus
    betait = zeros(nI, nT)
    # wages
    wit = zeros(nI,nT)
    # wage rank
    wit_q = zeros(nI,nT)

    # time step size
    dt = 1.0
    # prob approximations 
    # offer arrival
    lambda_dt = 1 - exp(-lambda * dt)
    # separation
    delta_dt = 1 - exp(-delta * dt)
    # prod shock
    chi_dt = 1 - exp(-chi * dt)
end
##

# Rouwenhorst discretization of AR(1)
function rouwenhorst(n, mu, sigma, rho)
    # Get width of grid
    width = sqrt((n-1)* sigma^2/(1-rho^2)) # This is correct, though both Greg and Kopecky have typos in their note
    
    # Equi-space grid
    grid = LinRange(mu - width, mu + width, n)

    # Initialize for 2-case
    p = (1 + rho)/2

    if n == 1
        error("nz = 1 not allowed")
    else
        trans = [p 1-p; 1-p p]
        # Recursively find n-case
        if n > 2
            for i in 2:n-1
                # Per Rouwenhorst
                trans = p.*[trans zeros(i); zeros(i)' 0] + (1-p).*[zeros(i) trans; 0 zeros(i)'] + (1-p).*[zeros(i)' 0; trans zeros(i)] + p .* [0 zeros(i)'; zeros(i) trans]
                # Row sums should be 2 for {2,...,n-1}, and 1 for {1,n}
                trans ./= sum(trans, dims=2)
            end
        end
    end

    return grid, trans
end

# function kesten(n,mu,sigma)
function kesten(n, xi)

    # grid = exp.(LinRange(0,n*sigma,n))
    
    sigma = 1
    mu = -xi/2

    trans = zeros(n,n)
    for i in 1:n
        for j in 1:n
            if i == j && j == 1
                trans[i,j] += max(-mu,0)
                trans[i,j] += sigma/2
            elseif i == n && j== n
                trans[i,j] += max(mu,0)
                trans[i,j] += sigma/2
            elseif j == (i - 1)
                trans[i,j] = max(-mu,0)
                trans[i,j] += sigma/2
            elseif j == (i+1)
                trans[i,j] = max(mu,0)
                trans[i,j] += sigma/2
            end
        end
    end

    # trans -= I * (abs(mu) + sigma)

    # trans ./= -(abs(mu) + sigma)

    # trans += I

    trans ./= (abs(mu) + sigma)

    grid = LinRange(0,1,n)
    mean_raw = find_statdist(trans)'* grid
    # xi = -2*mu/(sigma^2)
    # grid *= (xi /(xi - 1)) / mean_raw
    statd = find_statdist(trans)
    function mean_err(x)
        g = exp.(exp(x) * grid)
        return statd' * g - xi/(xi - 1)
    end
    x = find_zero(mean_err, 0.0)
    
    grid = exp.(exp(x) * grid)

    return grid, trans
end

##

# Find stationary distribution of Markov chain with trans mat A
function find_statdist(A)

    maxit = 1000
    tol = 1e-9

    it = 0
    err = Inf
    dist = ones(size(A,1))'
    dist ./= sum(dist)
    while it < maxit && err > tol
        newdist = dist * A
        newdist ./= sum(newdist)

        err = maximum(abs.(newdist - dist))
        it += 1
        dist = newdist
    end

    # println("it = ", it, ", err = ", err)

    return vec(dist)
end

function wage(m, beta, zi)
    @unpack_seq_auction_firm_dyn m
    flow_val = b + rho * beta * S[zi]
    
    sep_val = - delta * beta * S[zi]
    
    raise_val = lambda * (max.(beta * S[zi], S) .- beta * S[zi])[1:zi]' * f[1:zi]
    
    poach_val = lambda * (1 - F[zi]) * (1 - beta) * S[zi]
    
    shock_val = chi * ztrans[zi, :]' * (beta * S .- beta * S[zi])
    
    # println("flow_val: ", flow_val)
    # println("sep_val: ", sep_val)
    # println("raise_val: ", raise_val)
    # println("poach_val: ", poach_val)
    # println("shock_val: ", shock_val)

    return flow_val - sep_val - raise_val - poach_val - shock_val
end

# Markov Chain Monte Carlo simulation
function mcmc!(m)
    @unpack_seq_auction_firm_dyn m

    Random.seed!(seed)

    # Firm TFP paths
    for j in 1:nJ 
        # println("sim j = ", j)

        # Initially stationdary dist
        zjt_ind[j,1] = 1 + sum(rand() .> F)
        zjt_val[j,1] = zgrid[zjt_ind[j,1]]

        for t in 2:nT
            if rand() < chi_dt
                # shock arrived
                zjt_ind[j,t] = 1 + sum(rand() .> cumsum(ztrans[zjt_ind[j,t-1], :]))
            else
                # no shock
                zjt_ind[j,t] = zjt_ind[j,t-1]
            end

            zjt_val[j,t] = zgrid[zjt_ind[j,t]]
        end
    end

    # println("firm sim done")

    # Worker employment paths
    for i in 1:nI
        # println("sim i = ", i)

        # Start everyone unemp, will do burn-in
        # Easier than trying to compute stat dist and assign to firms
        jit[i,1] = 0
        betait[i,1] = 0
        wit[i,1] = b
        zit_val[i,1] = b

        for t in 2:nT
            r = rand()

            shock_occur = false
            if jit[i,t-1] != 0
                if zjt_ind[jit[i,t-1],t-1] != zjt_ind[jit[i,t-1],t]
                    shock_occur = true
                end
            end

            if r < delta_dt
                # separation
                jit[i,t] = 0
                betait[i,t] = 0
                wit[i,t] = b
            # elseif shock_occur
            #     # shock occurred
            #     jit[i,t] = jit[i,t-1]
            #     betait[i,t] = betait[i,t-1]
            #     wit[i,t] = wage(m,betait[i,t],zjt_ind[jit[i,t],t])
            elseif shock_occur && r < (delta_dt + lambda_dt)
                # shock and offer
                j_ind = Int(ceil.(rand() * nJ))
                o_ind = zjt_ind[j_ind,t]

                if o_ind > zjt_ind[jit[i,t-1],t]
                    # new offer better -> accept
                    jit[i,t] = j_ind
                    betait[i,t] = S[zjt_ind[jit[i,t-1],t]] / S[o_ind]
                    wit[i,t] = wage(m, betait[i,t],zjt_ind[jit[i,t],t])
                else
                    # new offer worse -> reject (but could generate raise)
                    jit[i,t] = jit[i,t-1]
                    betait[i,t] = max(betait[i,t-1], S[o_ind] / S[zjt_ind[jit[i,t-1],t]])
                    wit[i,t] = wage(m,betait[i,t],zjt_ind[jit[i,t],t])
                end
            elseif shock_occur 
                # just shock
                jit[i,t] = jit[i,t-1]
                betait[i,t] = betait[i,t-1]
                wit[i,t] = wage(m,betait[i,t],zjt_ind[jit[i,t],t])
            elseif r < (delta_dt + lambda_dt)
                # just offer
                j_ind = Int(ceil.(rand() * nJ))
                o_ind = zjt_ind[j_ind,t]

                if jit[i,t-1] == 0
                    # unemp -> accept
                    jit[i,t] = j_ind
                    betait[i,t] = 0
                    wit[i,t] = wage(m,betait[i,t],zjt_ind[jit[i,t],t])
                elseif o_ind > zjt_ind[jit[i,t-1],t]
                    # new offer better -> accept
                    jit[i,t] = j_ind
                    betait[i,t] = S[zjt_ind[jit[i,t-1],t]] / S[o_ind]
                    wit[i,t] = wage(m,betait[i,t],zjt_ind[jit[i,t],t])
                else
                    # new offer worse -> reject (but could generate raise)
                    jit[i,t] = jit[i,t-1]
                    betait[i,t] = max(betait[i,t-1], S[o_ind] / S[zjt_ind[jit[i,t-1],t]])
                    wit[i,t] = wage(m,betait[i,t],zjt_ind[jit[i,t],t])
                end
            else
                # no change
                jit[i,t] = jit[i,t-1]
                betait[i,t] = betait[i,t-1]
                wit[i,t] = wit[i,t-1]
            end

            if jit[i,t] == 0
                zit_val[i,t] = b
            else
                zit_val[i,t] = zgrid[zjt_ind[jit[i,t], t]]
            end
        end
    end
    # println("worker sim done")

    # Assign to wage quantiles
    for t in 1:nT
        wit_q[:,t] = denserank(wit[:,t])
        wit_q[:,t] ./= maximum(wit_q[:,t])
    end


    # println("m test: ", m.zjt_ind)
    
    @pack_seq_auction_firm_dyn! m
end

function quartile_trans(m)
    @unpack_seq_auction_firm_dyn m

    # Burn in, since starting from unemp and empty firms
    burn_t = ceil(nT*(2/3))
    burn_t = Int(burn_t)

    zqt = zeros(nJ, nT - burn_t)
    for t in 1:(nT - burn_t)
        zqt[:,t] = denserank(zjt_ind[:,t])
        zqt[:,t] ./= maximum(zqt[:,t])
        zqt[:,t] .= ceil.(zqt[:,t] .* 4)
    end
    zqt = Int.(zqt)

    trans_zq = zeros(4,4)
    for j in 1:nJ
        for t in 2:(nT - burn_t)
            trans_zq[zqt[j,t-1], zqt[j,t]] += 1
            # trans_zq[m0.zjt_ind[j,t-1], m0.zjt_ind[j,t]] += 1
        end
    end

    trans_zq ./= sum(trans_zq, dims=2)
    # println("trans_zq:")
    # display(trans_zq)

    # Want to exclude unemp from calc
    wqt = zeros(nI, nT - burn_t)
    for t in 1:(nT - burn_t)
        nob_inds = wit[:,t] .!= b
        # nob_inds = 1:nI
        if sum(nob_inds) > 0
            wqt[nob_inds,t] = denserank(wit[nob_inds,t])
            wqt[nob_inds,t] ./= maximum(wqt[nob_inds,t])
            wqt[nob_inds,t] .= ceil.(wqt[nob_inds,t] .* 4)
        end
    end
    wqt = Int.(wqt)

    trans_wq = zeros(4,4)
    for i in 1:nI
        for t in 2:(nT - burn_t)
            if wqt[i,t-1] > 0 && wqt[i,t] > 0
                trans_wq[wqt[i,t-1], wqt[i,t]] += 1
            end
        end
    end

    trans_wq ./= sum(trans_wq, dims=2)
    # println("trans_wq:")
    # display(trans_wq)

    return trans_zq, trans_wq
end

function plot_firm_size(m; add_to_plot=false)
    @unpack_seq_auction_firm_dyn m

    j_size= zeros(nJ, nT)
    for j in 1:nJ
        # println("size j = ", j)
        j_size[j,:] = sum(jit .== j, dims=1) #./ nI
    end

    kde_size = kde(vec(log.(j_size[j_size .!= 0])), bandwidth = 0.5)
    # kde_size = kde(vec(log.(j_size[j_size .!= 0])))
    # kde_size = kde(vec(j_size[j_size .!= 0]), bandwidth = 100.0)

    if add_to_plot
        plt = plot!(kde_size.x, kde_size.density, 
        label = L"\lambda="*string(lambda)*", "*L"\chi="*string(chi))
    else
        plt = plot(kde_size.x, kde_size.density, xlabel = "n",
        label = L"\lambda="*string(lambda)*", "*L"\chi="*string(chi))
    end
    display(plt)
    return plt
end

##

function job_trans(m)
    @unpack_seq_auction_firm_dyn m

    uu = 0
    ue = 0
    eu = 0
    ee_stay = 0
    ee_move = 0
    for i in 1:nI
        for t in 2:nT
            if jit[i,t-1] == 0
                if jit[i,t] == 0
                    uu += 1
                else
                    ue += 1
                end
            elseif jit[i,t] == 0
                eu += 1
            else
                if jit[i,t-1] == jit[i,t]
                    ee_stay += 1
                else
                    ee_move += 1
                end
            end
        end
    end

    uu /= nI*nT
    ue /= nI*nT
    eu /= nI*nT
    ee_stay /= nI*nT
    ee_move /= nI*nT

    return uu, ue, eu, ee_stay, ee_move
end

function earn_ineq(m)
    @unpack_seq_auction_firm_dyn m

    std_earn = std(log.(wit[wit .!= b]))
    # kde_1 = kde(log.(wit[wit .!= b]))
    # plt = plot(kde_1.x, kde_1.density, label = "earn_1")
    # display(plt)

    # Build 5-year avg
    wit_5 = zeros(nI,nT-4)
    for i in 1:nI
        for t in 1:(nT-4)
            tmid = t + 2
            denom = 0
            for toff in -2:2
                if wit[i,tmid + toff] != b
                    denom += 1
                    wit_5[i,t] += wit[i,tmid + toff]
                end
            end
            if denom > 0
                wit_5[i,t] /= denom
            else
                wit_5[i,t] = b
            end
        end
    end

    std_earn_5 = std(log.(wit_5[wit_5 .!= b]))
    # kde_5 = kde(log.(wit_5[wit_5 .!= b]))
    # plt = plot(kde_5.x, kde_5.density, label = "earn_5")
    # display(plt)

    return std_earn, std_earn_5
end

##

function rank_ar1(m)
    @unpack_seq_auction_firm_dyn m

    zj_rank = zeros(nJ,nT)
    for t in 1:nT
        zj_rank[:,t] .= denserank(zjt_ind[:,t])
        zj_rank[:,t] ./= maximum(zj_rank[:,t])
    end

    lp_rank = zeros(nJ*(nT-1))
    lp_rank_lag = zeros(nJ*(nT-1))
    for j in 1:nJ
        for t in 2:nT
            lp_rank[(j-1)*(nT-1) + t - 1] = zj_rank[j,t]
            lp_rank_lag[(j-1)*(nT-1) + t - 1] = zj_rank[j,t-1]
        end
    end
    # println("lp_rank: ", lp_rank)

    df = DataFrame((lp_rank = lp_rank, lp_rank_lag = lp_rank_lag))
    r_lp = reg(df, @formula(lp_rank ~ lp_rank_lag))


    w_rank = zeros(nI,nT)
    for t in 1:nT
        w_rank[:,t] .= denserank(wit[:,t])
        w_rank[:,t] ./= maximum(w_rank[:,t])
    end

    earn_rank = zeros(nI*(nT-1))
    earn_rank_lag = zeros(nI*(nT-1))
    for i in 1:nI
        for t in 2:nT
            earn_rank[(i-1)*(nT-1) + t - 1] = w_rank[i,t]
            earn_rank_lag[(i-1)*(nT-1) + t - 1] = w_rank[i,t-1]
        end
    end
    
    df = DataFrame((earn_rank = earn_rank, earn_rank_lag = earn_rank_lag))
    r_earn = reg(df, @formula(earn_rank ~ earn_rank_lag))

    return r_lp, r_earn
end

function binscatter_wt(m; add_to_plot=false)
    @unpack_seq_auction_firm_dyn m

    temp_wt_1 = zeros(nI, nT-1)
    temp_wt = zeros(nI, nT-1)

    wt_1 = vec(wit[:,1:(end-1)])
    wt = vec(wit[:,2:end])

    for t in 1:(nT-1)
        temp_wt_1[:,t] .= denserank(wit[:,t])
        temp_wt_1[:,t] ./= maximum(temp_wt_1[:,t])

        temp_wt[:,t] .= denserank(wit[:,t+1])
        temp_wt[:,t] ./= maximum(temp_wt[:,t])
    end

    wt_1 = vec(temp_wt_1[:,1:(end-1)])
    wt = vec(temp_wt[:,2:end])

    # zqt[:,t] = denserank(m0.zjt_ind[:,t])
    # zqt[:,t] ./= maximum(zqt[:,t])

    nq = 100
    # line45 = LinRange(minimum(wt), maximum(wt), nq)
    line45 = LinRange(1/nq, 1, nq)
    wqt_1 = zeros(nq)
    wqt = zeros(nq)
    for iq in 1:nq
        if iq == 1
            wqt_1[iq] = mean(wt_1[wt_1 .<= line45[iq]])
            wqt[iq] = mean(wt[wt_1 .<= line45[iq]])
        else
            wqt_1[iq] = mean(wt_1[wt_1 .> line45[iq-1] .&& wt_1 .<= line45[iq]])
            wqt[iq] = mean(wt[wt_1 .> line45[iq-1] .&& wt_1 .<= line45[iq]])
        end
    end

    if add_to_plot
        plt = scatter!(wqt_1, wqt)
        # plt = plot!(wqt_1, wqt)
    else
        plt = scatter(wqt_1, wqt)
        # plt = plot(wqt_1, wqt)
        plot!(wqt_1, wqt_1, label = "45")
    end
    display(plt)

    return wqt_1, wqt, plt
end

# wqt_10, wqt0, plt0 = binscatter_wt(m0)
# wqt_1lambda, wqtlambda, pltlambda = binscatter_wt(mlambda; add_to_plot=true)
# wqt_1chi, wqtchi, pltchi = binscatter_wt(mchi; add_to_plot=true)

# scatter(wqt_10, wqtlambda - wqt0)
# scatter!(wqt_10, wqtchi - wqt0)
# plt = plot(wqt_10, wqt0 - wqt_10)
# plot!(wqt_1lambda, wqtlambda - wqt_1lambda)
# plot!(wqt_1chi, wqtchi - wqt_1chi)
# display(plt)

function switch_mean_wage_change(m; add_to_plot=false)
    @unpack_seq_auction_firm_dyn m
    
    dws = ones(nI,nT) * -1e6
    j2j = 0
    for i in 1:nI
        for t in 2:nT
            if (jit[i,t-1] != jit[i,t]) && (jit[i,t-1] > 0) && (jit[i,t] > 0)
                dws[i,t] = wit_q[i,t] - wit_q[i,t-1]
                j2j += 1
            end
        end
    end

    dws = dws[dws .!= -1e6]
    
    if isempty(dws)
        dws = [1e6]
    end

    # if add_to_plot
    #     plt = density!(dws, label = "switch wage_q change")
    # else
    #     plt = density(dws, label = "switch wage_q change")
    # end
    # display(plt)

    # return sum(dws)/j2j
    return mean(dws)
end

# switch_mean_wage_change(m0)
# switch_mean_wage_change(mlambda; add_to_plot=true)
# switch_mean_wage_change(mchi; add_to_plot=true)

##

function stay_mean_wage_change(m; add_to_plot=false)
    @unpack_seq_auction_firm_dyn m
    
    dws = ones(nI,nT) * -1e6
    stay = 0
    for i in 1:nI
        for t in 2:nT
            if (jit[i,t-1] == jit[i,t]) && (jit[i,t-1] > 0) && (jit[i,t] > 0)
                dws[i,t] = wit_q[i,t] - wit_q[i,t-1]
                stay += 1
            end
        end
    end

    dws = dws[dws .!= -1e6]

    # if add_to_plot
    #     plt = density!(dws, label = "stay wage_q change")
    # else
    #     plt = density(dws, label = "stay wage_q change")
    # end
    # display(plt)

    # return sum(dws)/stay
    return mean(dws)
end

function stay_share(m)
    @unpack_seq_auction_firm_dyn m
    
    stay = 0
    # kills divide by zero errors
    emp = 1
    for i in 1:nI
        for t in 2:nT
            if (jit[i,t-1] > 0) && (jit[i,t] > 0)
                emp += 1
                if (jit[i,t-1] == jit[i,t])
                    stay += 1
                end
            end
        end
    end
    println("stay = ", stay, ", emp = ", emp)
    return stay / emp
end

function switch_share(m)
    @unpack_seq_auction_firm_dyn m
    
    switch = 0
    # kills divide by zero errors
    emp = 1
    for i in 1:nI
        for t in 2:nT
            if (jit[i,t-1] > 0) && (jit[i,t] > 0)
                emp += 1
                if (jit[i,t-1] != jit[i,t])
                    switch += 1
                end
            end
        end
    end
    println("switch = ", switch, ", emp = ", emp)
    return switch / emp
end

function ue_share(m)
    @unpack_seq_auction_firm_dyn m
    
    ue = 0
    # kills divide by zero errors
    emp = 1
    for i in 1:nI
        for t in 2:nT
            if (jit[i,t] > 0)
                emp += 1
                if (jit[i,t-1] == 0)
                    ue += 1
                end
            end
        end
    end
    println("ue = ", ue, ", emp = ", emp)
    return ue / emp

end

function unemp(m)
    @unpack_seq_auction_firm_dyn m
    
    unemp = 0
    # kills divide by zero errors
    for i in 1:nI
        for t in 1:nT
            if (jit[i,t] == 0)
                unemp += 1
            end
        end
    end
    println("unemp = ", unemp, ", tot = ", nI*nT)
    return unemp / (nI*nT)
end

# j2j_mean_wage_change(m0)
# j2j_mean_wage_change(mlambda)
# j2j_mean_wage_change(mchi)
# stay_mean_wage_change(m0)
# stay_mean_wage_change(mlambda; add_to_plot=true)
# stay_mean_wage_change(mchi; add_to_plot=true)
# println("baseline: ", stay_share(m0))
# stay_share(mlambda)
# stay_share(mchi)
# ue_share(m0)
# switch_share(m0)
# unemp(m0)
##
# Crank out a MCMC, first just to see distribution of firm size, J2J, etc.
# m0 = seq_auction_firm_dyn(nJ=1000, nI=1000, nT = 1000,
# chi = 10.0, lambda = 1.0, delta = 0.01, 
# n_z = 5, mu_z = 10, sigma_z = 2.0, 
# rho_z = 0.0,
# )

# m0 = seq_auction_firm_dyn(nJ=10, nI=1000, nT = 1000,
# chi = 0.05, lambda = 0.1, delta = 0.05, 
# n_z = 11, mu_z = -0.01, sigma_z = 0.1, 
# b = 0.9)

# ##
# mcmc(m0)

# ##
# # Check density of firms
# burn_t = round(Int64, m0.nT * 2/3)
# kde_f = kde(vec(m0.zjt_ind[:, burn_t:end]), bandwidth = 0.6)
# plt = plot(kde_f.x, kde_f.density)
# histogram!(vec(m0.zjt_ind[:,burn_t:end]), normalize=:pdf)
# plot!(m0.f)
# display(plt)
# ##

# # Plot size path of each firm
# j_size= zeros(m0.nJ, m0.nT)
# for j in 1:m0.nJ
#     println("size j = ", j)
#     j_size[j,:] = sum(m0.jit .== j, dims=1) ./ m0.nI
# end
# unemp_size = (1 .- sum(j_size, dims=1))'

# # smooth out the paths
# window = 100
# unemp_smooth = unemp_size[1:end-window+1]
# j_smooth = j_size[:, 1:end-window+1]
# for t in 2:window
#     j_smooth .+= j_size[:, t:end-window+t]
#     unemp_smooth .+= unemp_size[t:end-window+t]
# end
# unemp_smooth ./= window
# j_smooth ./= window

# plt = plot(unemp_size, label = "unemp")
# for j in 1:m0.nJ
#     println("plot size j: ", j)
#     plot!(j_size[j,:], label = j)
# end

# display(plt)

# # plt = plot(unemp_smooth, label = "unemp")
# plt = plot()
# for j in 1:m0.nJ
#     println("plot size j: ", j)
#     plot!(j_smooth[j,:], label = j)
# end

# display(plt)

# ##
# # J2J flows
# j2j = zeros(m0.nJ+1, m0.nJ+1)
# for i in 1:m0.nI
#     for t in 2:m0.nT
#         j2j[m0.jit[i,t-1] + 1, m0.jit[i,t] + 1] += 1
#     end
# end
# j2j ./= sum(j2j, dims=2)

# ##
# # Now want to see transitions between z levels
# jit_zind = zeros(m0.n_z + 1, m0.n_z + 1)
# for i in 1:m0.nI
#     for t in 2:m0.nT
#         if m0.jit[i,t-1] == 0 && m0.jit[i,t] == 0
#             jit_zind[1,1] += 1
#         elseif m0.jit[i,t-1] == 0
#             jit_zind[1, m0.zjt_ind[m0.jit[i,t],t] + 1] += 1
#         elseif m0.jit[i,t] == 0
#             jit_zind[m0.zjt_ind[m0.jit[i,t-1],t] + 1, 1] += 1
#         else
#             jit_zind[m0.zjt_ind[m0.jit[i,t-1],t-1] + 1, m0.zjt_ind[m0.jit[i,t],t] + 1] += 1
#         end
#     end
# end

# emp_trans = jit_zind[2:end, 2:end]
# emp_trans ./= sum(emp_trans, dims=2)

# unemp_out = jit_zind[1,2:end]
# unemp_out ./= sum(unemp_out)

# unemp_in = jit_zind[2:end,1]
# unemp_in ./= sum(unemp_in)

# plot(m0.f, label = "f")
# plot!(unemp_out, label = "out")
# plot!(unemp_in, label = "in")

# ##
# # Want to see dist of z across workers
# j_inds = m0.jit .!= 0
# zs = zeros(size(m0.jit))
# for i in 1:m0.nI
#     for t in 1:m0.nT
#         if m0.jit[i,t] == 0
#             zs[i,t] = 0
#         else
#             zs[i,t] = m0.zjt_ind[m0.jit[i,t],t]
#         end
#     end
# end
# histogram(vec(zs), normalize=:pdf)

# ##
# # Checking basic wage patterns
# nb = 1000
# bets = LinRange(0,1,nb)
# ws =[wage(m0, bets[bi], zi) for zi in 1:m0.n_z, bi in 1:nb]
# plt = plot(legend=false)
# for zi in 1:m0.n_z
#     plot!(bets, ws[zi,:], label = "zi = "*string(zi))
# end
# display(plt)

# ##
# # wage Distributions (account for burn-in!)
# burn_t = 700
# kde_w = kde(vec(m0.wit[:,burn_t:1000]))
# # kde_w = kde(log.(vec(m0.wit[:,burn_t:1000])))
# # kde_w = kde(vec(m0.wit[:,burn_t:1000]), bandwidth = 0.1)
# # kde_w = kde(vec(m0.wit[500,burn_t:end]))
# plot(kde_w.x, kde_w.density)

# ##
# # Now want to think about wage distributions within firms
# j = 50
# for t in 1:10:m0.nT
#     jinds = m0.jit[:,t] .== j
#     ws = m0.wit[jinds,t]
#     println("t: ", t, ", length(ws): ", length(ws))
#     if length(ws) > 0
#         kde_w = kde(ws)
#         plt = plot(kde_w.x, kde_w.density)
#         display(plt)
#     end
# end


# ##
# m0 = seq_auction_firm_dyn(nJ=100, nI=1000, nT = 10000,
# chi = 0.3, lambda = 0.1, delta = 0.02, 
# n_z = 11, mu_z = -0.05, sigma_z = sqrt(0.05), 
# b = 0.95,
# dt = 1.0)

# ##
# mcmc(m0)


# ##
# # Get earnings quartile and LP quartile transition matrices
# # How do they both change when messing with frictions or LP process?
# burn_t = 0
# zqt = zeros(m0.nJ, m0.nT - burn_t)
# for t in 1:(m0.nT - burn_t)
#     zqt[:,t] = denserank(m0.zjt_ind[:,t])
#     zqt[:,t] ./= maximum(zqt[:,t])
#     zqt[:,t] .= ceil.(zqt[:,t] .* 4)
#     # zqt[:,t] = ceil.((denserank(m0.zjt_ind[:,t]) ./ m0.nJ) .* 4)
# end
# zqt = Int.(zqt)

# trans_zq = zeros(4,4)
# for j in 1:m0.nJ
#     for t in 2:(m0.nT - burn_t)
#         trans_zq[zqt[j,t-1], zqt[j,t]] += 1
#         # trans_zq[m0.zjt_ind[j,t-1], m0.zjt_ind[j,t]] += 1
#     end
# end

# trans_zq ./= sum(trans_zq, dims=2)
# println("trans_zq:")
# display(trans_zq)

# # Want to exclude unemp from calc
# wqt = zeros(m0.nI, m0.nT - burn_t)
# for t in 1:(m0.nT - burn_t)
#     # nob_inds = m0.wit[:,t] .!= m0.b
#     nob_inds = 1:m0.nI
#     if sum(nob_inds) > 0
#         wqt[nob_inds,t] = denserank(m0.wit[nob_inds,t])
#         wqt[nob_inds,t] ./= maximum(wqt[nob_inds,t])
#         wqt[nob_inds,t] .= ceil.(wqt[nob_inds,t] .* 4)
#     end
# end
# wqt = Int.(wqt)

# trans_wq = zeros(4,4)
# for i in 1:m0.nI
#     for t in 2:(m0.nT - burn_t)
#         trans_wq[wqt[i,t-1], wqt[i,t]] += 1
#     end
# end

# trans_wq ./= sum(trans_wq, dims=2)
# println("trans_wq:")
# display(trans_wq)

# ##
# # Baseline
# # m0 = seq_auction_firm_dyn(nJ=100, nI=100000, nT = 100,
# # chi = 0.1, lambda = 0.5, delta = 0.02, 
# # n_z = 11, mu_z = -0.05, sigma_z = 0.05, 
# # b = 0.95,
# # dt = 0.1)
# m0 = seq_auction_firm_dyn(nJ=1000, nI=10000, nT = 1000,
# chi = 0.5, lambda = 0.2, delta = 0.05, 
# n_z = 100, mu_z = -0.05, sigma_z = sqrt(0.05), 
# b = 0.95,
# dt = 1.0)
# mcmc(m0)

# # Decrease lambda
# mlambda = seq_auction_firm_dyn(nJ=1000, nI=10000, nT = 1000,
# chi = 0.5, lambda = 0.1, delta = 0.05, 
# n_z = 100, mu_z = -0.05, sigma_z = sqrt(0.05), 
# b = 0.95,
# dt = 1.0)
# mcmc(mlambda)


# # # Decrease chi
# mchi = seq_auction_firm_dyn(nJ=1000, nI=10000, nT = 1000,
# chi = 0.25, lambda = 0.2, delta = 0.05, 
# n_z = 100, mu_z = -0.05, sigma_z = sqrt(0.05), 
# b = 0.95,
# dt = 1.0)
# mcmc(mchi)


# ##
# # See quartile transition
# z0_trans, w0_trans = quartile_trans(m0)
# println("baseline z_trans: ")
# display(z0_trans)
# println()
# println("baseline w_trans: ")
# display(w0_trans)

# zlambda_trans, wlambda_trans = quartile_trans(mlambda)
# println("lambda z_trans: ")
# display(zlambda_trans)
# println()
# println("lambda w_trans: ")
# display(wlambda_trans)

# zchi_trans, wchi_trans = quartile_trans(mchi)
# println("chi z_trans: ")
# display(zchi_trans)
# println()
# println("chi w_trans: ")
# display(wchi_trans)


# # quartile_tran(mlambda)
# # quartile_trans(mchi)

# ##
# # Firm size distribution
# plot_firm_size(m0)
# plot_firm_size(mlambda; add_to_plot=true)
# plt = plot_firm_size(mchi; add_to_plot=true)
# # savefig(plt, "firm_size.svg")
# ##
# # J2J stayer vs. mover
# uu, ue, eu, ee_stay, ee_move = job_trans(m0)
# println("baseline J2J: ", ee_move / m0.dt)
# println("baseline J2J (normalized): ", ee_move/(eu + ee_stay + ee_move) / m0.dt)
# uu, ue, eu, ee_stay, ee_move = job_trans(mlambda)
# println("lambda J2J: ", ee_move / m0.dt)
# println("lambda J2J (normalized): ", ee_move/(eu + ee_stay + ee_move) / m0.dt)
# uu, ue, eu, ee_stay, ee_move = job_trans(mchi)
# println("chi J2J: ", ee_move / m0.dt)
# println("chi J2J (normalized): ", ee_move/(eu + ee_stay + ee_move) / m0.dt)

# ##
# # cross-section vs 5-year std earn
# std_1, std_5 = earn_ineq(m0)
# println("baseline: std = ", std_1, ", std_5 = ", std_5)
# std_1, std_5 = earn_ineq(mlambda)
# println("lambda: std = ", std_1, ", std_5 = ", std_5)
# std_1, std_5 = earn_ineq(mchi)
# println("chi: std = ", std_1, ", std_5 = ", std_5)

# ##
# # rank ar(1)
# r_lp, r_earn = rank_ar1(m0)
# println("baseline AR(1)s: ")
# println(r_lp)
# println()
# println(r_earn)

# r_lp, r_earn = rank_ar1(mlambda)
# println("lambda AR(1)s: ")
# println(r_lp)
# println()
# println(r_earn)

# r_lp, r_earn = rank_ar1(mchi)
# println("chi AR(1)s: ")
# println(r_lp)
# println()
# println(r_earn)
# ##
# # Mean wage change J2J
# dw = j2j_mean_wage_change(m0)
# println("baseline dw (J2J): ", dw)
# dw = j2j_mean_wage_change(mlambda; add_to_plot=true)
# println("lambda dw (J2J): ", dw)
# dw = j2j_mean_wage_change(mchi; add_to_plot=true)
# println("chi dw (J2J): ", dw)
# println()
# # stayer change
# dw = stay_mean_wage_change(m0)
# println("baseline dw (stayer): ", dw)
# dw = stay_mean_wage_change(mlambda; add_to_plot=true)
# println("lambda dw (stayer): ", dw)
# dw = stay_mean_wage_change(mchi; add_to_plot=true)
# println("chi dw (stayer): ", dw)

# ##
# println("baseline stay share: ", stay_share(m0))
# println("lambda stay share: ", stay_share(mlambda))
# println("chi stay share: ", stay_share(mchi))


# ##
# # Wage distribution
# # density(vec(log.(m0.wit[m0.wit .!= m0.b])))
# burn_t = 900
# w_burn = m0.wit[:, burn_t:end]
# w_burn = w_burn[w_burn .!= m0.b]
# kde_w = kde(vec(log.(w_burn)), bandwidth = 0.05)
# plot(kde_w.x, kde_w.density, label = "baseline")

# w_burn = mlambda.wit[:, burn_t:end]
# w_burn = w_burn[w_burn .!= mlambda.b]
# kde_w = kde(vec(log.(w_burn)), bandwidth = 0.05)
# plot!(kde_w.x, kde_w.density, label = "lambda")

# w_burn = mchi.wit[:, burn_t:end]
# w_burn = w_burn[w_burn .!= mchi.b]
# kde_w = kde(vec(log.(w_burn)), bandwidth = 0.05)
# plot!(kde_w.x, kde_w.density, label = "chi")

# ##
# function plot_wage_dist(m; name="unlabelled")
#     burn_t = Int(floor(m.nT * 2/3))
#     w_burn = m.wit[:, burn_t:end]
#     w_burn = w_burn[w_burn .!= m.b]
#     kde_w = kde(vec(log.(w_burn)), bandwidth = 0.05)
#     return plot(kde_w.x, kde_w.density, label = name)
# end

# plt = plot_wage_dist(mcal)
# display(plt)
# ##
# # Joint distribution over wages and LP
# burn_t = 900
# w_burn = m0.wit[:, burn_t:end]
# w_burn = w_burn[w_burn .!= m0.b]
# ws = vec(log.(w_burn))
# z_burn = m0.zit_val[:, burn_t:end]
# z_burn = z_burn[z_burn .!= m0.b]
# zs = vec(z_burn)

# kzw = kde((zs, ws), bandwidth = (0.01, 0.005))
# plt = surface(kzw.y, kzw.x, kzw.density, camera = (20,40), 
# color=reverse(cgrad(:RdYlBu_11)), colorbar= :false)
# display(plt)
# plt = heatmap(kzw.x, kzw.y, kzw.density', color=reverse(cgrad(:RdYlBu_11)))
# display(plt)

# # distribution over z (firm-level)
# # kz = kde(log.(vec(m0.zjt_val)), bandwidth = 0.01)
# # plt = plot(kz.x, kz.density, label = "firms")

# # # distribution over z (worker-level)
# # kz = kde(log.(zs), bandwidth = 0.01)
# # plot!(kz.x, kz.density, label = "workers")

# # display(plt)


# ##
# # Binscatter wage transitions
# wt_1 = vec(m0.wit[:,1:(end-1)])
# wt = vec(m0.wit[:,2:end])
# n = 10000
# inds = rand(1:length(wt),n)

# nq = 10
# q_cuts = zeros(nq)
# for iq in 1:nq
#     q_cuts[iq] = quantile(wt_1, iq/nq)
# end
# wqt_1 = zeros(nq)
# wqt = zeros(nq)
# for iq in 1:nq
#     if iq == 1
#         wqt_1[iq] = mean(wt_1[wt_1 .<= q_cuts[iq]])
#         wqt[iq] = mean(wt[wt_1 .<= q_cuts[iq]])
#     else
#         wqt_1[iq] = mean(wt_1[wt_1 .> q_cuts[iq-1] .&& wt_1 .<= q_cuts[iq]])
#         wqt[iq] = mean(wt[wt_1 .> q_cuts[iq-1] .&& wt_1 .<= q_cuts[iq]])
#     end
# end




# # zqt[:,t] = denserank(zjt_ind[:,t])
# #         zqt[:,t] ./= maximum(zqt[:,t])
# #         zqt[:,t] .= ceil.(zqt[:,t] .* 4)
# # scatter(wt_1[inds], wt[inds])

# # plot!(wt_1[inds], wt_1[inds])
# scatter(wqt_1, wqt)
# line45 = LinRange(minimum(wt), maximum(wt), 100)
# plot!(line45, line45)

# ##
# # Alternatively, just do equal-sized binning?
# nq = 10
# line45 = LinRange(minimum(wt), maximum(wt), nq)
# wqt_1 = zeros(nq)
# wqt = zeros(nq)
# for iq in 1:nq
#     if iq == 1
#         wqt_1[iq] = mean(wt_1[wt_1 .<= line45[iq]])
#         wqt[iq] = mean(wt[wt_1 .<= line45[iq]])
#     else
#         wqt_1[iq] = mean(wt_1[wt_1 .> line45[iq-1] .&& wt_1 .<= line45[iq]])
#         wqt[iq] = mean(wt[wt_1 .> line45[iq-1] .&& wt_1 .<= line45[iq]])
#     end
# end

# scatter(wqt_1, wqt)
# plot!(line45, line45)

# ##
# # Checking my understanding on quantile rank regression
# n = 10000
# x = LinRange(0,1,n)
# y = circshift(x,500)
# # y = shuffle(x)

# df = DataFrame((x = x, y = y))

# r = reg(df, @formula(y ~ x))
# display(r)

# mean(x - y)

##
# Start of calibration approach (will want to do full in separate file)

function mom_err(x, moms)
    # Flow to unemployment
    # e_to_u = moms[1]
    # stayer share
    # stayer = moms[2]

    # u_to_e = moms[2]

    # switch = moms[2]

    # earn_p = moms[2]

    # xi = moms[3]

    # Drift
    chi = moms[1]
    # Pareto
    xi = moms[2]
    # e_to_u
    delta = moms[3]
    # stay mean wage change
    # stay_dwq = moms[4]
    # or switch mean wage change
    switch_dwq = moms[4]

    # Guess parameters
    # delta = exp(x[1])
    # lambda = exp(x[2])

    lambda = exp(x)
    # Initialize and solve model
    m = seq_auction_firm_dyn(lambda = lambda, delta = delta,
    chi = chi, xi = xi)

    mcmc!(m)

    # compute moments
    # err = Inf* ones(2)
    # # flow to unemp
    # err[1] = m.delta_dt - e_to_u
    # # err[2] = stay_share(m) - stayer
    # err[2] = switch_share(m) - switch
    # err[2] = m.lambda_dt - u_to_e
    # err[2] = rank_ar1(m)[2].coef[2] - earn_p
    # err[2] = 0.0
    # err = [switch_share(m) - switch]
    # println("err: ", err)
    # # println("e_to_u model: ", m.delta_dt, ", moment: ", e_to_u)
    # # println("stay model: ", stay_share(m), ", moment: ", stayer)
    # # println("switch model: ", switch_share(m), ", moment: ", switch)
    # println("lambda: ", m.lambda, ", delta: ", m.delta)
    # println("lambda_dt: ", m.lambda_dt, ", delta_dt: ", m.delta_dt)

    # err = stay_dwq - stay_mean_wage_change(m)
    err = switch_dwq - switch_mean_wage_change(m)
    println("lambda: ", lambda, ", err: ", err)
    # so can see immediate update
    flush(stdout)
    return err
end

##
# lds = LinRange(-10, 10, 10)
# errs = zeros(length(lds))
# for i in 1:length(lds)
#     errs[i] = mom_err([lds[i], -2.0], [0.1, 0.7, 0.05])[1]
# end

# plt = plot(lds, errs)
# display(plt)

##
# lls = LinRange(-20, 10, 20)
# errs = zeros(length(lls))
# for i in 1:length(lls)
#     errs[i] = mom_err([-3.0, lls[i]], [0.1, 0.7, 1.1])[2]
# end

# plt = plot(lls, errs)
# display(plt)

# ##
# m0 = seq_auction_firm_dyn(delta = 0.05, lambda = 0.2, xi = 1.1, n_z = 15)
# mcmc(m0)

##

function calibrate_model(moms)
    # # Flow to unemployment
    # e_to_u = moms[1]
    # # stayer share
    # stayer = moms[2]

    # sig_z = moms[3]
    xi = moms[3]

    # Drift
    chi = moms[1]
    # Pareto
    xi = moms[2]
    # e_to_u
    delta = moms[3]
    # stay mean wage change
    stay_dwq = moms[4]

    # function mom_err(x)
    #     # Guess parameters
    #     delta = exp(x[1])
    #     lambda = exp(x[2])
    #     # Initialize and solve model
    #     m = seq_auction_firm_dyn(delta = delta, lambda = lambda, 
    #     sigma_z = sig_z)
    #     mcmc(m)
    #     # compute moments
    #     err = Inf* ones(2)
    #     # flow to unemp
    #     err[1] = m.delta_dt - e_to_u
    #     err[2] = stay_share(m) - stayer
    #     println("err: ", err)
    #     return err
    # end

    # init_guess = [0.05, 0.2]
    # init_guess = [log(0.05), log(0.2)]
    init_lambda = log(0.2)
    # min error (or zero, if identified)
    # sol = nlsolve(x -> mom_err(x,moms), init_guess,
    # show_trace = true, iterations = 10,
    # # autodiff = :forward,
    # # method = :broyden,
    # # ftol = m.tol,
    # # rtol = m.tol
    # )

    # # delta = exp(sol.zero[1])
    # # lambda = exp(sol.zero[2])
    # lambda = exp(sol.zero[1])
    # m = seq_auction_firm_dyn(lambda = lambda, 
    # xi = xi, n_z = 15, dt = 1.0)

    # sol = find_zero(x -> mom_err(x, moms), log(0.05); xatol=1e-6)
    lambda = exp(find_zero(x -> mom_err(x, moms), init_lambda))

    m = seq_auction_firm_dyn(lambda = lambda, delta = delta,
    chi = chi, xi = xi)

    # For quick debugging
    # m = seq_auction_firm_dyn(lambda = exp(init_lambda), delta = delta,
    # chi = chi, xi = xi)

    # Need to fix so can allow to adapt to chi input

    mcmc!(m)

    return m
end

# mcal = calibrate_model([0.1, 0.02, 1.1])



##
# Need to plot moments as functions of parameters to 
# make sure errors cross zero
# lls = LinRange(-4.5, 1, 100)
# ls = exp.(lls)
# errs = zeros(length(ls))
# j2jcs = zeros(length(ls))
# staycs = zeros(length(ls))
# for i in 1:length(lls)
#     # errs[i] = mom_err([lls[i]], [0.1, 0.02, 1.1])[1]

#     m0 = seq_auction_firm_dyn(delta = 0.05, lambda = ls[i], 
#     chi = 0.1, xi = 1.1, n_z = 15, nI = 100, nT = 100, nJ = 100)
#     mcmc(m0)

#     j2jcs[i] = j2j_mean_wage_change(m0)
#     staycs[i] = stay_mean_wage_change(m0)
# end

# plt = plot(ls, j2jcs, label = "j2j changes")
# display(plt)
# plt = plot(ls, staycs, label = "stay changes")
# display(plt)

# ##
# lcs = LinRange(-4.5, 1, 20)
# cs = exp.(lcs)
# errs = zeros(length(cs))
# j2jcs = zeros(length(cs))
# staycs = zeros(length(cs))
# for i in 1:length(cs)
#     # errs[i] = mom_err([lls[i]], [0.1, 0.02, 1.1])[1]

#     m0 = seq_auction_firm_dyn(delta = 0.05, lambda = 0.2, 
#     chi = cs[i], xi = 1.1, n_z = 15, nI = 1000, nT = 100, nJ = 1000)
#     mcmc(m0)

#     j2jcs[i] = j2j_mean_wage_change(m0)
#     staycs[i] = stay_mean_wage_change(m0)
# end

# plt = plot(cs, j2jcs, label = "j2j changes")
# display(plt)
# plt = plot(cs, staycs, label = "stay changes")
# display(plt)

# ##
# # Example moment to hit (not from any real data):
# # LP drift = -0.2, LP Pareto = 4.0
# # e_to_u (so delta) = 0.15
# # change in earnings rank for stayers 0.003 (try with -0.003 also)
# # change in earnings rank for switchers -0.02 (try with +0.02?)

# # mcal = calibrate_model([0.2, 4.0, 0.15, 0.003])
# # println("0.003 lambda: ", mcal.lambda)
# # mcal = calibrate_model([0.2, 4.0, 0.15, -0.003])
# # println("-0.003 lambda: ", mcal.lambda)

# # mcal = calibrate_model([0.2, 4.0, 0.15, -0.02])
# # println("-0.02 lambda: ", mcal.lambda)
# # mcal = calibrate_model([0.2, 4.0, 0.15, 0.02])
# # println("0.02 lambda: ", mcal.lambda)

# # Want to try some more sample calibrations
# # mcal = calibrate_model([0.1, 3.0, 0.1, -0.004])
# # println("-0.004 lambda: ", mcal.lambda)
# # mcal = calibrate_model([0.1, 10.0, 0.2, -0.02])
# # println("lambda: ", mcal.lambda)
# # mcal = calibrate_model([0.23, 3.0, 0.2, 0.06])
# # println("lambda: ", mcal.lambda)
# # mcal = calibrate_model([0.06, 9.0, 0.15, 0.08])
# # println("lambda: ", mcal.lambda)

# mcal0 = calibrate_model([0.2, 3.0, 0.15, -0.02])
# # # mcallambda = calibrate_model([0.2, 3.0, 0.15, -0.02]) 
# # # mcalchi = calibrate_model([0.15, 3.0, 0.15, -0.02])
# mcal1 = calibrate_model([0.15, 3.0, 0.15, -0.01])

# mcallambda = deepcopy(mcal0)
# mcallambda.lambda = mcal1.lambda
# mcallambda.lambda_dt = 1 - exp(-mcallambda.lambda * mcallambda.dt)
# mcmc(mcallambda)

# mcalchi = deepcopy(mcal0)
# mcalchi.chi = mcal1.chi
# mcalchi.chi_dt = 1 - exp(-mcalchi.chi * mcalchi.dt)
# mcmc(mcalchi)


# mom0 = switch_mean_wage_change(mcal0)
# mom1 = switch_mean_wage_change(mcal1)
# momlambda = switch_mean_wage_change(mcallambda)
# momchi = switch_mean_wage_change(mcalchi)
# println("switch mean wage change")
# println("total: ", mom1 - mom0)
# println()
# println("m1-mlambda: ", mom1 - momlambda)
# println("mlambda-m0: ", momlambda - mom0)
# println()
# println("m1-mchi: ", mom1 - momchi)
# println("mchi-m0: ", momchi - mom0)

# println("-------------------")

# mom0 = stay_mean_wage_change(mcal0)
# mom1 = stay_mean_wage_change(mcal1)
# momlambda = stay_mean_wage_change(mcallambda)
# momchi = stay_mean_wage_change(mcalchi)
# println("stay mean wage change")
# println("total: ", mom1 - mom0)
# println()
# println("m1-mlambda: ", mom1 - momlambda)
# println("mlambda-m0: ", momlambda - mom0)
# println()
# println("m1-mchi: ", mom1 - momchi)
# println("mchi-m0: ", momchi - mom0)







