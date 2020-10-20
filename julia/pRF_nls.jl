# pRF fitting in julia - using LsqFit.jl
#
# a worked example of doing the pRF fitting in julia.
# also trying to establish how fast it is + how easy to extend.
#
# ds 2020-10

# dependencies
using NIfTI, MAT 
using Plots
using Colors
using LinearAlgebra
using ImageFiltering # convolution (more choice in algos)
using StatsBase # mean, etc
using Distributions # Normal, for quick bootstrapping
using LsqFit # curve_fit lmfit

# helper functions
include("prf_nls_helpers.jl")

# set gr backend for plotting
gr() 

# for now: data in parent folder
tseries_file = "../tSeries_10_31_18.mat"
stim_image_file = "../pRFStimImage_example.mat"

# unpack the dict that is being read in:
ts = matread(tseries_file)["tSeries"]
stimImage = matread(stim_image_file)["pRFStimImage"] 

# display an image from stim_image, 3d video is easy
# just one line of code:
img = permutedims(stimImage["im"], (2, 1, 3));
heatmap(img[:,:,39],c=:hot)

# vectors / extent of x, y 
x = stimImage["x"][:,1,1];
y = stimImage["y"][1,:,1];

# mesh versions
mx = stimImage["x"];
my = stimImage["y"];

# Parameters that control the pRF 
# (the "region of visual space that drives the response")
x0 = -2.0;
y0 = +2.0;
sigma_xy = 3.0;
amp = 100.0; # an additional scale factor.

# calculate
p = [x0; y0; sigma_xy; amp];
pRFhat = mgauss(p, mx, my);

# heatmap the pRF
# display - a transpose is needed on the Z image / pRFhat aspect_ratio
# as per plotting conventions for heatmap.
# the following should be equivalent to what happens in mrTools.
hm_ = heatmap(x,y,transpose(pRFhat), 
        c=:hot,
        yflip=true,
        aspect_ratio=1,
        xaxis="y (deg)", yaxis="x (deg)") # use jet1 for matlab look!
hline!([0.0], c=:white, label="")
vline!([0.0], c=:white, label="")
savefig(hm_, "fig-pRFimage.png") # really best to use SVG/PDF

# get a haemodynamic response function
TR = 1.5; # acq TR
t_hrf = 0:TR:30;
hrf = fmribHRF(t_hrf);
# msw - markerstrokewidth, mc - markercolor, ms = markersize
hrf_ = plot(t_hrf, hrf, 
            lw=3, msw=3, m=:circle, ms=8, mc=:white, c=:black, 
            xlabel="Time (s)", ylabel="response (a.u.)",
            label="", 
            title="haemodynamic response")
savefig(hrf_, "fig-hrf-plot.png") # really best to use SVG/PDF

# details of the stim image sequence
t = transpose(stimImage["t"])
szStimIm = size(stimImage["im"])
nPixels = prod(szStimIm[1:2])
nT = szStimIm[3]

# unwrap the stim image into a nx * ny by t matrix.
# only needed once, so saves time by being available in this form.
@time stimImageUW = reshape(stimImage["im"], nPixels, nT)

# pRF response without HRF effect
@time r = calculate_prf_response(stimImageUW, pRFhat);
# try BenchmarkTools.jl and @btime macro (instead of @time) 
# to measure performance numbers a bit more accurately (including warmup)
# on my (slightly tired laptop... ~1.25ms w/o ) 

# ... and with HRF - doesn't take longer (@btime says :)
@time r_hrf = calculate_prf_response(stimImageUW, pRFhat, hrf);

# make a plot of stimulated pRF response
prf1_ = plot(t, r, c=:black, lw=2, label="pRF", title="simple, simulated pRF response")
plot!(t, r_hrf, c=:blue, lw=2, label="pRF w/ HRF", xlabel="Time (s)", ylabel="response")
savefig(prf1_, "fig-prf-example-plot.png") # really best to use SVG/PDF

## now set up the optimisation / fitting
#
# 1. COST FUNCTION
#
# define the cost function for nonlinear least squares
# nb. data, mx and my and hrf are "baked" into this closure
function cost_func(x)
    measured = [ts...;]
    model = x[4] .* calculate_prf_response(stimImageUW, mgauss(x[1:3], mx, my), hrf) 
    # residual = abs.(measured .- model)
    residual = (measured .- model)
    return residual
end


# 2. OPTIMIZE / run NLS
#
# starting parameters and bounds
p0 = [-3.0, -3.0, 1.0, 10000.0]; # starting vals
lb = [-30.0, -30.0, 0.5, 1.0];
ub = [+30.0, +30.0, 8.0, +Inf]

# Levenberg-Marquardt as implemented in LsqFit.
@time res = LsqFit.lmfit(cost_func, p0, Float64[]; 
                         autodiff=:forward,
                         lower=lb,
                         upper=ub)
res.converged ? println("✓ lmfit converged...") :  println("\:x:k lmfit didn't converge...")                       

# 2b. ELABORATED VERSION / using gradient w/ autodiff
#
# there may be a bit more accuracy and speed if we also try to use
# gradient (forwardiff / numerical, but could think about explicitly defining
# jacobian?)
@time resad = LsqFit.lmfit(cost_func, p0, Float64[], 
            autodiff=:forward)
resad.param

# 3. CHECK and PLOT results
#
# the solution parameters are here
res.param
# LsqFit.jl includes fumctions for CI and SE.
ci = confidence_interval(res, 0.01)
se = standard_errors(res)

# this function takes a param vector and turns it into a timeseries prediction
# cf the cost_func... this makes plotting and bootstrapping responses easier.
function calcFit(x) 
    x[4] .* calculate_prf_response(stimImageUW, mgauss(x[1:3], mx, my), hrf) 
end

# this is taking the SE of each parameter and making the 
# extreme prediction. NB! ignores any interactions / covariance in the 
# parameters in the model... but easy to visualise
seMinus = calcFit(res.param .- se) 
sePlus = calcFit(res.param .+ se) 

perror_ = plot(ts, c=:black, lw=2, label="", title="Naive approach to error bars, pHat ± 1 SE ")
plot!(calcFit(res.param), c=:pink, ribbon=(seMinus, sePlus), label="+/- 1 SE")
plot!(calcFit(res.param), c=:red, lw=2, label="NLS fit")
savefig(perror_, "fig-prf-tc+error-plot.png") # really best to use SVG/PDF

# 4. ELABORATED version of error bars / bootstrapping
#
# make a multivariate distribution of params 
d = Distributions.MultivariateNormal(res.param, se)

println("get a 1000 samples of params")
@time d_boostrapped = rand(d, 1000)

println("get a 1000 boostrap replicates of timeseries")
@time boots = mapslices(calcFit, d_boostrapped; dims=1)

# and plot those boostrap replicates
pboot_ = plot(t, boots, la=0.1, c=:gray, label="")
plot!(t, ts, c=:black, lw=2, label="", title="Boostrapped error bars, pHat, 1000 boostraps ")
plot!(t, calcFit(res.param), c=:red, lw=2, label="NLS fit")
savefig(perror_, "fig-prf-tc+boot-plot.png") # really best to use SVG/PDF



## a way to plot pRF shapes (parametric plots)
# r1 = res.param .+ se
# r2 = res.param .- se
# xₜ(t) = r1[3] .* cos(t) .+ r1[1] # , r2[3] .* cos(t) .+ r2[1]) 
# yₜ(t) = r1[3] .* sin(t) .+ r1[2] # , r2[3] .* sin(t) .+ r2[2])

# plot(xₜ, yₜ, 0, 2π, leg=false, fill=(0, :orange),        xlim=(-15, +15), ylim=(-15, +15), aspect_ratio=1)

# https://stackoverflow.com/a/56174228/4961292 

