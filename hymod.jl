using DataFrames
using CSV
using Pkg

using ClearStacktrace

input = "/home/iff/research/dev/nsgaii/data/hymod/hymod_input.csv"

df = DataFrame(CSV.File(input))



names(df) = ["date","precip","TURC","dis"]

df = select(df, :Date => :date, :"rainfall[mm]" => :precip, :"TURC [mm d-1]" => :TURC,:"Discharge[ls-1]" => :dis) # rename columns

df


df = filter(:dis => x -> !(ismissing(x) || isnothing(x) || isnan(x)), df)

cmax  = 412.33
bexp  = 0.1725
alpha = 0.8127
Ks    = 0.0404
Kq    = 0.5592

Precip = df.precip
Pet = df.TURC



function hymod(
    Precip::Array{Float64,1},
    PET::Array{Float64,1},
    cmax::Float64,
    bexp::Float64,
    alpha::Float64,
    Rs::Float64,
    Rq::Float64)

    # HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL

    x_loss = 0.0
    # Initialize slow tank state
    x_slow = 2.3503 / (Rs * 22.5)
    #x_slow = 0  # --> works ok if calibration data starts with low discharge
    # Initialize state(s) of quick tank(s)
    x_quick = zeros(3)
    t = 1
    outflow = []
    output = []

    # START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS
    while t <= length(Precip)


         Pval = Precip[t]
         PETval = PET[t]
    #     # Compute excess precipitation and evaporation
         ER1, ER2, x_loss = excess(x_loss, cmax, bexp, Pval, PETval)
    #     # Calculate total effective rainfall
         ET = ER1 + ER2
    #     #  Now partition ER between quick and slow flow reservoirs
         UQ = alpha * ET
         US = (1 - alpha) * ET
    #     # Route slow flow component with single linear reservoir
         x_slow, QS = linres(x_slow, US, Rs)
    #     # Route quick flow component with linear reservoirs
         inflow = UQ

         for i in 1:1:3
    #         # Linear reservoir

             x_quick[i], outflow = linres(x_quick[i], inflow, Rq)
             inflow = outflow
         end

    #     # Compute total flow for timestep
         output = push!(output,QS + outflow)

         t +=1

    end

    return output


end


function power(X,Y)
    X=abs(X) # Needed to capture invalid overflow with netgative values
    return X^Y
end

function linres(x_slow,inflow,Rs::Float64)
    # Linear reservoir
    x_slow = (1 - Rs) * x_slow + (1 - Rs) * inflow
    outflow = (Rs / (1 - Rs)) * x_slow
    return x_slow,outflow
end

function excess(x_loss,cmax,bexp,Pval,PETval)
    # this function calculates excess precipitation and evaporation
    xn_prev = x_loss
    ct_prev = cmax * (1 - power((1 - ((bexp + 1) * (xn_prev) / cmax)), (1 / (bexp + 1))))
    # Calculate Effective rainfall 1
    ER1 = max((Pval - cmax + ct_prev), 0.0)
    Pval = Pval - ER1
    dummy = min(((ct_prev + Pval) / cmax), 1)
    xn = (cmax / (bexp + 1)) * (1 - power((1 - dummy), (bexp + 1)))

    # Calculate Effective rainfall 2
    ER2 = max(Pval - (xn - xn_prev), 0)

    # Alternative approach
    evap = (1 - (((cmax / (bexp + 1)) - xn) / (cmax / (bexp + 1)))) * PETval  # actual ET is linearly related to the soil moisture state
    xn = max(xn - evap, 0)  # update state

    return ER1,ER2,xn


end

using BenchmarkTools

@benchmark hymod(Precip,Pet,cmax,bexp,alpha,Ks,Kq)
