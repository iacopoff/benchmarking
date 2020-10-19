'''
Copyright (c) 2015 by Tobias Houska
This file is part of Statistical Parameter Estimation Tool (SPOTPY).
:author: Tobias Houska
:paper: Houska, T., Kraft, P., Chamorro-Chavez, A. and Breuer, L.: 
SPOTting Model Parameters Using a Ready-Made Python Package, 
PLoS ONE, 10(12), e0145180, doi:10.1371/journal.pone.0145180, 2015.

'''
# import time
# def timer_func(func):

#    def function_timer(*args, **kwargs):
#         start = time.time()
#         value = func(*args, **kwargs)
#         end = time.time()
#         runtime = end - start
#         msg = "{func} took {time} microseconds to complete its execution."
#         print(msg.format(func = func.__name__,time = runtime*10**6))
#         return value
#    return function_timer

# import pandas as pd


# input = "/home/iff/research/dev/nsgaii/data/hymod/hymod_input.csv"


# df = pd.read_csv(input,sep=";")
# df.columns = ["date","precip","TURC","dis"]
# df.dropna(inplace=True)
# df

# cmax  = 412.33
# bexp  = 0.1725
# alpha = 0.8127
# Ks    = 0.0404
# Kq    = 0.5592

# Precip = df.precip.values
# Pet = df.TURC.values



@timer_func
def hymod(Precip, PET, cmax,bexp,alpha,Rs,Rq):
    """
    See https://www.proc-iahs.net/368/180/2015/piahs-368-180-2015.pdf for a scientific paper.

    :param cmax:
    :param bexp:
    :param alpha:
    :param Rs:
    :param Rq:
    :return: Dataset of water in hymod (has to be calculated in litres)
    :rtype: list
    """

    # HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
    x_loss = 0.0
    # Initialize slow tank state
    x_slow = 2.3503 / (Rs * 22.5)
    x_slow = 0  # --> works ok if calibration data starts with low discharge
    # Initialize state(s) of quick tank(s)
    x_quick = [0,0,0]
    t = 0
    outflow = []
    output = []
    # START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS

    while t <= len(Precip)-1:
        Pval = Precip[t]
        PETval = PET[t]
        # Compute excess precipitation and evaporation
        ER1, ER2, x_loss = excess(x_loss, cmax, bexp, Pval, PETval)
        # Calculate total effective rainfall
        ET = ER1 + ER2
        #  Now partition ER between quick and slow flow reservoirs
        UQ = alpha * ET
        US = (1 - alpha) * ET
        # Route slow flow component with single linear reservoir
        x_slow, QS = linres(x_slow, US, Rs)
        # Route quick flow component with linear reservoirs
        inflow = UQ

        for i in range(3):
            # Linear reservoir
            x_quick[i], outflow = linres(x_quick[i], inflow, Rq)
            inflow = outflow

        # Compute total flow for timestep
        output.append(QS + outflow)
        t = t+1


    return output

def power(X,Y):
    X=abs(X) # Needed to capture invalid overflow with netgative values
    return X**Y


def linres(x_slow,inflow,Rs):
    # Linear reservoir
    x_slow = (1 - Rs) * x_slow + (1 - Rs) * inflow
    outflow = (Rs / (1 - Rs)) * x_slow
    return x_slow,outflow


def excess(x_loss,cmax,bexp,Pval,PETval):
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


# if __name__ == '__main__':
   
#     hymod(Precip, Pet, cmax,bexp,alpha,Ks,Kq)