""" mainscript """
using ReservesTypes
using JLD
@load "debugmodels.jld"
solvereservesmodel!(basemodel,true,iternum,iternum+100,1)

