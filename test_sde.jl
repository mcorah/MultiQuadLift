#=
Test code plotting basic Brownian motion
=#

using SDE
using PyPlot

mu = 1

tf = 1
dt = 1e-3
R = 1

X,T = SDE.em(x->0,x->mu, 0, tf, dt, R)
figure()
plot(T',X')
