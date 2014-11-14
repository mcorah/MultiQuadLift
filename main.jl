using Quadrotor
using PyPlot
using SDE

tmax = 2e1
dt = 1e-4
R = 1
init = [0.2,0.2,0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,1.0]


I = 1.0
M = 1.0
wn_rp = 9.0
wn_z = 5.0
wn_xy = 2.0
mu = 0.1

mass_params = MassParams(M,I)

gains = criticallyDamped(mass_params, wn_xy, wn_z, wn_rp)

params = QuadrotorParams(mass_params, gains, mu)

system = createQuadrotorSystem(params)

@time X,T = SDE.em(system, init, tmax, dt, R, sample_rate = 100)

figure()
temp = X'
plot3D(temp[:,1],temp[:,2],temp[:,3])
