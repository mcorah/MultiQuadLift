using Quadrotor
using PyPlot
using SDE

tmax = 2e1
dt = 1e-4
R = 1

I = 1.0
M = 1.0
wn_rp = 9.0
wn_z = 5.0
wn_xy = 2.0
mu = 0.1

mass_params = MassParams(M,I)

gains = critically_damped(mass_params, wn_xy, wn_z, wn_rp)

params = QuadrotorParams(mass_params, gains, mu)

system = create_system(params, [0.3,0.3,0.3])

system_matrix = system_dynamics(system)
noise_matrix = system_noise(system)
init = init_vals(system)

model = SDE.Model(system_matrix, noise_matrix)

@time X,T = SDE.em(model, init, tmax, dt, R, sample_rate = 100)

figure()
temp = X'
plot3D(temp[:,1],temp[:,2],temp[:,3])
