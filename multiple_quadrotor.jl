using Quadrotor
using PyPlot
using SDE

tmax = 1e1
dt = 4e-4
R = 1

I = 1.0
M = 1.0
wn_rp = 9.0
wn_z = 5.0
wn_xy = 2.0
mu = 0.2

dims = [3,3,3]
dist = 0.5

mass_params = MassParams(M,I)

gains = critically_damped(mass_params, wn_xy, wn_z, wn_rp)

quad_params = QuadrotorParams(mass_params, gains, mu)

multi_params = MultiAgentParams(dims, dist, quad_params; relative=true)

system = create_multi_agent_system(multi_params)

system_matrix = system_dynamics(system)
noise_matrix = system_noise(system)
init = init_vals(system)

model = SDE.Model(system_matrix, noise_matrix)

@time X,T = SDE.em(model, init, tmax, dt, R, sample_rate = 100)

figure()
len = size(Quadrotor.linear_dynamics(),1)
for i=1:prod(dims)
  temp = X[1+(i-1)*len:3+(i-1)*len, :]'
  plot3D(temp[:,1],temp[:,2],temp[:,3])
end
