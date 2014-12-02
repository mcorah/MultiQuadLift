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
mu = 0.1

dims = [1,1,10]
dist = 1

mass_params = MassParams(M,I)

gains = critically_damped(mass_params, wn_xy, wn_z, wn_rp)

quad_params = QuadrotorParams(mass_params, gains, mu)

multi_params = MultiAgentParams(dims, dist, quad_params; relative=true)

system = create_system(multi_params)

system_matrix = system_dynamics(system)
noise_matrix = system_noise(system)
init = init_vals(system)

model = SDE.Model(system_matrix, noise_matrix)

@time X,T = SDE.em(model, init, tmax, dt, R, sample_rate = 100)

figure()
poss = get_states(pos_states, system, X)

for i = 1:size(poss,3)
  temp = poss[:,:,i]'
  plot3D(temp[:,1],temp[:,2],temp[:,3])
end
