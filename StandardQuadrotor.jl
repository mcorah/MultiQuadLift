module StandardQuadrotor

using Quadrotor
using SDE

export Result
export multi_agent_system, relative_multi_agent_system, payload_system,
  relative_payload_system, run_test

# saving command, timeseries, times
type Result
  data
  command
  times
end

tmax = 1e1
dt = 4e-4
sample_rate = 100
R = 1

I = 1.0
M = 1.0
wn_rp = 9.0
wn_z = 5.0
wn_xy = 2.0
mu = 0.1

dist = 1

spring_constant = 1e2
weight_fraction = 0.7

function multi_agent_system(dims; relative=false)
  mass_params = MassParams(M,I)

  gains = critically_damped(mass_params, wn_xy, wn_z, wn_rp)

  quad_params = QuadrotorParams(mass_params, gains, mu)

  multi_params = MultiAgentParams(dims, dist, quad_params; relative=relative)

  payload_params = PayloadSystemParams(weight_fraction, spring_constant,
                                       multi_params)

  create_system(payload_params)
end

function payload_system(dims; relative=false)
  mass_params = MassParams(M,I)

  gains = critically_damped(mass_params, wn_xy, wn_z, wn_rp)

  quad_params = QuadrotorParams(mass_params, gains, mu)

  multi_params = MultiAgentParams(dims, dist, quad_params; relative=relative)

  payload_params = PayloadSystemParams(weight_fraction, spring_constant,
                                       multi_params)

  create_system(payload_params)
end

relative_multi_agent_system(dims) = multi_agent_system(dims; relative=true)

relative_payload_system(dims) = payload_system(dims; relative=true)

function simulate_system(system::System)
  init = init_vals(system)
  system_matrix = system_dynamics(system)
  noise_matrix = system_noise(system)

  model = SDE.Model(system_matrix, noise_matrix)
  X,T = SDE.em(model, init, tmax, dt, R, sample_rate = sample_rate)
  
  positions = get_states(pos_states, system, X)
  command = get_states(pos_states, system, init)
  Result(positions, command, T)
end

function run_test(system_function::Function, params)
  system = system_function(params)
  simulate_system(system)
end

end
