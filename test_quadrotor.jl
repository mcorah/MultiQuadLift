using DynamicSystem
using Quadrotor
using SDE
using PyPlot
s = System()
arr = SystemArray(s)
set_specification(s,arr)

for i=1:2
  for j=1:2
    for k=1:2
      temp_system = System()
      temp_system.parent = s

      q = QuadrotorSpecification()
      mass_params = MassParams(1,1)
      gains = critically_damped(mass_params, 2, 5, 9)
      quad_params = QuadrotorParams(mass_params, gains, 0.1)
      q.params = quad_params
      q.desired_position = 1.0 * [i,j,k]

      controller = QuadrotorController(GlobalEstimator())
      q.traits = [controller]

      set_specification(temp_system,q)
      push(arr,temp_system)
    end
  end
end
  
dyn = system_dynamics(s)
noise = system_noise(s)
println("dynamics")
show(dyn); println()
println("noise")
show(noise)

init = init_vals(s)
system = SDE.Model(dyn, noise)

@time X,T = SDE.em(system, init, 2e1, 1e-4, 1, sample_rate = 100)

figure()
poss = get_states(pos_states, s, X)

for i = 1:size(poss,3)
  temp = poss[:,:,i]'
  plot3D(temp[:,1],temp[:,2],temp[:,3])
end
