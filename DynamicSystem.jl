module DynamicSystem

using Util

# Types
export System, Specification, Trait, SystemArray

# methods
export push, height, noise_width, state_index, noise_index, init_vals,
  set_specification, num_state, num_noise, state_indices, noise_indices,
  system_dynamics, noise_dynamics, local_dynamics_matrix, local_noise_matrix,
  state_eye, get_states

#=
Abstract types
=#
abstract AbstractSystem
abstract Trait

#=
Specification
=#
abstract Specification

getSystem(x::Specification) = x.container

setSystem(x::Specification, s::AbstractSystem) = x.container = s

height(x::Specification) = height(getSystem(x))

noise_width(x::Specification) = noise_width(getSystem(x))

state_index(x::Specification) = state_index(getSystem(x))

noise_index(x::Specification) = noise_index(getSystem(x))

function local_dynamics_matrix(matrix, x::Specification)
  [spzeros(num_state(x), state_index(x)-1) matrix spzeros(num_state(x), height(x)-state_index(x)-num_state(x)+1)]
end

function local_noise_matrix(matrix, x::Specification)
  [spzeros(num_state(x), noise_index(x)-1) matrix spzeros(num_state(x), noise_width(x)-noise_index(x)-num_noise(x)+1)]
end

function state_eye(state_fun::Function, spec::Specification)
  spec_start = state_index(spec)
  width = height(spec)
  states = state_fun(spec)
  mat = spzeros(length(states), width)
  range = ((spec_start-1) + states[1]) : ((spec_start-1) + states[end])
  mat[:, range] = speye(length(states))
  mat
end

function get_states(accessor::Function, specification::Specification, array)
  range = accessor(specification) + state_index(specification) - 1
  array[range, :]
end


#=
System
=#
type System <: AbstractSystem
  state_offset::Integer
  noise_offset::Integer
  specification::Specification
  parent::System
  System(s = 0, n = 0) = new(s,n)
end

function set_specification(sys::System, spec::Specification)
  sys.specification = spec
  setSystem(spec,sys)
end

num_state(x::System) = num_state(x.specification)

num_noise(x::System) = num_noise(x.specification)

height(x::System) = (isdefined(x, :parent) ? height(x.parent) : num_state(x)+1)

noise_width(x::System) = (isdefined(x, :parent) ? noise_width(x.parent) : num_noise(x))

state_index(x::System) = x.state_offset + (isdefined(x, :parent) ? state_index(x.parent) : 1)

noise_index(x::System) = x.noise_offset + (isdefined(x, :parent) ? noise_index(x.parent) : 1)

state_indices(x::System) = state_index(x) : state_index(x) + num_state(x) - 1

noise_indices(x::System) = noise_index(x) : noise_index(x) + num_noise(x) - 1

system_dynamics(system::System) = system_dynamics(system.specification)

noise_dynamics(system::System) = noise_dynamics(system.specification)

init_vals(system::System) = init_vals(system.specification)

function get_states(accessor::Function, system::System, array)
  get_states(accessor, system.specification, array)
end


#=
SystemArray
=#
type SystemArray <: Specification
  members::Array{System,1}
  container::System
  SystemArray() = new([])
  SystemArray(c::System) = new(System[],c)
end

num_state(x::SystemArray) = mapreduce(num_state,+,0,x.members)

num_noise(x::SystemArray) = mapreduce(num_noise,+,0,x.members)

function push(x::SystemArray, s::System)
  s.state_offset = mapreduce(num_state,+,0,x.members)
  s.noise_offset = mapreduce(num_noise,+,0,x.members)
  x.members = [x.members,s]
end

function system_dynamics(x::SystemArray)
  d = spzeros(num_state(x),height(x))
  for i = 1:length(x.members)
    member = x.members[i]
    d[state_indices(member),:] = system_dynamics(member)
  end
  d
end

function noise_dynamics(x::SystemArray)
  d = spzeros(num_state(x),num_noise(x))
  for i = 1:length(x.members)
    member = x.members[i]
    d[state_indices(member),:] = noise_dynamics(member)
  end
  d
end

function get_states(accesor::Function, specification::SystemArray, array)
  cat(3, map((x) -> get_states(accesor, x, array), specification.members)...)
end

init_vals(specification::SystemArray) = vcat(map(init_vals, specification.members)...)

end
