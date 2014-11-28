module DynamicSystem

using Util

# Types
export System, Specification, Trait, SystemArray

# methods
export push, height, noise_width, state_index, noise_index, init_vals,
  set_specification, num_state, num_noise, state_indices, noise_indices,
  subsystem_dynamics, subsystem_noise, construct_local_dynamic_matrix,
  construct_local_noise_matrix, state_eye, get_states, system_dynamics,
  system_noise, subsystem_init_vals, constant_block, state_zeros

#=
Abstract types
=#
abstract AbstractSystem
abstract Trait

#=
Specification
=#
abstract Specification

get_system(x::Specification) = x.container

set_system(x::Specification, s::AbstractSystem) = x.container = s

height(x::Specification) = height(get_system(x))

noise_width(x::Specification) = noise_width(get_system(x))

state_index(x::Specification) = state_index(get_system(x))

noise_index(x::Specification) = noise_index(get_system(x))

function construct_local_dynamic_matrix(matrix, x::Specification)
  [spzeros(num_state(x), state_index(x)-1) matrix spzeros(num_state(x), height(x)-state_index(x)-num_state(x)+1)]
end

function construct_local_noise_matrix(matrix, x::Specification)
  [spzeros(num_state(x), noise_index(x)-1) matrix spzeros(num_state(x), noise_width(x)-noise_index(x)-num_noise(x)+1)]
end

function state_zeros(state_fun::Function, spec::Specification)
  width = height(spec)
  states = state_fun(spec)
  spzeros(length(states), width)
end

function state_eye(state_fun::Function, spec::Specification)
  mat = state_zeros(state_fun, spec)
  spec_start = state_index(spec)
  states = state_fun(spec)
  range = ((spec_start-1) + states[1]) : ((spec_start-1) + states[end])
  mat[:, range] = speye(length(states))
  mat
end

function constant_block(c::Vector, spec::Specification)
  width = height(spec)
  mat = spzeros(length(c), width)
  mat[:, end] = c
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
  set_system(spec,sys)
end

num_state(x::System) = num_state(x.specification)

num_noise(x::System) = num_noise(x.specification)

height(x::System) = (isdefined(x, :parent) ? height(x.parent) : num_state(x)+1)

noise_width(x::System) = (isdefined(x, :parent) ? noise_width(x.parent) : num_noise(x))

state_index(x::System) = x.state_offset + (isdefined(x, :parent) ? state_index(x.parent) : 1)

noise_index(x::System) = x.noise_offset + (isdefined(x, :parent) ? noise_index(x.parent) : 1)

state_indices(x::System) = state_index(x) : state_index(x) + num_state(x) - 1

noise_indices(x::System) = noise_index(x) : noise_index(x) + num_noise(x) - 1

subsystem_dynamics(system::System) = subsystem_dynamics(system.specification)

subsystem_noise(system::System) = subsystem_noise(system.specification)

subsystem_init_vals(system::System) = subsystem_init_vals(system.specification)

init_vals(s) = [subsystem_init_vals(s), 1]

function get_states(accessor::Function, system::System, array)
  get_states(accessor, system.specification, array)
end

function system_dynamics(system::System)
  subsystem = subsystem_dynamics(system)
  [subsystem; spzeros(1, size(subsystem,2))]
end

function system_noise(system::System)
  subsystem = subsystem_noise(system)
  [subsystem; spzeros(1, size(subsystem,2))]
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
  s.parent = x.container
  x.members = [x.members,s]
end

function push(x::SystemArray, s::Specification)
  system = System()
  set_specification(system, s)
  push(x, system)
end

function subsystem_dynamics(x::SystemArray)
  d = spzeros(num_state(x),height(x))
  for i = 1:length(x.members)
    member = x.members[i]
    d[state_indices(member),:] = subsystem_dynamics(member)
  end
  d
end

function subsystem_noise(x::SystemArray)
  d = spzeros(num_state(x),num_noise(x))
  for i = 1:length(x.members)
    member = x.members[i]
    if num_noise(member) > 0
      d[state_indices(member),:] = subsystem_noise(member)
    end
  end
  d
end

function get_states(accessor::Function, specification::SystemArray, array)
  cat(3, map((x) -> get_states(accessor, x, array), specification.members)...)
end

subsystem_init_vals(specification::SystemArray) = vcat(map(subsystem_init_vals, specification.members)...)

end
