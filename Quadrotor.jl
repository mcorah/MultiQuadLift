module Quadrotor

using SDE
using Util
importall DynamicSystem

export Gains, MassParams, QuadrotorParams, MultiAgentParams, State,
  QuadrotorController, QuadrotorSpecification, Payload, noise_dynamics,
  Estimator, GlobalEstimator, RelativeEstimator, PayloadSystemParams, System,
  Specification

export critically_damped, generate_controller, create_system,
  num_state, num_noise, pos_states, dpos_states, att_states, datt_states,
  system_dynamics, system_noise, eval_estimator, init_vals, get_states

#=
Abstract types
=#
abstract Estimator

#=
Global parameters and matrices
=#

g = 9.80665

linear_dynamics(thrust = g) =
  [
    zeros(3,3) eye(3) zeros(3,4);
    zeros(2,6) [0 thrust; -thrust 0] zeros(2,2);
    zeros(1,10);
    zeros(2,8) eye(2);
    zeros(2,10)
  ]
   
control_mapping = 
  [
    zeros(5,3);
    1 0 0;
    zeros(2,3);
    0 1 0;
    0 0 1;
  ]

noise_mapping = 
  [
    zeros(3,3);
    eye(3,3);
    zeros(4,3);
  ]

#=
Parameter types
=#

type MassParams
  M
  I
end

type Gains
  Kp::Array{Float64,2}
  Kdp::Array{Float64,2}
  Ka::Array{Float64,2}
  Kda::Array{Float64,2}
end

type QuadrotorParams
  mass::MassParams
  gains::Gains
  mu
  payload_fraction
  QuadrotorParams(m, gains, mu, p=0) = new(m, gains, mu, p)
end

type MultiAgentParams
  dim
  dist
  quadrotor::QuadrotorParams
  relative # true or false
  MultiAgentParams{T}(dim::Array{T}, dist, quadrotor; relative = false) =
    new(dim, dist, quadrotor, relative) 
end
MultiAgentParams(dim::Real, x...) =  MultiAgentParams(dim, dim, dim, x...)

type PayloadSystemParams
  weight_fraction::Real
  spring_constant::Real
  multi_agent::MultiAgentParams
  PayloadSystemParams(w, s, m) = new(w, s, m)
end

type State
  p
  dp
  a
  da
end

#=
Controller design
=#

function critically_damped(params::MassParams, wn_xy, wn_z, wn_rp)
  prop(M,wn)  = M * wn^2
  deriv(M,wn) = M * wn
  gain_mat(f,M,W) = diagm(map(x->f(M,x), W))

  Kp = gain_mat(prop, params.M, [wn_xy, wn_xy, wn_z])
  Kdp = gain_mat(deriv, params.M, [wn_xy, wn_xy, wn_z])
  Ka = gain_mat(prop, params.I, [wn_rp, wn_rp])
  Kda = gain_mat(deriv, params.I, [wn_rp, wn_rp])

  Gains(Kp, Kdp, Ka, Kda)
end

function generate_controller(params::QuadrotorParams, pos::State,
  command::State; desired_thrust = g)
  gains = params.gains
  mass = params.mass
  acc_des = 1/mass.M .* 
    ( (gains.Kp*(command.p - pos.p) + gains.Kdp*(command.dp - pos.dp)) )

  att_des = 1/desired_thrust .* [0 -1.0 0; 1.0 0 0] * acc_des

  acc_att = 1/mass.I *
    ( gains.Ka * (att_des-pos.a) + gains.Kda * (-pos.da) )

  U = [
       acc_des[3,:];
       acc_att
      ]

  control_mapping*U
end


#=
Specification subtypes
=#
#=
  QuadrotorSpecification
=#
type QuadrotorSpecification <: Specification
  desired_position::Vector
  desired_thrust::Real
  container::System
  params::QuadrotorParams
  traits::Array{Trait,1}

  QuadrotorSpecification() = new([0,0,0],g)
end

num_state(::QuadrotorSpecification) = 10

num_noise(::QuadrotorSpecification) = 3

pos_states(::QuadrotorSpecification) = 1:3

dpos_states(::QuadrotorSpecification) = 4:6

att_states(::QuadrotorSpecification) = 7:8

datt_states(::QuadrotorSpecification) = 9:10

function subsystem_dynamics(quad::QuadrotorSpecification)
  mat::Array
  if isdefined(quad, :traits)
    mat = construct_local_dynamic_matrix(linear_dynamics(quad.desired_thrust), quad) + sum(map((x)->subsystem_dynamics(x, quad), quad.traits))
  else
    mat = construct_local_dynamic_matrix(linear_dynamics(quad.desired_thrust), quad)
  end
  mat[dpos_states(quad)[3],end] += quad.desired_thrust - g
  mat
end

subsystem_noise(quad::QuadrotorSpecification) = construct_local_noise_matrix((quad.params.mu / quad.params.mass.M) * noise_mapping, quad)

subsystem_init_vals(quad::QuadrotorSpecification) = [quad.desired_position,zeros(7)]

get_mass(s::QuadrotorSpecification) = s.params.mass

#=
  Payload
=#
type Payload <: Specification
  params::MassParams
  traits::Array{Trait,1}
  container::System
  Payload() = new()
end

num_state(::Payload) = 6
num_noise(::Payload) = 0
pos_states(::Payload) = 1:3
dpos_states(::Payload) = 4:6
att_states(::Payload) = 0:-1
datt_states(::Payload) = 0:-1
subsystem_init_vals(payload::Payload) = zeros(6)
get_mass(p::Payload) = p.params

function subsystem_dynamics(payload::Payload)
  gravity = [0, 0, -g]
  gravitation_acceleration = constant_block(gravity, payload)
  base_dynamics = [state_eye(dpos_states, payload); gravitation_acceleration]

  if isdefined(payload, :traits)
    base_dynamics + sum(map((x)->subsystem_dynamics(x, payload), payload.traits))
  else
    base_dynamics
  end
end

subsystem_noise(payload::Payload) = spzeros(num_state(payload),0)

#=
Trait subtypes
=#
#=
  Link
=#
type Link <: Trait
  spring_constant::Real
  length::Real
  source_offset::Vector
  target_offset::Vector
  target::Specification
  Link(sc, l) = new(sc, l, zeros(3), zeros(3))
end

function subsystem_dynamics(link::Link, specification)
  source_initial = subsystem_init_vals(specification)[pos_states(specification)]
  target_initial = subsystem_init_vals(link.target)[pos_states(link.target)]
  initial_relative = (target_initial + link.target_offset) -
                     (source_initial + link.source_offset)
  link_direction = initial_relative / norm(initial_relative)

  f0 = link.spring_constant * (norm(initial_relative) - link.length)
  perpendicular_constant = f0 / norm(initial_relative)

  source_position = state_eye(pos_states, specification) +
    constant_block(link.source_offset, specification)
  target_position = state_eye(pos_states, link.target) +
    constant_block(link.target_offset, link.target)

  relative_position = target_position - source_position

  unstretched_relative = constant_block(link.length * link_direction, specification)
  deviation = relative_position - unstretched_relative

  tangential_deviation = link_direction * (link_direction' * deviation)
  perpendicular_deviation = deviation - tangential_deviation

  force = link.spring_constant * tangential_deviation +
          perpendicular_constant * perpendicular_deviation
  acceleration = force / get_mass(specification).M

  mat = spzeros(num_state(specification), height(specification))
  mat[dpos_states(specification), :] = acceleration
  mat
end

#=
QuadrotorController
=#
type QuadrotorController <: Trait
  position_estimator::Estimator
  dposition_estimator::Estimator
  attitude_estimator::Estimator
  dattitude_estimator::Estimator
  QuadrotorController(x::Estimator) = new(x, x, x, x)
  QuadrotorController(x::Array{Estimator}) = new(x[1], x[2], x[3], x[4])
  QuadrotorController(a, b, c, d) = new(a, b, c, d)
end

function subsystem_dynamics(controller::QuadrotorController, specification::QuadrotorSpecification)
  estimators = [controller.position_estimator, controller.dposition_estimator, controller.attitude_estimator, controller.dattitude_estimator]
  states = (pos_states, dpos_states, att_states, datt_states)
  estimator_matrices = map((x) -> eval_estimator(x..., specification), zip(estimators,states))
  pos = State(estimator_matrices...)
  command = State(map((x)->spzeros(size(x)...), estimator_matrices)...)
  command.p[:,end] = specification.desired_position

  sparse(generate_controller(specification.params, pos, command;
    desired_thrust = specification.desired_thrust))
end

#=
Estimator subtypes
=#
#=
GlobalEstimator
=#
type GlobalEstimator <: Estimator
end
function eval_estimator(estimator::GlobalEstimator, state_fun::Function, specification::Specification)
  state_eye(state_fun, specification)
end

#=
RelativeEstimator
=#
type RelativeEstimator <: Estimator
  neigbors::Array{Specification,1}
end

function eval_estimator(estimator::RelativeEstimator, state_fun::Function, specification::Specification)
  position = state_eye(state_fun, specification)

  function relative_estimate(target)
    target_position = state_eye(state_fun, target)
    relative = position - target_position
    target_desired = constant_block(target.desired_position, target)
    target_desired + relative
  end

  estimate = spzeros(size(position)...)
  estimate[:,:] = sum(map(relative_estimate,estimator.neigbors)) / length(estimator.neigbors)
  estimate
end

#=
Top level functions
=#

function create_system(params::QuadrotorParams, set_point=[0,0,0])
  s = System()
  system_array = SystemArray()
  set_specification(s, system_array)
  quad = QuadrotorSpecification()
  quad.desired_position = set_point
  quad.params = params
  quad.traits = [QuadrotorController(GlobalEstimator())]
  set_specification(s, quad)
  s
end

function create_system(params::MultiAgentParams)
  system = System()
  system_array = SystemArray()
  set_specification(system, system_array)
  dim = params.dim
  dist = params.dist
  quad_params = params.quadrotor

  num_quad = prod(dim)

  index(x) = index_dim(x, dim)

  quadrotors = [(q = QuadrotorSpecification(); q.params = quad_params; q) for i=1:num_quad]

  for i = 1:num_quad
    push(system_array, quadrotors[i])
  end

  for i = 1:dim[1]
    for j = 1:dim[2]
      for k = 1:dim[3]
        quadrotor = quadrotors[index([i,j,k])]
        quadrotor.desired_position = dist * [i,j,k]

        controller = QuadrotorController(GlobalEstimator())

        function dim_quads(dim_num)
          quad_index = [i,j,k]
          # switch to absolute position at ends
          quads = Specification[]

          this_dim = zeros(Int, 3)
          this_dim[dim_num] = 1

          if quad_index[dim_num] != 1 
            push!(quads, quadrotors[index(quad_index - this_dim)])
          end
          if quad_index[dim_num] != dim[dim_num]
            push!(quads, quadrotors[index(quad_index + this_dim)])
          end
          quads
        end

        if k != 1
          estimator = RelativeEstimator(vcat(map(dim_quads, 1:3)...))
          controller.position_estimator = estimator
        end
        
        quadrotor.traits = [controller]
      end
    end
  end
  
  system
end

function create_system(params::PayloadSystemParams)
  system = System()
  system_array = SystemArray()
  set_specification(system, system_array)
  weight_fraction = params.weight_fraction
  spring_constant = params.spring_constant

  multi_params = params.multi_agent
  dim = multi_params.dim
  dist = multi_params.dist

  quad_params = multi_params.quadrotor

  num_quad = prod(dim)

  index(x) = index_dim(x, dim)

  quadrotors = [(q = QuadrotorSpecification(); q.params = quad_params;
    q.desired_thrust = (1 + weight_fraction) * g; q) for i=1:num_quad]

  payload_mass = MassParams(num_quad * weight_fraction * quad_params.mass.M, 0)
  payload = Payload()
  payload.params = payload_mass
  payload.traits = Trait[]

  for i = 1:num_quad
    push(system_array, quadrotors[i])
  end
  push(system_array, payload)

  delta_per_quad = weight_fraction * quad_params.mass.M * g / spring_constant

  # controller and initial position
  for i = 1:dim[1]
    for j = 1:dim[2]
      for k = 1:dim[3]
        quadrotor = quadrotors[index([i,j,k])]

        # static spring deviation is proportional to sum of quadrotors at and
        # above the current
        # full delta adds in deviation for each lower quadrotor
        delta = delta_per_quad * (k * dim[3] - 1/2 * (k * (k-1)))
        quadrotor.desired_position = dist * [i, j, k] + [0, 0, delta]

        controller = QuadrotorController(GlobalEstimator())

        function dim_quads(dim_num)
          quad_index = [i,j,k]
          # switch to absolute position at ends
          quads = Specification[]

          this_dim = zeros(Int, 3)
          this_dim[dim_num] = 1

          if quad_index[dim_num] != 1 
            push!(quads, quadrotors[index(quad_index - this_dim)])
          end
          if quad_index[dim_num] != dim[dim_num]
            push!(quads, quadrotors[index(quad_index + this_dim)])
          end
          quads
        end

        if k != 1
          estimator = RelativeEstimator(vcat(map(dim_quads, 1:3)...))
          controller.position_estimator = estimator
        end

        quadrotor.traits = [controller]
      end
    end
  end

  # add links
  for i = 1:dim[1]
    for j = 1:dim[2]
      for k = 1:dim[3]
        quadrotor = quadrotors[index([i, j, k])]
        source_link = Link(spring_constant, dist)
        target_link = Link(spring_constant, dist)

        quadrotor.traits = [quadrotor.traits, source_link]

        target_link.target =  quadrotor

        if k == 1
          payload_offset = dist * [i, j, 0]
          source_link.target_offset = payload_offset
          target_link.source_offset = payload_offset
          payload.traits = [payload.traits, target_link]
          source_link.target = payload
        else
          target_quadrotor = quadrotors[index([i, j, k-1])]
          source_link.target = target_quadrotor
          target_quadrotor.traits = [target_quadrotor.traits, target_link]
        end
      end
    end
  end

  system
end

end
