module Quadrotor

using SDE
using Util
importall DynamicSystem

export Gains, MassParams, QuadrotorParams, MultiAgentParams, State,
  QuadrotorController, QuadrotorSpecification, Payload, noise_dynamics,
  Estimator, GlobalEstimator, RelativeEstimator

export critically_damped, generate_controller, create_quadrotor_system,
  create_multi_agent_system, num_state, num_noise, pos_states, dpos_states,
  att_states, datt_states, system_dynamics, system_noise, eval_estimator,
  init_vals

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
  link # false or spring constant
  MultiAgentParams{T}(dim::Array{T}, dist, quadrotor; relative = false,
    link = false) = new(dim, dist, quadrotor, relative, link) 
end
MultiAgentParams(dim::Real, x...) =  MultiAgentParams(dim, dim, dim, x...)

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

function generate_controller(params::QuadrotorParams, pos::State, command::State)
  gains = params.gains
  mass = params.mass
  acc_des = 1/mass.M .* 
    ( (gains.Kp*(command.p - pos.p) + gains.Kdp*(command.dp - pos.dp)) )

  att_des = 1/g .* [0 -1.0 0; 1.0 0 0] * acc_des

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
  container::System
  params::QuadrotorParams
  desired_position::Vector
  traits::Array{Trait,1}

  QuadrotorSpecification() = new()
end

num_state(::QuadrotorSpecification) = 10

num_noise(::QuadrotorSpecification) = 3

pos_states(::QuadrotorSpecification) = 1:3

dpos_states(::QuadrotorSpecification) = 4:6

att_states(::QuadrotorSpecification) = 7:8

datt_states(::QuadrotorSpecification) = 9:10

function subsystem_dynamics(quad::QuadrotorSpecification)
  if isdefined(quad, :traits)
    construct_local_dynamic_matrix(linear_dynamics(), quad) + sum(map((x)->subsystem_dynamics(x, quad), quad.traits))
  else
    construct_local_dynamic_matrix(linear_dynamics(), quad)
  end
end

subsystem_noise(quad::QuadrotorSpecification) = construct_local_noise_matrix((quad.params.mu / quad.params.mass.M) * noise_mapping, quad)

subsystem_init_vals(quad::QuadrotorSpecification) = [quad.desired_position,zeros(7)]

#=
  Payload
=#
type Payload <: Specification
  params::MassParams
  traits::Array{Trait,1}
end
num_state(::Payload) = 6
num_noise(::Payload) = 0
pos_states(::Payload) = 1:3
dpos_states(::Payload) = 4:6
att_states(::Payload) = 0:-1
datt_states(::Payload) = 0:-1
subsystem_init_vals(payload::Payload) = zeros(6)

#=
Trait subtypes
=#
#=
  Link
=#
type Link <: Trait
  spring_constant::Real
  link_vector::Vector
  source_offset::Vector
  target_offset::Vector
  target::Specification
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

  sparse(generate_controller(specification.params, pos, command))
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
  x::Array{Specification,1}
  y::Array{Specification,1}
  z::Array{Specification,1}
end

function eval_estimator(estimator::RelativeEstimator, state_fun::Function, specification::Specification)
  position = state_eye(state_fun, specification)

  function relative_estimate(target)
    target_position = state_eye(state_fun, target)
    relative = position - target_position
    target_desired = spzeros(size(relative)...)
    target_desired[:,end] = target.desired_position
    target_desired + relative
  end

  desired = specification.desired_position
  dims = Array[estimator.x, estimator.y, estimator.z]

  estimate = spzeros(size(position)...)
  for i=1:3
    temp = sum(map(relative_estimate,dims[i])) / length(dims[i])
    estimate[i,:] = temp[i,:]
  end
  estimate
end

#=
Top level functions
=#

function create_quadrotor_system(params::QuadrotorParams, set_point=[0,0,0])
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

function create_multi_agent_system(params::MultiAgentParams)
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
          this_dim = zeros(Int, 3)
          this_dim[dim_num] = 1
          indices = map(index, Array[quad_index+this_dim, quad_index-this_dim])
          quadrotors[indices]
        end

        if params.relative != false && !any([i,j,k].==1) && !any([i,j,k] .== dim)
          estimator = RelativeEstimator(map(dim_quads, 1:3)...)
          controller.position_estimator = estimator
        end
        quadrotor.traits = [controller]
      end
    end
  end
  
  system
end

end
