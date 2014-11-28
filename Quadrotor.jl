module Quadrotor

using SDE
using Util
importall DynamicSystem

export Gains, MassParams, QuadrotorParams, MultiAgentParams, State,
  QuadrotorController, QuadrotorSpecification, Payload, noise_dynamics,
  Estimator, GlobalEstimator, RelativeEstimator

export critically_damped, generate_controller, create_quadrotor_system,
  create_multi_agent_system, num_state, num_noise, pos_states, dpos_states,
  att_states, datt_states, system_dynamics, local_noise_matrix, eval_estimator,
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
  quadrotor
  MultiAgentParams(dim::Real, dist, quadrotor) = new(dim * ones(Int,3), dist,
    quadrotor) 
  MultiAgentParams{T}(dim::Array{T}, dist, quadrotor) = new(dim, dist, quadrotor) 
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

function generate_controller(params::QuadrotorParams, pos::State, command::State)
  gains = params.gains
  mass = params.mass
  println(size(gains.Kp))
  println(size(command.p))
  println(size(pos.p))
  println(size(gains.Kdp))
  println(size(command.dp))
  println(size(pos.dp))
  acc_des = 1/mass.M .* 
    ( (gains.Kp*(command.p - pos.p) + gains.Kdp*(command.dp - pos.dp)) )

  att_des = 1/g .* [0 -1.0 0; 1.0 0 0] * acc_des

  show(size(gains.Ka))
  show(size(att_des))
  show(size(pos.a))
  show(size(gains.Kda))
  show(size(pos.da))
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

function system_dynamics(quad::QuadrotorSpecification)
  if isdefined(quad, :traits)
    local_dynamics_matrix(linear_dynamics(), quad) + sum(map((x)->system_dynamics(x, quad), quad.traits))
  else
    local_dynamics_matrix(linear_dynamics(), quad)
  end
end

noise_dynamics(quad::QuadrotorSpecification) = local_noise_matrix((quad.params.mu / quad.params.mass.M) * noise_mapping, quad)

init_vals(quad::QuadrotorSpecification) = [quad.desired_position,zeros(7)]

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
init_vals(payload::Payload) = zeros(6)

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

function system_dynamics(controller::QuadrotorController, specification::QuadrotorSpecification)
  estimators = [controller.position_estimator, controller.dposition_estimator, controller.attitude_estimator, controller.dattitude_estimator]
  states = (pos_states, dpos_states, att_states, datt_states)
  estimator_matrices = map((x) -> eval_estimator(x..., specification), zip(estimators,states))
  pos = State(estimator_matrices...)
  println(pos)
  command = State(map((x)->spzeros(size(x)...), estimator_matrices)...)
  println(command)
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
    target_position = sateEye(state_fun, target_position)
    relative = position - target_position
    target_desired = spzeros(size(relative)...)
    target_desired[:,end] = target.desired_position
    target_desired + relative
  end

  desired = specification.desired_position
  dims = [estimator.x, estimator.y, estimator.z]

  estimate = spzeros(position)
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
  p = [eye(3) zeros(3,8)]
  dp = [zeros(3,3) eye(3) zeros(3,5)]
  a = [zeros(2,6) eye(2) zeros(2,3)]
  da = [zeros(2,8) eye(2) zeros(2,1)]
  pos = State(p, dp, a, da)

  p = [zeros(3,10) set_point]
  dp = zeros(3,11)
  a = zeros(2,11)
  da = zeros(2,11)
  command = State(p, dp, a, da)

  controller = generate_controller(params, pos, command)

  println("Sizes: $(size(linear_dynamics())), $(size(controller))")
  dynamics = block_diag(linear_dynamics(),0) + [controller; zeros(1,11)]
  noise = params.mu/params.mass.M * [noise_mapping; 0 0 0]

  SDE.Model(dynamics, noise)
end

function create_multi_agent_system(params::MultiAgentParams)
  dim = params.dim
  dist = params.dist
  quad_params = params.quadrotor

  l_quad = size(linear_dynamics(), 1)
  num_quad = prod(dim)
  num_mat = num_quad * l_quad + 1

  index(x) = index_dim(x, dim)

  system_dynamics = zeros(num_mat, num_mat)
  noise_dynamics = zeros(num_mat, num_quad * size(noise_mapping,2))
  init = zeros(num_mat)
  init[end] = 1
  for i = 1:dim[1]
    for j = 1:dim[2]
      for k = 1:dim[3]
        ind = index([i,j,k])
        set_point = dist * ([i,j,k] - 1)

        init[(l_quad*(ind-1)+1):(l_quad*(ind-1)+3)] = set_point

        p = [zeros(3,(ind-1)*(l_quad)) eye(3) zeros(3,7+(num_quad-ind)*l_quad+1)]
        dp = [zeros(3,(ind-1)*(l_quad)+3) eye(3) zeros(3,4+(num_quad-ind)*l_quad+1)]
        a = [zeros(2,(ind-1)*(l_quad)+6) eye(2) zeros(2,2+(num_quad-ind)*l_quad+1)]
        da = [zeros(2,(ind-1)*(l_quad)+8) eye(2) zeros(2,(num_quad-ind)*l_quad+1)]
        pos = State(p, dp, a, da)

        p = [zeros(3,num_mat-1) set_point]
        dp = zeros(3,num_mat)
        a = zeros(2,num_mat)
        da = zeros(2,num_mat)
        command = State(p, dp, a, da)

        controller = generate_controller(quad_params, pos, command)
        quad_dynamics = [zeros(l_quad,l_quad*(ind-1)) linear_dynamics() zeros(l_quad, l_quad * (num_quad-ind)+1)] + controller
        system_dynamics[1+(ind-1)*l_quad:ind*l_quad, :] = quad_dynamics
        noise_dynamics[1+(ind-1)*l_quad:ind*l_quad, 1+(ind-1)*size(noise_mapping,2):ind*size(noise_mapping,2)] = quad_params.mu/quad_params.mass.M * noise_mapping
      end
    end
  end
  
  (SDE.Model(sparse(system_dynamics), sparse(noise_dynamics)), init)

end

end
