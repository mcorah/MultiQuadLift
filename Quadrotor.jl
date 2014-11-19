module Quadrotor

export Gains, MassParams, QuadrotorParams, criticallyDamped, generateController,
       MultiAgentParams, createQuadrotorSystem, createMultiAgentSystem

using SDE
using Util

#=
Global parameters and matrices
=#

g = 9.80665

linear_dynamics =
  [
    zeros(3,3) eye(3) zeros(3,4);
    zeros(2,6) [0 g; -g 0] zeros(2,2);
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
end

type MultiAgentParams
  dim
  dist
  quadrotor
  MultiAgentParams(dim::Real, dist, quadrotor) = new(dim * ones(Int,3), dist,
    quadrotor) 
  MultiAgentParams{T}(dim::Array{T}, dist, quadrotor) = new(dim, dist,
    quadrotor) 
end

type State
  p
  dp
  a
  da
end

#=
Helper functions
=#


indexDim(coord, dim) = dot((coord-1), cumprod([1;dim][1:end-1])) + 1

function pos(ind, height, width, block)
end

function criticallyDamped(params::MassParams, wn_xy, wn_z, wn_rp)
  prop(M,wn)  = M * wn^2
  deriv(M,wn) = M * wn
  gain_mat(f,M,W) = diagm(map(x->f(M,x), W))

  Kp = gain_mat(prop, params.M, [wn_xy, wn_xy, wn_z])
  Kdp = gain_mat(deriv, params.M, [wn_xy, wn_xy, wn_z])
  Ka = gain_mat(prop, params.I, [wn_rp, wn_rp])
  Kda = gain_mat(deriv, params.I, [wn_rp, wn_rp])

  Gains(Kp, Kdp, Ka, Kda)
end

function generateController(params::QuadrotorParams, pos::State, command::State)
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
Top level functions
=#

function createQuadrotorSystem(params::QuadrotorParams, set_point=[0,0,0])
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

  controller = generateController(params, pos, command)

  println("Sizes: $(size(linear_dynamics)), $(size(controller))")
  dynamics = block_diag(linear_dynamics,0) + [controller; zeros(1,11)]
  noise = params.mu/params.mass.M * [noise_mapping; 0 0 0]

  SDE.Model(dynamics, noise)
end

function createMultiAgentSystem(params::MultiAgentParams)
  dim = params.dim
  dist = params.dist
  quad_params = params.quadrotor

  l_quad = size(linear_dynamics, 1)
  num_quad = prod(dim)
  num_mat = num_quad * l_quad + 1

  index(x) = indexDim(x, dim)

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

        controller = generateController(quad_params, pos, command)
        quad_dynamics = [zeros(l_quad,l_quad*(ind-1)) linear_dynamics zeros(l_quad, l_quad * (num_quad-ind)+1)] + controller
        system_dynamics[1+(ind-1)*l_quad:ind*l_quad, :] = quad_dynamics
        noise_dynamics[1+(ind-1)*l_quad:ind*l_quad, 1+(ind-1)*size(noise_mapping,2):ind*size(noise_mapping,2)] = quad_params.mu/quad_params.mass.M * noise_mapping
      end
    end
  end
  
  (SDE.Model(sparse(system_dynamics), sparse(noise_dynamics)), init)

end

end
