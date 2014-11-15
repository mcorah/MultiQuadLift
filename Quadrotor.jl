module Quadrotor

export Gains, MassParams, QuadrotorParams, criticallyDamped, generateController,
       createQuadrotorSystem, createMultiAgentSystem

using SDE
using Util

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

type State
  p
  dp
  a
  da
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

function createMultiAgentSystem(params)
  
end

end
