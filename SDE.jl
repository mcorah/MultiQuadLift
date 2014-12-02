module SDE

type Model
  f::Function
  g::Function
  Model(f::AbstractArray{Float64}, g::AbstractArray{Float64}) = new((t,x) -> f*x, (t,x) -> g)
  Model(f::Function, g::Function) = new(f, g)
end

#=
Euler-Maruyama method
dX(t,X) = f(t,X(t))dt + g(t,X(t))dW(t)
t0 <= t <= tf
Brownian steps occur at a rate dt
R Brownian steps per step of integration

Based on examples from:
An Algorithmic Introduction to Numerical Simulation of Stochastic Differential
Equations by Desmond Higham
=#
em(f, g, X0, tf, dt, R; sample_rate = 1) =
  em(Model(f,g), X0, tf, dt, R, sample_rate)
function em(model::Model, X0, tf, dt, R; sample_rate = 1)
  T = [0:dt*sample_rate:tf]'

  N::Int64 = floor(tf/dt) + 1

  #The number of random variables is not given explicitly
  nW = size(model.g(0,X0),2)

  dW = sqrt(dt) * randn(nW,N)

  L = N/R
  Dt = R*dt

  Xs = zeros(length(X0),length(T))
  Xs[:,1] = X0
  sample = 1
  X = deepcopy(X0)

  Wi = zeros(size(dW,1))
  inter_sample = 0
  t = 0

  for i = 2:L
    t += Dt
    
    range = R*(i-1)+1:R*i
    
    for j = 1:size(dW,1)
      Wi[j] = 0
      @simd for k = range
        @inbounds Wi[j] += dW[j,k]
      end
    end

    X[:] = X + model.f(t,X)*Dt + model.f(t,X)*Wi

    inter_sample += 1
    if inter_sample == sample_rate
      inter_sample = 0
      sample += 1
      Xs[:,sample] = X
    end
  end

  return Xs,T
end

em_matrix(f, g, X0, tf, dt, R; sample_rate = 1) =
  em_matrix(Model(f,g), X0, tf, dt, R, sample_rate)
function em_matrix(model::Model, X0, tf, dt, R; sample_rate = 1)
  T = [0:dt*sample_rate:tf]'

  N::Int64 = floor(tf/dt) + 1

  #The number of random variables is not given explicitly
  nW = size(model.g(0,X0),2)

  dW = sqrt(dt) * randn(nW,N)

  L = N/R
  Dt = R*dt

  Xs = zeros(length(X0),length(T))
  Xs[:,1] = X0
  sample = 1
  X = deepcopy(X0)

  Wi = zeros(size(dW,1))
  inter_sample = 0
  t = 0

  F = model.f(0, X0)
  G = model.g(0, X0)
  for i = 2:L
    t += Dt
    
    range = R*(i-1)+1:R*i
    
    for j = 1:size(dW,1)
      Wi[j] = 0
      @simd for k = range
        @inbounds Wi[j] += dW[j,k]
      end
    end

    X[:] = X + F*Dt + G*Wi

    inter_sample += 1
    if inter_sample == sample_rate
      inter_sample = 0
      sample += 1
      Xs[:,sample] = X
    end
  end

  return Xs,T
end

end
