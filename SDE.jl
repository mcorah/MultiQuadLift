module SDE

type Model
  f::Function
  g::Function
  Model(f::Matrix{Float64}, g::Matrix{Float64}) = new(x -> f*x, x -> g)
  Model(f::Function, g::Function) = new(f, g)
end

#=
Euler-Maruyama method
dX(t) = f(X(t))dt + g(X(t))dW(t)
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
  nW = size(model.g(X0),2)

  dW = sqrt(dt) * randn(nW,N)

  L = N/R
  Dt = R*dt

  Xs = zeros(length(X0),length(T))
  Xs[:,1] = X0
  sample = 1
  X = X0

  Wi = zeros(size(dW,1),1)
  for i = 2:L
    Wi[:] = mapslices(sum, dW[:, R*(i-1)+1:R*i], 2)
    X[:] = X + model.f(X)*Dt + model.g(X)*Wi

    if (i-1) % sample_rate == 0
      sample += 1
      Xs[:,sample] = X
    end
  end

  return Xs,T
end

end
