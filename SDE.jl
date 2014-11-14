module SDE
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
function em(f, g, X0, tf, dt, R)
  T = [0:dt:tf]'

  N = length(T)

  #The number of random variables is not given explicitly
  nW = size(g(X0),2)

  dW = sqrt(dt) * randn(nW,N)

  L = N/R
  Dt = R*dt

  Xs = zeros(length(X0),N)
  Xs[:,1] = X0

  Wi = zeros(size(dW,1),1)
  for i = 2:L
    Wi[:] = mapslices(sum, dW[:, R*(i-1)+1:R*i], 2)
    Xs[:,i] = Xs[:,i-1] + f(Xs[:,i-1])*Dt + g(Xs[:,i-1])*Wi
  end

  return Xs,T
end

type Model
  f::Function
  g::Function
  Model(fn::Matrix{Float64}, gn::Matrix{Float64}) = new(x -> fn*x, x -> gn)
end
end
