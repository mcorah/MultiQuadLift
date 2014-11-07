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

    for i = 2:L
      X = Xs[:,i-1]
      Wi = mapslices(sum, dW[:, R*(i-1)+1:R*i], 2)
      X = X + f(X)*Dt + g(X)*Wi
      Xs[:,i] = X
    end

    return Xs,T
  end
end
