module AnalyzeQuadrotor

using PyPlot
using TestQuadrotor
using StandardQuadrotor
using Util

export plot_test, compute_error, compute_variance, compute_dav, compute_ler

function plot_test(result::Result)
  figure()
  for i = 1:size(result.data,3)
    temp = result.data[:,:,i]'
    plot3D(temp[:,1],temp[:,2],temp[:,3])
  end
end

function compute_error(result::Result, dim)
  num_quads = prod(dim)
  error = copy(result.data[:,:,1:num_quads])
  for i=1:size(error,2)
    error[:,i,:] -= result.command[:,1,1:num_quads]
  end
  error
end

function compute_variance(deviation)
  result = 0
  temp = zeros(size(deviation, 1) * size(deviation, 3))
  for i = 1:size(deviation, 2)
    temp[:] = deviation[:,i,:][:]
    result += dot(temp,temp)
  end
  result /= size(deviation, 3)
end

function compute_dav(error)
  aves = mean(error, 3)
  out = copy(error)
  for i = 1:size(error, 2)
    ave = aves[:,i,i]
    for j = 1:size(error, 3)
      out[:,i,j] -= ave
    end
  end
  out
end

function compute_ler(error, dim)
  out = copy(error)
  for i = 1:size(error, 2) # times
    for j = 1:size(error, 3) # positions
      coord = coord_dim(j, dim)
      neighbors = compute_neighbors(coord, dim)
      local_mean = mean(error[:,i,neighbors],3)
      out[:,i,j] -= local_mean
    end
  end
  out
end

end
