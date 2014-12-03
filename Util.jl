module Util

export block_diag, index_dim, coord_dim, compute_neighbors, cross_matrix

block_diag(a) = a
function block_diag(a,b,c...)
  block_ab = [a zeros(size(a,1), size(b,2));
              zeros(size(b,1), size(a,2)) b]
  block_diag(block_ab,c...)
end

index_dim(coord, dim) = dot((coord-1), cumprod([1;dim][1:end-1])) + 1

function coord_dim(index, dim)
  temp = index - 1
  sizes = cumprod([1;dim][1:end-1])
  coord = zeros(dim)
  for i = size(dim,1):-1:1
    coord[i] = div(temp, sizes[i]) + 1
    temp %= sizes[i]
  end
  coord
end

function compute_neighbors(coord, dim)
  neighbors = Int64[]
  n_coord = zeros(coord)
  for i = 1:length(dim)
    if coord[i] !=  1
      n_coord[:] = coord
      n_coord[i] -= 1
      push!(neighbors, index_dim(n_coord, dim))
    end

    if coord[i] !=  dim[i]
      n_coord[:] = coord
      n_coord[i] += 1
      push!(neighbors, index_dim(n_coord, dim))
    end
  end
  neighbors
end

cross_matrix(v::Vector) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

end
