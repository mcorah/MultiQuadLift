module Util

export block_diag, index_dim, cross_matrix

block_diag(a) = a
function block_diag(a,b,c...)
  block_ab = [a zeros(size(a,1), size(b,2));
              zeros(size(b,1), size(a,2)) b]
  block_diag(block_ab,c...)
end

index_dim(coord, dim) = dot((coord-1), cumprod([1;dim][1:end-1])) + 1

cross_matrix(v::Vector) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

end
