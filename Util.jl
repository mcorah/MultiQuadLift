module Util
  block_diag(a) = a
  function block_diag(a,b,c...)
    block_ab = [a zeros(size(a,1), size(b,2));
                zeros(size(b,1), size(a,2)) b]
    block_diag(block_ab,c...)
  end

  export block_diag
end
