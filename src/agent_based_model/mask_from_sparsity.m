function A_masked = mask_from_sparsity(sparsity, A)
%

[m, n]             = size(A);
total_num          = m * n;
theMask            = ones(m, n);
if sparsity > 0 && sparsity <= 1
  theMask          = theMask(:);
  total_sparse_num = round(total_num * sparsity);
% choose #sparse of them to be zero
  p_mask           = randperm(total_num, total_sparse_num);
  theMask(p_mask)  = 0;
  theMask          = reshape(theMask, m, n);
elseif sparsity < 0
  error('');
end
A_masked           = theMask .* A;
end