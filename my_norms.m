function y = my_norms(x, dim)
y = sqrt( sum( x .* conj( x ), dim ) );