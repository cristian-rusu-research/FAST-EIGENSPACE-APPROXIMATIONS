close all
clear
clc

%% random positive definite matrix
n = 128;
S = randn(n);
S = S'*S;

%% size of the approximation
g = round(n*log2(n));

%% perform the eigenvalue decomposition, get the true spectrum
[V, D] = eig(S);
[vals, inds] = sort(diag(D));
D = D(inds, inds);
S = V*D*V';

%% update the spectrum?
update_spectrum = 1;
%% only polish the result, after the initialization?
only_polish = 1;

%% call Algorithm 1, for symmetric matrices
[positions, values, approx_error, tus, Ubar] = orthogonal_approximation_for_symmetric(S, diag(D), g, update_spectrum, only_polish);

%% save results
save(['random psd n = ' num2str(n) ' g = ' num2str(g) '.mat']);
