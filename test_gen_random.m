close all
clear
clc

%% random matrix
n = 64;
C = randn(n);

%% size of the approximation
m = round(n*log2(n));

%% update the spectrum?
update_spectrum = 1;
%% only polish the result, after the initialization?
only_polish = 1;

%% call Algorithm 1, for general matrices
[positions, values, approx_error, tus, cbar] = general_approximation(C, diag(C), m, update_spectrum, only_polish);

%% save results
save(['random general n = ' num2str(n) ' m = ' num2str(m) '.mat']);
