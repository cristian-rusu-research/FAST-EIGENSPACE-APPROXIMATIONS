function [positions, values, approx_error, tus, c] = general_approximation(C, c, m, updateSpectrum, onlyPolish)
%% Demo code for paper Construction of fast approximate eigenspaces: fast graph Fourier transforms

%% Input:
% C - a symmetric matrix of size dxd
% c - the spectrum of C (or an estimation of it)
% m - the number of T transformations to use for the approximation of
% the eigenspace

%% Output:
% The m T transformations:
% positions - the two indices (i,j) where the transformation operates
% values - the four values of the transformations
% approx_error - the approximation error, as defined in the paper
% tus - the total running time
% c - the spectrum (updated if updateSpectrum = 1)

tic;
[d, ~] = size(C);

%% basic sanity check
if (d <= 1) || (m < 1)
    positions = []; values = []; tus = toc;
    return;
end

%% make sure we have a positive integer
m = round(m);

%% vector that will store the indices (i,j) and the values of the transformations for each of the g Givens transformations
positions = zeros(2, m);
values = zeros(1, m);

%% number of iterations
K = 1;
Kinit = 1;

for k = 1:Kinit
    %% compute all scores
    B = diag(c);
    [scores, best_values] = get_Ttransforms_initialization(B, C, d);

    approx_error = norm(C - B, 'fro')^2;
    Tbar = eye(d);
    Tbarinv = eye(d);
    %% initialization of each Givens transformation
    for kk = 1:m
        %% check where the maximum scores is, to find the optimum indices
        [~, index_nuc] = min(scores(:));
        [i_nuc, j_nuc] = ind2sub([d d], index_nuc);

        %% compute the optimum T transformation on the optimum indices
        if (i_nuc == j_nuc)
            % we have a scaling
            B(i_nuc, :) = best_values(i_nuc, i_nuc)*B(i_nuc, :);
            Tbar(i_nuc, :) = best_values(i_nuc, i_nuc)*Tbar(i_nuc, :);
            B(:, i_nuc) = B(:, i_nuc)/best_values(i_nuc, i_nuc);
            Tbarinv(:, i_nuc) = Tbarinv(:, i_nuc)/best_values(i_nuc, i_nuc);
        else
            % we have a type-1 shear
            if (j_nuc > i_nuc)
                B(i_nuc, :) = B(i_nuc, :) + best_values(i_nuc, j_nuc)*B(j_nuc, :);
                Tbar(i_nuc, :) = Tbar(i_nuc, :) + best_values(i_nuc, j_nuc)*Tbar(j_nuc, :);
                B(:, j_nuc) = -best_values(i_nuc, j_nuc)*B(:, i_nuc) + B(:, j_nuc);
                Tbarinv(:, j_nuc) = -best_values(i_nuc, j_nuc)*Tbarinv(:, i_nuc) + Tbarinv(:, j_nuc);
            else
                % we have a type-2 shear
                B(i_nuc, :) = best_values(i_nuc, j_nuc)*B(j_nuc, :) + B(i_nuc, :);
                Tbar(i_nuc, :) = best_values(i_nuc, j_nuc)*Tbar(j_nuc, :) + Tbar(i_nuc, :);
                B(:, j_nuc) = B(:, j_nuc) - best_values(i_nuc, j_nuc)*B(:, i_nuc);
                Tbarinv(:, j_nuc) = Tbarinv(:, j_nuc) - best_values(i_nuc, j_nuc)*Tbarinv(:, i_nuc);
            end
        end

        %% save the Givens transformation
        positions(1, kk) = i_nuc;
        positions(2, kk) = j_nuc;
        values(:, kk) = best_values(i_nuc, j_nuc);

        %% update all the scores
        if (kk < m)
            [scores, best_values] = get_Ttransforms_initialization(B, C, d);
        end

        approx_error = [approx_error norm(C - B, 'fro')^2];
    end

    if (updateSpectrum)
        Z = zeros(d*d, d);
        for kk = 1:d
            Z(:, kk) = kron(Tbarinv(kk, :).', Tbar(:, kk));
        end
        c = Z\vec(C);
        
        [~, msgid] = lastwarn;
        if ~isempty(msgid)
            if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                c = pinv(Z)*vec(C);
            end
        end
        
        clear Z;
    end

    approx_error = [approx_error norm(C - Tbar*diag(c)*Tbarinv, 'fro')^2];
end

%% iterative process to refine the initialization
for k = 1:K
    A = eye(d);
    Ainv = eye(d);
    for tt = 1:m
        i_nuc = positions(1, tt);
        j_nuc = positions(2, tt);

        if (i_nuc == j_nuc)
            % we have a scaling
            A(i_nuc, :) = values(:, tt)*A(i_nuc, :);
            Ainv(:, i_nuc) = Ainv(:, i_nuc)/values(:, tt);
        else
            % we have a shear
            A(i_nuc, :) = A(i_nuc, :) + values(:, tt)*A(j_nuc, :);
            Ainv(:, j_nuc) = -values(:, tt)*Ainv(:, i_nuc) + Ainv(:, j_nuc);
        end
    end

    B = diag(c);
    
    for kk = 1:m
        i_nuc = positions(1, kk);
        j_nuc = positions(2, kk);

        if (i_nuc == j_nuc)
            % we have a scaling
            A(:, i_nuc) = A(:, i_nuc)/values(:, kk);
            Ainv(i_nuc, :) = Ainv(i_nuc, :)*values(:, kk);
        else
            % we have a shear
            A(:, j_nuc) = A(:, j_nuc) - values(:, kk)*A(:, i_nuc);
            Ainv(i_nuc, :) = +values(:, kk)*Ainv(j_nuc, :) + Ainv(i_nuc, :);
        end
        
        if (onlyPolish == 0)
            [scores, best_values] = get_Ttransforms_iterations(A, B, C, Ainv, d);
            
            %% check where the maximum scores is, to find the optimum indices
            [~, index_nuc] = min(scores(:));
            [i_nuc, j_nuc] = ind2sub([d d], index_nuc);
            
            %% compute the optimum T transformation on the optimum indices
            if (i_nuc == j_nuc)
                % we have a scaling
                B(i_nuc, :) = best_values(i_nuc, i_nuc)*B(i_nuc, :);
                B(:, i_nuc) = B(:, i_nuc)/best_values(i_nuc, i_nuc);
            else
                % we have a shear
                B(i_nuc, :) = B(i_nuc, :) + best_values(i_nuc, j_nuc)*B(j_nuc, :);
                B(:, j_nuc) = -best_values(i_nuc, j_nuc)*B(:, i_nuc) + B(:, j_nuc);
            end
            
            positions(1, kk) = i_nuc;
            positions(2, kk) = j_nuc;
            values(:, kk) = best_values(i_nuc, j_nuc);
        else
            [scores, best_values] = get_Ttransforms_iterations_only_polish(A, B, C, Ainv, d, i_nuc, j_nuc);
            
            if (i_nuc == j_nuc)
                % we have a scaling
                B(i_nuc, :) = best_values*B(i_nuc, :);
                B(:, i_nuc) = B(:, i_nuc)/best_values;
                values(:, kk) = best_values;
            else
                [~, index_nuc] = min(scores(:));
                [iii, jjj] = ind2sub([2 2], index_nuc);
                
                % we have a shear
                B(i_nuc, :) = B(i_nuc, :) + best_values(iii, jjj)*B(j_nuc, :);
                B(:, j_nuc) = -best_values(iii, jjj)*B(:, i_nuc) + B(:, j_nuc);
                
                values(:, kk) = best_values(iii, jjj);
            end
        end

        approx_error = [approx_error norm(C - A*B*Ainv, 'fro')^2];
    end
    
    if (updateSpectrum)
        A = eye(d);
        Ainv = eye(d);
        for tt = 1:m
            i_nuc = positions(1, tt);
            j_nuc = positions(2, tt);

            if (i_nuc == j_nuc)
                % we have a scaling
                A(i_nuc, :) = values(:, tt)*A(i_nuc, :);
                Ainv(:, i_nuc) = Ainv(:, i_nuc)/values(:, tt);
            else
                % we have a shear
                A(i_nuc, :) = A(i_nuc, :) + values(:, tt)*A(j_nuc, :);
                Ainv(:, j_nuc) = -values(:, tt)*Ainv(:, i_nuc) + Ainv(:, j_nuc);
            end
        end

        Z = zeros(d*d, d);
        for kk = 1:d
            Z(:, kk) = kron(Ainv(kk, :).', A(:, kk));
        end
        c = Z\vec(C);
        
        [~, msgid] = lastwarn;
        if ~isempty(msgid)
            if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                c = pinv(Z)*vec(C);
            end
        end
        
        clear Z;
        B = diag(c);
    end

    approx_error = [approx_error norm(C - A*B*Ainv, 'fro')^2];
end

approx_error = approx_error/norm(C,'fro')^2*100;

%% time everything
tus = toc;
