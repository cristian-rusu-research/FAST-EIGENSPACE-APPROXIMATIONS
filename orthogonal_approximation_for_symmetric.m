function [positions, values, approx_error, tus, Ubar, s] = orthogonal_approximation_for_symmetric(S, s, g, updateSpectrum, onlyPolish)
%% Demo code for paper Construction of fast approximate eigenspaces: fast graph Fourier transforms

%% Input:
% S - a symmetric matrix of size dxd
% s - the spectrum of S (or an estimation of it)
% g - the number of generalized Givens transformations to use for the
% approximation of the eigenspace
% updateSpectrum - 0/1 if the spectrum gets updates in the iterations
% onlyPolish - 0/1 update the values of the generalized Givens
% transformations but keep the indices found in the initialization phase
% (much faster algorithm, but less accurate approximation)

%% Output:
% The g generalized Givens transformations:
% positions - the two indices (i,j) where the transformation operates
% values - the four values of the transformations
% approx_error - the approximation error, as defined in the paper
% tus - the total running time
% Ubar - the explicit eigenspace
% s - the spectrum (updated if updateSpectrum = 1)

tic;
[d, ~] = size(S);

%% basic sanity check
if (d <= 1) || (g < 1)
    positions = []; values = []; tus = toc;
    return;
end
if norm(S-S', 'fro') >= 10e-7
    error('S has to be symmetric');
end

%% make sure we have a positive integer
g = round(g);

%% keep track of the approximation error
approx_error = [];

%% vector that will store the indices (i,j) and the values of the transformations for each of the g Givens transformations
positions = zeros(2, g);
values = zeros(4, g);

%% number of iterations
K = 2;

%% compute all scores
scores = zeros(d);
for i = 1:d
    for j = i+1:d
        aux = S(i,i)-S(j,j);
        scores(i, j) = 1/2*(aux+sqrt(aux^2+4*S(i,j)^2))*(s(j)-s(i));
    end
end

%% initialization of each Givens transformation
A = S;
for kk = g:-1:1
    %% check where the maximum scores is, to find the optimum indices
    [~, index_nuc] = max(scores(:));
    [i_nuc, j_nuc] = ind2sub([d d], index_nuc);

    %% compute the optimum orthogonal transformation on the optimum indices
    [Vv, Dd] = eig(A([i_nuc j_nuc], [i_nuc j_nuc]));
    [~, inds] = sort(diag(Dd));
    Vv = Vv(:, inds);
    GG = Vv';

    %% save the Givens transformation
    positions(1, kk) = i_nuc;
    positions(2, kk) = j_nuc;
    values(:, kk) = vec(GG);

    %% update the working matrix
    A = applyGTransformOnRightTransp(A, i_nuc, j_nuc, values(:, kk));
    A = applyGTransformOnLeft(A, i_nuc, j_nuc, values(:, kk));

    %% update the scores only for the coordinates that were selected, everything else is the same
    for i = [i_nuc j_nuc]
        for j = i+1:d
            aux = A(i,i)-A(j,j);
            scores(i, j) = 1/2*(aux+sqrt(aux^2 + 4*A(i,j)^2))*(s(j)-s(i));
        end
    end

    for j = [i_nuc j_nuc]
        for i = 1:j-1
            aux = A(i,i)-A(j,j);
            scores(i, j) = 1/2*(aux+sqrt(aux^2 + 4*A(i,j)^2))*(s(j)-s(i));
        end
    end

    approx_error = [approx_error norm(A - diag(s), 'fro')^2];
end

if (updateSpectrum)
    s = diag(A);
end

Ubar = eye(d);
for k = 1:g
    aux = values(2, k);
    values(2, k) = values(3, k);
    values(3, k) = aux;
    
    Ubar = applyGTransformOnLeft(Ubar, positions(1, k), positions(2, k), values(:, k));
end

%% iterative process to refine the initialization
for k = 1:K
    A = S;
    B = diag(s);
    for kk = 1:g
        B = applyGTransformOnLeft(B, positions(1, kk), positions(2, kk), values(:, kk));
        B = applyGTransformOnRightTransp(B, positions(1, kk), positions(2, kk), values(:, kk));
    end
    
    for kk = g:-1:1
        B = applyGTransformOnLeftTransp(B, positions(1, kk), positions(2, kk), values(:, kk));
        B = applyGTransformOnRight(B, positions(1, kk), positions(2, kk), values(:, kk));
   
        if (onlyPolish == 0)
            %% compute all the scores from scratch
            Wa = my_norms(A, 1).^2;
            Wb = my_norms(B, 1).^2;
            Z = A*B;
            V = A.*B;

            w_norm_squared_original = sum(Wa) + sum(Wb) - 2*sum(sum(V));
            scores = inf(d);
            for i = 1:d
                for j = i+1:d
                    % correction
                    w_norm_squared = w_norm_squared_original;
                    w_norm_squared = w_norm_squared - Wa(i) - Wb(i) + 2*V(i,i);
                    w_norm_squared = w_norm_squared - Wa(j) - Wb(j) + 2*V(j,j);
                    for k1 = [i j]
                        for k2 = k1+1:d
                            w_norm_squared = w_norm_squared + 4*V(k1,k2);
                        end
                    end
                    for k1 = setdiff(1:d, [i j])
                        for k2 = [i j]
                            if (k2 > k1)
                                w_norm_squared = w_norm_squared + 4*V(k1,k2);
                            end
                        end
                    end

                    %% first solution
                    P = [Wa(i)+Wb(i)-2*V(i,i)+Wa(j)+Wb(j)-2*V(j,j)-2*(2*V(i,j)) -2*A(i,i)*B(i,j)+2*A(i,j)*B(i,i)-2*A(i,j)*B(j,j)+2*A(j,j)*B(i,j); 0 Wa(i)+Wb(j)-2*A(i,i)*B(j,j)+Wa(j)+Wb(i)-2*A(j,j)*B(i,i)-2*(-2*V(i,j))];
                    P(2,1) = P(1,2);
                    gg = [2*(-Z(i,i)+V(i,i)+2*V(i,j)-Z(j,j)+V(j,j)); -2*Z(i,j)+2*A(i,j)*B(j,j)+2*A(i,i)*B(i,j)+2*Z(j,i)-2*A(j,j)*B(i,j)-2*A(i,j)*B(i,i)];
                    my_lambda = eigs([P*P-gg*gg' zeros(2); zeros(2) eye(2)], [2*P -eye(2); eye(2) zeros(2)], 1, 'SR');
                    x = real(-(P-my_lambda*eye(2)) \ gg);
                    
                    [~, msgid] = lastwarn;
                    if ~isempty(msgid)
                        if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                            x = real(pinv(-(P-my_lambda*eye(2)))*gg);
                        end
                    end
                    
                    scores(i, j) = w_norm_squared+2*x'*gg+x'*P*x;

                    %% second solution
                    P = [Wa(i)+Wb(i)-2*V(i,i)+Wa(j)+Wb(j)-2*V(j,j)-2*(-2*V(i,j)) -2*A(i,i)*B(i,j)-2*A(i,j)*B(i,i)+2*A(i,j)*B(j,j)+2*A(j,j)*B(i,j); 0 Wa(i)+Wb(j)-2*A(i,i)*B(j,j)+Wa(j)+Wb(i)-2*A(j,j)*B(i,i)+2*(-2*V(i,j))];
                    P(2,1) = P(1,2);
                    gg = [2*(-Z(i,i)+V(i,i)+Z(j,j)-V(j,j)); -2*Z(i,j)+2*A(i,j)*B(j,j)+2*A(i,i)*B(i,j)-2*Z(j,i)+2*A(j,j)*B(i,j)+2*A(i,j)*B(i,i)];
                    my_lambda = eigs([P*P-gg*gg' zeros(2); zeros(2) eye(2)], [2*P -eye(2); eye(2) zeros(2)], 1, 'SR');
                    x = real(-(P-my_lambda*eye(2)) \ gg);
                    
                    [~, msgid] = lastwarn;
                    if ~isempty(msgid)
                        if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                            x = real(pinv(-(P-my_lambda*eye(2)))*gg);
                        end
                    end
                    
                    scores(j, i) = w_norm_squared+2*x'*gg+x'*P*x;
                end
            end
            
            %% check where the maximum scores is, to find the optimum indices
            [~, index_nuc] = min(scores(:));
            [i_nuc, j_nuc] = ind2sub([d d], index_nuc);

            i = min(i_nuc, j_nuc);
            j = max(i_nuc, j_nuc);
            
            %% compute the optimum orthogonal transformation on the optimum indices
            if (j_nuc > i_nuc)
                P = [Wa(i)+Wb(i)-2*V(i,i)+Wa(j)+Wb(j)-2*V(j,j)-2*(2*V(i,j)) -2*A(i,i)*B(i,j)+2*A(i,j)*B(i,i)-2*A(i,j)*B(j,j)+2*A(j,j)*B(i,j); 0 Wa(i)+Wb(j)-2*A(i,i)*B(j,j)+Wa(j)+Wb(i)-2*A(j,j)*B(i,i)-2*(-2*V(i,j))];
                P(2,1) = P(1,2);
                gg = [2*(-Z(i,i)+V(i,i)+2*V(i,j)-Z(j,j)+V(j,j)); -2*Z(i,j)+2*A(i,j)*B(j,j)+2*A(i,i)*B(i,j)+2*Z(j,i)-2*A(j,j)*B(i,j)-2*A(i,j)*B(i,i)];
                my_lambda = eigs([P*P-gg*gg' zeros(2); zeros(2) eye(2)], [2*P -eye(2); eye(2) zeros(2)], 1, 'SR');
                x = -(P-my_lambda*eye(2)) \ gg;
                
                [~, msgid] = lastwarn;
                if ~isempty(msgid)
                    if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                        x = real(pinv(-(P-my_lambda*eye(2)))*gg);
                    end
                end
                
                GG = [x(1) x(2); -x(2) x(1)];
            else
                P = [Wa(i)+Wb(i)-2*V(i,i)+Wa(j)+Wb(j)-2*V(j,j)-2*(-2*V(i,j)) -2*A(i,i)*B(i,j)-2*A(i,j)*B(i,i)+2*A(i,j)*B(j,j)+2*A(j,j)*B(i,j); 0 Wa(i)+Wb(j)-2*A(i,i)*B(j,j)+Wa(j)+Wb(i)-2*A(j,j)*B(i,i)+2*(-2*V(i,j))];
                P(2,1) = P(1,2);
                gg = [2*(-Z(i,i)+V(i,i)+Z(j,j)-V(j,j)); -2*Z(i,j)+2*A(i,j)*B(j,j)+2*A(i,i)*B(i,j)-2*Z(j,i)+2*A(j,j)*B(i,j)+2*A(i,j)*B(i,i)];
                my_lambda = eigs([P*P-gg*gg' zeros(2); zeros(2) eye(2)], [2*P -eye(2); eye(2) zeros(2)], 1, 'SR');
                x = -(P-my_lambda*eye(2)) \ gg;
                
                [~, msgid] = lastwarn;
                if ~isempty(msgid)
                    if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                        x = real(pinv(-(P-my_lambda*eye(2)))*gg);
                    end
                end
                    
                GG = [x(1) x(2); x(2) -x(1)];
            end
        else
            i = positions(1, kk); i_nuc = i;
            j = positions(2, kk); j_nuc = j;
            
            %% compute all the scores from scratch
            Wa = my_norms(A(:, [i j]), 1).^2;
            Wb = my_norms(B(:, [i j]), 1).^2;
            Z = A([i j], :)*B(:, [i j]);
            V = A([i j], [i j]).*B([i j], [i  j]);
            
            %% compute the optimum orthogonal transformation on the optimum indices
            P = [Wa(1)+Wb(1)-2*V(1,1)+Wa(2)+Wb(2)-2*V(2,2)-2*(2*V(1,2)) -2*A(i,i)*B(i,j)+2*A(i,j)*B(i,i)-2*A(i,j)*B(j,j)+2*A(j,j)*B(i,j); 0 Wa(1)+Wb(2)-2*A(i,i)*B(j,j)+Wa(2)+Wb(1)-2*A(j,j)*B(i,i)-2*(-2*V(1,2))];
            P(2,1) = P(1,2);
            gg = [2*(-Z(1,1)+V(1,1)+2*V(1,2)-Z(2,2)+V(2,2)); -2*Z(1,2)+2*A(i,j)*B(j,j)+2*A(i,i)*B(i,j)+2*Z(2,1)-2*A(j,j)*B(i,j)-2*A(i,j)*B(i,i)];
            my_lambda = eigs([P*P-gg*gg' zeros(2); zeros(2) eye(2)], [2*P -eye(2); eye(2) zeros(2)], 1, 'SR');
            x = -(P-my_lambda*eye(2)) \ gg;
            
            [~, msgid] = lastwarn;
            if ~isempty(msgid)
                if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                    x = real(pinv(-(P-my_lambda*eye(2)))*gg);
                end
            end
            
            GG1 = [x(1) x(2); -x(2) x(1)];
            val1 = norm(gg + P*x, 'fro');

            P = [Wa(1)+Wb(1)-2*V(1,1)+Wa(2)+Wb(2)-2*V(2,2)-2*(-2*V(1,2)) -2*A(i,i)*B(i,j)-2*A(i,j)*B(i,i)+2*A(i,j)*B(j,j)+2*A(j,j)*B(i,j); 0 Wa(1)+Wb(2)-2*A(i,i)*B(j,j)+Wa(2)+Wb(1)-2*A(j,j)*B(i,i)+2*(-2*V(1,2))];
            P(2,1) = P(1,2);
            gg = [2*(-Z(1,1)+V(1,1)+Z(2,2)-V(2,2)); -2*Z(1,2)+2*A(i,j)*B(j,j)+2*A(i,i)*B(i,j)-2*Z(2,1)+2*A(j,j)*B(i,j)+2*A(i,j)*B(i,i)];
            my_lambda = eigs([P*P-gg*gg' zeros(2); zeros(2) eye(2)], [2*P -eye(2); eye(2) zeros(2)], 1, 'SR');
            x = -(P-my_lambda*eye(2)) \ gg;
            
            [~, msgid] = lastwarn;
            if ~isempty(msgid)
                if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                    x = real(pinv(-(P-my_lambda*eye(2)))*gg);
                end
            end
            
            GG2 = [x(1) x(2); x(2) -x(1)];
            val2 = norm(gg + P*x, 'fro');
        end
        
        %% save the Givens transformation
        if (abs(norm(x) -  1) <= 10e-5)
            positions(1, kk) = i;
            positions(2, kk) = j;
            
            if (onlyPolish == 0)
                values(:, kk) = vec(GG);
            else
                if (val1 < val2)
                    values(:, kk) = vec(GG1);
                else
                    values(:, kk) = vec(GG2);
                end
            end
        end
        
        %% update the working matrices
        A = applyGTransformOnLeftTransp(A, positions(1, kk), positions(2, kk), values(:, kk));
        A = applyGTransformOnRight(A, positions(1, kk), positions(2, kk), values(:, kk));

        approx_error = [approx_error norm(A - B, 'fro')^2];
    end
    
    if (updateSpectrum)
        s = diag(B);
    end
end

%% the explicit approximation, can be avoided
Ubar = eye(d);
for k = 1:g
    Ubar = applyGTransformOnLeft(Ubar, positions(1, k), positions(2, k), values(:, k));
end

%% normalize errors
approx_error = approx_error/norm(S,'fro')^2*100;

%% time everything
tus = toc;
