function [scores, best_values] = get_Ttransforms_iterations_only_polish(A, B, C, D, d, istar, jstar)
%% Implementation of Theorem 3 from the paper, but we do not update the
%% indices i_k, j_k of the T-transformations, just their values a_k

w = vec(C);
for i = 1:d
    for j = 1:d
        w = w - B(j, i)*kron(D(i,:)', A(:,j));
    end
end

V = D*D';
H = A'*A;
J = A'*C;

%% for one shear, use the summation formulas from the paper
if (istar == jstar)
    g1 = 0;
    g2 = 0;
    for k = 1:d
        if (k == istar)
            continue;
        end
        for i = 1:d
            for j = 1:d

                g1 = g1 - B(j, i)*B(istar, k)*V(k,i)*H(istar,j);
                g2 = g2 - B(j, i)*B(k,istar)*V(istar,i)*H(j,k);
            end

            g1 = g1 + B(istar, k)*D(k,i)*J(istar, i);
            g2 = g2 + B(k,istar)*D(istar,i)*J(k, i);
        end
    end

    r1 = 0;
    for k = 1:d
        for i = 1:d
            if (k == istar || i == istar)
                continue;
            end

            r1 = r1 + B(istar, k)*B(istar, i)*V(k,i)*H(istar,istar);
        end
    end

    r2 = 0;
    for j = 1:d
        for k = 1:d
            if (j == istar || k == istar)
                continue;
            end

            r2 = r2 + B(j, istar)*B(k, istar)*V(istar,istar)*H(k,j);
        end
    end

    r12 = 0;
    for j = 1:d
        for k = 1:d
            if (j == istar || k == istar)
                continue;
            end

            r12 = r12 + B(j, istar)*B(istar, k)*V(k,istar)*H(istar,j);
        end
    end

    q = roots([r1 (-r1-r12-g1) 0 r12+r2+g2 -r2]);
    indices = find(abs(q-conj(q))<10e-3);
    q = q(indices);
    val = zeros(1, length(q));
    index = 0;
    for a = q'
        index = index + 1;
        val(index) = 1/a^2*((1-a)^2*r2 - 2*a*(a-1)^2*r12 + a^2*(a-1)^2*r1 - 2*a^2*(a-1)*g1 - 2*a*(1-a)*g2 + a^2*(w'*w));
    end
    
    scores = 0;
    best_values = 0;
    if (isempty(val))
        scores(1,1) = inf;
        best_values(1,1) = 0;
    else
        [the_value, the_index] = min(val);
        scores(1,1) = the_value;
        best_values(1,1) = q(the_index);
    end
    
    return;
end

scores = inf(2);
best_values = inf(2);

%% for the other shear use matrix-vector operations
% for istar = 1:d
%     for jstar = istar+1:d
        g1 = 0;
        for k = 1:d
            for i = 1:d
                for j = 1:d
                    g1 = g1 + B(j, i)*B(istar, k)*V(k,i)*H(jstar,j);
                    g1 = g1 - B(j, i)*B(k, jstar)*V(istar,i)*H(k,j);
                end

                g1 = g1 - B(istar, k)*D(k,i)*J(jstar, i);
                g1 = g1 + B(k, jstar)*D(istar,i)*J(k, i);
            end
        end

        g2 = 0;
        for i = 1:d
            for j = 1:d
                g2 = g2 - B(j, i)*B(istar, jstar)*V(istar,i)*H(jstar,j);
            end

            g2 = g2 + B(istar, jstar)*D(istar,i)*J(jstar, i);
        end

        r1 = 0;
        for i = 1:d
            for j = 1:d
                r1 = r1 + B(i, jstar)*B(j, jstar)*V(istar,istar)*H(i,j);
                r1 = r1 + B(istar, i)*B(istar, j)*V(i,j)*H(jstar,jstar);
                r1 = r1 - 2*B(istar,i)*B(j, jstar)*V(istar,i)*H(jstar, j);
            end
        end

        r2 = B(istar, jstar)*B(istar, jstar)*V(istar,istar)*H(jstar,jstar);

        r12 = 0;
        for i = 1:d
            r12 = r12 - B(istar, i)*B(istar, jstar)*V(i,istar)*H(jstar,jstar);
            r12 = r12 + B(i, jstar)*B(istar, jstar)*V(istar,istar)*H(i, jstar);
        end
        
        q = roots([4*r2 6*r12 2*(r1+2*g2) 2*g1]);
        indices = find(abs(q-conj(q))<10e-3);
        q = q(indices);
        val = zeros(1, length(q));
        index = 0;
        for a = q'
            index = index + 1;
            val(index) = a^4*r2 + 2*a^3*r12 + a^2*r1 + 2*a*g1 + 2*a^2*g2 + w'*w;
        end
        if (isempty(val))
            scores(2,1) = inf;
            best_values(2,1) = 0;
        else
            [the_value, the_index] = min(val);
            scores(2,1) = the_value;
            best_values(2,1) = q(the_index);
        end
%     end
% end

% for istar = 1:d
%     for jstar = istar+1:d
        z1 = zeros(d^2, 1);
        for i = 1:d
            z1 = z1 - B(jstar, i)*kron(D(i,:)', A(:,istar));
            z1 = z1 + B(i, istar)*kron(D(jstar,:)', A(:,i));
        end
        z2 = B(jstar, istar)*kron(D(jstar,:)', A(:,istar));
        
        g1 = z1'*w;
        g2 = z2'*w;
        r1 = z1'*z1;
        r2 = z2'*z2;
        r12 = z1'*z2;
        
        q = roots([4*r2 6*r12 2*(r1+2*g2) 2*g1]);
        indices = find(abs(q-conj(q))<10e-5);
        q = q(indices);
        val = zeros(1, length(q));
        index = 0;
        for a = q'
            index = index + 1;
            val(index) = a^4*r2 + 2*a^3*r12 + a^2*r1 + 2*a*g1 + 2*a^2*g2 + w'*w;
        end
        if (isempty(val))
            scores(1,2) = inf;
            best_values(1,2) = 0;
        else
            [the_value, the_index] = min(val);
            scores(1,2) = the_value;
            best_values(1,2) = q(the_index);
        end
%     end
% end
