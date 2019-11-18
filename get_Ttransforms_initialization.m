function [scores, best_values] = get_Ttransforms_initialization(B, C, d)
%% Implementation of Theorem 3 from the paper
scores = zeros(d);
best_values = zeros(d);

E = C-B;

RR = diag(E*B');
PP = diag(E'*B);
TT = diag(E.*B);
NN = my_norms(B, 2).^2;
MM = my_norms(B, 1).^2;
BB = diag(B).^2;
for i = 1:d
    q = roots([BB(i)-NN(i) -BB(i)+NN(i)+RR(i)-TT(i) 0 BB(i)-MM(i)-PP(i)+TT(i) -BB(i)+MM(i)]);
    indices = find(abs(q-conj(q))<10e-5);
    q = q(indices);

    val = zeros(1, length(q));
    index = 0;
    for a = q'
        index = index + 1;
        val(index) = -2*(a-1)*RR(i) - 2*(1/a-1)*PP(i) + 2*(a-1)^2/a*TT(i) + (a-1)^2*NN(i) + (1/a-1)^2*MM(i) - (a-1)^2*(a^2+1)/a^2*BB(i);
    end
    
    if (isempty(val))
        scores(i,i) = inf;
        best_values(i,i) = 1;
    else
        [the_value, the_index] = min(val);
        scores(i,i) = the_value;
        best_values(i,i) = q(the_index);
    end
end

RR = E*B';
PP = E'*B;
for i = 1:d
    for j = i+1:d
        x = B(i,i);
        t = B(i,j);
        z = B(j,i);
        y = B(j,j);
        cc = E(i,j);
        dd = E(j,i);

        q = roots([2*z^2 2*(x-y)*z MM(i)+NN(j)-2*x*y+2*cc*z PP(j,i)-RR(i,j)]);
        indices = find(abs(q-conj(q))<10e-5);
        q = q(indices);
        val = zeros(1, length(q));
        index = 0;
        for a = q'
            index = index + 1;
            val(index) = a^2*NN(j)-a^2*BB(j)+a^2*MM(i)-a^2*BB(i)+(-a*B(i, i) - a^2*B(j, i) + a*B(j, j))^2-2*a*RR(i,j) + 2*a*PP(j,i) + 2*E(i,j)*a^2*B(j, i);
        end
        if (isempty(val))
            scores(i,j) = inf;
            best_values(i,j) = 0;
        else
            [the_value, the_index] = min(val);
            scores(i,j) = the_value;
            best_values(i,j) = q(the_index);
        end
        
        q = roots([2*t^2 3*t*(y-x) 2*(dd*t+MM(j)+NN(i)-2*x*y) PP(i,j)-RR(j,i)]);
        indices = find(abs(q-conj(q))<10e-5);
        q = q(indices);
        val = zeros(1, length(q));
        index = 0;
        for a = q'
            index = index + 1;
            val(index) = a^2*NN(i)-a^2*BB(i)+a^2*MM(j)-a^2*BB(j)+(-a*B(j, j)-a^2*B(i, j)+a*B(i, i))^2-2*a*RR(j,i)+2*a*PP(i,j)+2*E(j,i)*a^2*B(i, j);
        end
        if (isempty(val))
            scores(j,i) = inf;
            best_values(j,i) = 0;
        else
            [the_value, the_index] = min(val);
            scores(j,i) = the_value;
            best_values(j,i) = q(the_index);
        end
    end
end
