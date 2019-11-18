function X = applyGTransformOnRightInverse(X, i, j, values)
a = X(:, i); b = X(:, j);
BB = [values(1) values(3); values(2) values(4)];
BB = inv(BB);
values = BB(:);

X(:, i) = values(1)*a + values(2)*b;
X(:, j) = values(3)*a + values(4)*b;
