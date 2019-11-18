function X = applyGTransformOnRightTransp(X, i, j, values)
a = X(:, i); b = X(:, j);
X(:, i) = values(1)*a + values(3)*b;
X(:, j) = values(2)*a + values(4)*b;
