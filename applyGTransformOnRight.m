function X = applyGTransformOnRight(X, i, j, values)
a = X(:, i); b = X(:, j);
X(:, i) = values(1)*a + values(2)*b;
X(:, j) = values(3)*a + values(4)*b;
