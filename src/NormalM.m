function [f, grad, dir, h, h_inv] = NormalM(X, U, V, M, L)

[n, d] = size(X);
k = length(U) / n;
U = reshape(U, n, k);
V = reshape(V, d, k);

[regU, regV, regV2] = regNormal(n, d, k);
[f, O] = objFunc(X, U, V, M, L, regU, regV, regV2); 
f = f + 0.5*norm(O,'fro')^2;

%compute gradient
if nargout >= 2
	grad = -X'*ones(n, 1) + O'*ones(n, 1);
end	

%compute dir (-h_inv * grad)
if nargout >= 3
	dir = (-1) * eye(d) \ grad;
end

%compute hessian
if nargout >= 4
	h = eye(d);
end

%compute the inv hessian
if nargout == 5
	h_inv =  eye(d);
end

end
