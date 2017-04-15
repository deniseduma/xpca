function [f, grad, dir, h, h_inv] = BernoulliV(X, U, V, M, L)

[n, d] = size(X);
k = length(U) / n;
U = reshape(U, n, k);
V = reshape(V, d, k);

[regU, regV, regV2] = regBernoulli(n, d, k);
[f, O] = objFunc(X, U, V, M, L, regU, regV, regV2); 
f = f + trace(log(1 + exp(O)));

%compute gradient
if nargout >= 2
	sig = sigmoid(O);
	grad = -X'*U + sig'*U + regV2*L*V;
	grad = grad(:);
end

%compute dir (-h_inv * grad)
if nargout >= 3
	dir = zeros(size(grad));
	sig2 = sig .* (1-sig);
	%D = diag(reshape(sig2, n*d, 1));
	D = reshape(sig2, n*d, 1);
	for index=1:d
		%A = U' * D((index-1)*n + 1:index*n, (index-1)*n + 1:index*n) * U + c*a2*regV*eye(k);
		A = U' * diag(D((index-1)*n + 1:index*n)) * U + c*a2*regV*eye(k);
		b = grad((index-1)*k + 1:index*k);
		dir((index-1)*k + 1:index*k) = (-1) * A \ b;
	end	


%compute hessian
if nargout >= 4
	h = kron(eye(d), U') * D * kron(eye(d), U) + c*a2*regV*eye(d*k) + c*a3*kron(eye(k), L) + c*a3*kron(eye(k), L');
end

%compute inverse hessian
if nargout == 5
	h_inv = eye(d*k);
	for index=1:d
		B = U' * D((index-1)*n + 1:index*n, (index-1)*n + 1:index*n) * U + c*a2*regV*eye(k);
		h_inv((index-1)*k + 1:index*k, (index-1)*k + 1:index*k) = (B \ eye(size(B)));	
	end
end

end
