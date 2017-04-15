function [f, grad, dir, h, h_inv] = PoissonV(X, U, V, M, L)

[n, d] = size(X);
k = length(U) / n;
U = reshape(U, n, k);
V = reshape(V, d, k);

[regU, regV, regV2] = regPoisson(n, d, k);
[f, O] = objFunc(X, U, V, M, L, regU, regV, regV2); 
f = f + trace(exp(O));

%compute gradient
if nargout >= 2
	exponent = exp(O);
	if ~isfinite(exponent)
		error('[PoissonV]exponent is infinite or NaN!');
	end
	grad = -X'*U + exponent'*U + regV2*L*V;
	grad = grad(:);
end

%compute dir (-h_inv * grad)
if nargout >= 3
	dir = zeros(size(grad));
	D = diag(reshape(exponent, n*d, 1));
	for index=1:d
		A = U' * D((index-1)*n + 1:index*n, (index-1)*n + 1:index*n) * U + c*a2*regV*eye(k);
		b = grad((index-1)*k + 1:index*k);
		dir((index-1)*k + 1:index*k) = (-1) * A \ b;
	end	

%compute hessian
if nargout >= 4
	h = kron(eye(d), U') * D * kron(eye(d), U) + c*a2*regV*eye(d*k) + c*a3*regV2*kron(eye(k), L) + c*a3*regV2*kron(eye(k), L');
end

%compute the inv hessian
if nargout == 5
	h_inv = eye(d*k);
	for index=1:d
		B = U' * D((index-1)*n + 1:index*n, (index-1)*n + 1:index*n) * U + c*a2*regV*eye(k);
		if rank(B)~=k
			B
    			display(['rank of B ' num2str(rank(B))]);
			error('[PoissonV]Rank of B is not k!');
		end
		r1 = (index-1) * k + 1; r2 = index * k;
		c1 = (index-1) * k + 1; c2 = index * k;
		h_inv(r1:r2, c1:c2) = (B \ eye(size(B)));	
	end
end

end
