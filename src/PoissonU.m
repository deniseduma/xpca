function [f, grad, dir, h, h_inv] = PoissonU(X, U, V, M, L)

[n, d] = size(X);
k = length(U) / n;
U = (reshape(U, k, n))';
V = (reshape(V, k, d))';

[regU, regV, regV2] = regPoisson(n, d, k);
[f, O] = objFunc(X, U, V, M, L, regU, regV, regV2); 
f = f + trace(exp(O));

%compute gradient
if nargout >= 2
	exponent = exp(O);
	if ~isfinite(exponent)
		error('[PoissonU]exponent is infinite or NaN!');
	end
	if ~isfinite(exponent*V)
		error('[PoissonU]exponent*V is infinite or NaN!');
	end
	grad = -X*V + exponent*V + regU*U;
	grad = grad(:);
end

%compute dir (-h_inv * grad)
if nargout >= 3
	dir = zeros(size(grad));
	D = diag(reshape(exponent', n*d, 1));
	for index=1:n
		A = V' * D((index-1)*d + 1:index * d, (index-1)*d + 1:index * d) * V + c*a1*regU*eye(k);
		if ~isfinite(A)
			A
			error('[PoissonU]A is infinite or NaN ');
		end
		b = grad((index-1)*k + 1:index*k);
		dir((index-1)*k + 1:index*k) = (-1) * A \ b;
		Z = dir((index-1)*k + 1:index*k);
		if ~isfinite(Z)
			Z
			error('[PoissonU]Z is infinite or NaN ');
		end
	end	
end

%compute hessian
if nargout >= 4
	h = kron(eye(n), V') * D * kron(eye(n),V) + c*a1*regU*eye(n*k);
	%DEBUG
	%display(['	[PoissonU]rank(hessian) ' num2str(rank(hessian))]);
end

%compute the inv hessian
if nargout == 5
	h_inv = eye(n*k);
	for index=1:n
		B = V' * D((index-1)*d + 1:index * d, (index-1)*d + 1:index * d) * V + c*a1*regU*eye(k);
		if rank(B)~=k
			B
    			display(['rank of B ' num2str(rank(B))]);
			error('[PoissonU]Rank of B is not k!');
		end
		r1 = (index-1) * k + 1; r2 = index * k;
		c1 = (index-1) * k + 1; c2 = index * k;
		h_inv(r1:r2, c1:c2) = (B \ eye(size(B)));	
	end
end

end
