function [f,grad, dir, h, h_inv] = PoissonM(X, U, V, M, L)

[n, d] = size(X);
k = length(U) / n;
U = reshape(U, n, k);
V = reshape(V, d, k);

[regU, regV, regV2] = regPoisson(n, d, k);
[f, O] = objFunc(X, U, V, M, L, regU, regV, regV2); 
f = f + trace(exp(O));

%compute gradient
if nargout >= 2
	exponent = exp(O)'*ones(n, 1);
	if ~isfinite(exponent)
		error('[PoissonM]exponent is infinite or NaN!');
	end
	grad = -X'*ones(n, 1) + exponent;
end	

%compute dir (-h_inv * grad)
if nargout >= 3
	D = diag(exponent);
	dir = (-1) * D \ grad;
end

%compute hessian
if nargout >= 4
	h = D;
end

%compute the inv hessian
if nargout == 5
	h_inv =  diag(1 ./ exponent);
end

end
