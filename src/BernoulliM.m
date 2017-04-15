function [f, grad, dir, h, h_inv] = BernoulliM(X, U, V, M, L)

[n, d] = size(X);
k = length(U) / n;
U = reshape(U, n, k);
V = reshape(V, d, k);

[regU, regV, regV2] = regBernoulli(n, d, k);
[f, O] = objFunc(X, U, V, M, L, regU, regV, regV2); 
f = f + trace(log(1 + exp(O)));

%compute gradient
if nargout >= 2
	sig = -X'*ones(n ,1) + sigmoid(O)'*ones(n, 1);
	grad = grad(:);
end	

%compute dir (-h_inv * grad)
if nargout >= 3
	sig2  = sig .* (1-sig);
	D = diag(sig2);
	dir = (-1) * D \ grad;
end

%compute hessian
if nargout >= 4
	h = D;
end

%compute the inv hessian
if nargout == 5
	h_inv =  diag(1 ./ sig2);
end

end
