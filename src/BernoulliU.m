function [f, grad, dir, h, h_inv] = BernoulliU(X, U, V, M, L)

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
	grad = -X*V + sig*V + regU*U;
	grad = grad(:);
end

%compute dir (-h_inv * grad)
if nargout >= 3
	dir = zeros(size(grad));
	sig2  = sig .* (1-sig);
	%D = diag(reshape(sig2', n*d, 1));
	D = reshape(sig2', n*d, 1);
	for index=1:n
		%A = V' * D((index-1)*d + 1:index * d, (index-1)*d + 1:index * d) * V + c*a1*regU*eye(k);
		A = V' * diag(D((index-1)*d + 1:index * d)) * V + c*a1*regU*eye(k);
		if rank(A) ~= k 
		        A
			disp(['rank(A) is ' num2str(rank(A))]);
			error('rank(A) is not k!');
		end
		b = grad((index-1)*k + 1:index*k);
		dir((index-1)*k + 1:index*k) = (-1) * A \ b;
	end	
end

%compute hessian
if nargout >= 4
	h = kron(eye(n), V') * D * kron(eye(n),V) + c*a1*regU*eye(n*k);
end

if nargout == 5
	h_inv = eye(n*k);
	for index=1:n
		B = V' * D((index-1)*d + 1:index*d, (index-1)*d + 1:index*d) * V + c*a1*regU*eye(k);
		if rank(B)~=k
    			display(['rank of B ' num2str(rank(B))]);
			error('[BernoulliU]Rank of B is not k!');
		end
		h_inv((index-1)*k + 1:index*k, (index-1)*k + 1:index*k) = (B \ eye(size(B)));	
	end
end

end
