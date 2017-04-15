function [f, grad, dir, h, h_inv] = NormalU(X, U, V, M, L)

[n, d] = size(X);
k = length(U) / n;
U = reshape(U, n, k);
V = reshape(V, d, k);

[regU, regV, regV2] = regNormal(n, d, k);
[f, O] = objFunc(X, U, V, M, L, regU, regV, regV2);
f = f + 0.5*norm(O,'fro')^2;

%compute gradient
if nargout >= 2
	grad = -X*V + O*V + regU*U;
	grad = grad(:);
end

%compute dir (-h_inv * grad)
if nargout >= 3
	VV = V' * V;
	dir = zeros(size(grad));
	for index=1:n
		b = grad((index-1)*k + 1:index*k);
		dir((index-1)*k + 1:index*k) = (-1) * VV \ b;
	end	
end

%compute hessian
if nargout >= 4
	h = kron(eye(n), VV) + c*a1*regU*eye(n*k);
end

%compute inv hessian
if nargout == 5
	VV_inv = (VV + c*a1*regU*eye(k)) \ eye(size(VV));
	h_inv =  kron(eye(n), VV_inv);
end 

end
