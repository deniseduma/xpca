function [f,grad, dir, h, h_inv] = NormalV(X, U, V, M, L)

[n, d] = size(X);
k = length(U) / n;
U = reshape(U, n, k);
V = reshape(V, d, k);

[regU, regV, regV2] = regNormal(n, d, k);
[f, O] = objFunc(X, U, V, M, L, regU, regV, regV2);
f = f + 0.5*norm(O,'fro')^2;

%compute gradient
if nargout >= 2
	grad = -X'*U + O'*U + regV2*L*V;
	grad = grad(:);
end

%compute dir (-h_inv * grad)
if nargout >= 3
	UU = U' * U;
	dir = zeros(size(grad));
	for index=1:d
		b = grad((index-1)*k + 1:index*k);
		dir((index-1)*k + 1:index*k) = (-1) * UU \ b;
	end	
end

%compute hessian
if nargout >= 4
	h = kron(eye(d), UU) + c*a2*regV*eye(d*k) + c*a3*kron(eye(k), L) + c*a3*kron(eye(k), L');
end

%compute inv hessian
if nargout == 5
	UU_inv = (UU + c*a2*regV*eye(k)) \ eye(size(UU));
	h_inv =  kron(eye(d), UU_inv);
end 

end
