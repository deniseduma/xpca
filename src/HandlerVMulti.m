function [f, grad, dir, h, h_inv] = HandlerVMulti(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3, which)

n = size(X1, 1);
d1 = size(X1, 2);
d2 = size(X2, 2);
d3 = size(X3, 2);
k = size(U, 1) / n;

regU = 1/(2*log(n));
regV1 = 1/(2*log(d1));
regV2 = 1/(2*log(d2));
regV3 = 1/(2*log(d3));

U = (reshape(U, k, n))';
V1 = (reshape(V1, k, d1))';
V2 = (reshape(V2, k, d2))';
V3 = (reshape(V3, k, d3))';

UV1 = U * V1' + ones(n, 1) * M1';
UV2 = U * V2' + ones(n, 1) * M2';
UV3 = U * V3' + ones(n, 1) * M3';

%obj function 
f = ... 
- trace(X1' * UV1) + ones(d1, 1)' * (1/2) * UV1' * UV1 * ones(d1, 1) + ...
- trace(X2' * UV2) + ones(n, 1)' * exp(UV2) * ones(d2, 1) + ...
- trace(X3' * UV3) + ones(n, 1)' * log(1 + exp(UV3)) * ones(d3, 1) + ...
c*a1*regU*1/2 * norm(U, 'fro')^2 + ... 
c*a2*regV1*1/2 * norm(V1, 'fro')^2 + ... 
c*a2*regV2*1/2 * norm(V2, 'fro')^2 + ... 
c*a2*regV3*1/2 * norm(V3, 'fro')^2; 
%+ c*a3*trace(U' * L * U);

%compute gradient
if nargout >= 2
	if which == 1
		grad = - X1' * U + UV1' * U + c*a2*regV1*V1;
		grad = reshape(grad', d1*k, 1);
	elseif which == 2
		grad = - X2' * U + exp(UV2)' * U + c*a2*regV2*V2;
		grad = reshape(grad', d2*k, 1);
	elseif which == 3
		grad = (-1) * X3' * U + sigmoid(UV3)' * U + c*a2*regV3*V3;
		grad = reshape(grad', d3*k, 1);
	end	
end

%%compute dir (-h_inv * grad)
%if nargout >= 3
%	UU = U' * U;
%	dir = zeros(size(grad));
%	for index=1:d
%		b = grad((index-1)*k + 1:index*k);
%		dir((index-1)*k + 1:index*k) = (-1) * UU \ b;
%	end	
%end
%
%%compute hessian
%if nargout >= 4
%	h = kron(eye(d), UU) + c*a2*regV*eye(d*k);
%end
%
%%compute inv hessian
%if nargout == 5
%	UU_inv = (UU + c*a2*regV*eye(k)) \ eye(size(UU));
%	h_inv =  kron(eye(d), UU_inv);
%end 

end
