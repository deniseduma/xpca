function [f, O] = objFunc(X, U, V, M, L, regU, regV, regV2) 
	
	[n, d] = size(X);
	k = size(U, 2);
	O = U*V' + ones(n, 1)*M';
	%f = -trace(X'*O) + 0.5*regU*norm(U,'fro')^2 + 0.5*regV*norm(V,'fro')^2 + 0.5*regV2*trace(V'*L*V);
	f = -trace(X'*O) + 0.5*regU*norm(U,'fro')^2 + 0.5*regV2*trace(V'*L*V);

end
