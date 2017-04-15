function [U, V, M, iter] = xpca_newton(dist, X, U, V, M, L, c, a1, a2, a3, tol, crt_axis)

[n, d] = size(X);
k = size(U, 1) / n;

iter=0; citer=1;
while (true) 
    	
	prevU = U;
	prevV = V;
	prevM = M;
       	
	iter = iter + 1;
	
	%keep U and V fixed and update M
	if dist == 'n'
		handler3 = @(M)NormalM(X, U, V, M, L, c, a1, a2, a3);
	elseif dist == 'p'
		handler3 = @(M)PoissonM(X, U, V, M, L, c, a1, a2, a3);
	elseif dist == 'b'
		handler3 = @(M)BernoulliM(X, U, V, M, L, c, a1, a2, a3);
	end
	%solve for M iteratively
        [M, fval, niter]= newton(handler3, prevM, citer, crt_axis);
        citer = citer + niter + 1;
	
	%solve for M directly
        %if dist == 'n'
		%FIXME
	%	M = prevM;
	%elseif dist == 'p'
		%FIXME
	%	M = prevM;
	%elseif dist == 'b'
	%	[M, fval, niter]= newton(handler3, prevM, citer, crt_axis);
        %	citer = citer + niter + 1;
	%%	MX = X';
	%%	M=(log(MX ./ (1 - MX))-(reshape(V, k, d))'*reshape(U, k, n)) * ones(n, 1);
	%	M'
	%end
       
	[f, g1]=handler3(prevM);
	[f, g2]=handler3(M);
	gM1=norm(g1, 'fro') * 1/d;
	gM2=norm(g2, 'fro') * 1/d;
        fprintf('\n\n[xpca_newton:M]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gM1, gM2); 			
        
	%keep M and V fixed and update U
	if dist == 'n'
		handler1 = @(U)NormalU(X, U, V, M, L, c, a1, a2, a3);
	elseif dist == 'p' 
		handler1 = @(U)PoissonU(X, U, V, M, L, c, a1, a2, a3);
	elseif dist == 'b'
		handler1 = @(U)BernoulliU(X, U, V, M, L, c, a1, a2, a3);
        end	
	%solve for U
	[U, fval, niter] = newton(handler1, prevU, citer, crt_axis);
	citer = citer + niter + 1;
	
	[f, g1]=handler1(prevU);
	[f, g2]=handler1(U);
	gU1=norm(g1, 'fro') * 1/(n*k);
	gU2=norm(g2, 'fro') * 1/(n*k);
	fprintf('\n\n[xpca_newton:U]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gU1, gU2); 			
	%keep U and M fixed and update V
	if dist == 'n'
		handler2 = @(V)NormalV(X, U, V, M, L, c, a1, a2, a3);
	elseif dist == 'p'	
		handler2 = @(V)PoissonV(X, U, V, M, L, c, a1, a2, a3);
	elseif dist == 'b'
		handler2 = @(V)BernoulliV(X, U, V, M, L, c, a1, a2, a3);
	end
	%solve for V
	[V, fval, niter]= newton(handler2, prevV, citer, crt_axis);
	citer = citer + niter + 1;
                     	
	[f, g1]=handler2(prevV);
	[f, g2]=handler2(V);
	gV1=norm(g1, 'fro') * 1/(d*k);
	gV2=norm(g2, 'fro') * 1/(d*k);
        fprintf('\n\n[xpca_newton:V]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gV1, gV2); 			
        %check convergence criteria for U and V
        diff_U = norm(U - prevU, 'fro') * 1/(n*k);
        diff_V = norm(V - prevV, 'fro') * 1/(d*k);
        diff_M = norm(M - prevM, 'fro') * 1/d;
        if diff_U < tol & diff_V < tol & diff_M < tol
       		fprintf('\n\n[xpca_newton]diff_U %.6f\n', diff_U);
       		fprintf('[xpca_newton]diff_V %.6f\n', diff_V);
       		fprintf('[xpca_newton]diff_M %.6f\n\n', diff_M);
       		break;
        end

       %if gU2 < tol & gV2 < tol
       %	break;
       %end;
	
end
end
