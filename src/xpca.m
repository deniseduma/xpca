function [U, V, M, ii] = xpca(dist, X, U, V, M, L, tol, maxit, crt_axis)

[n, d] = size(X);
k = length(U) / n;

%options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'HessUpdate', 'bfgs','CheckGradients',true,'FiniteDifferenceType','central');
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'HessUpdate', 'bfgs');
%options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);

iter=0;
gU = NaN;
gV = NaN;
gM = NaN;
for ii=1:maxit
    	
	prevU = U;
	prevV = V;
	prevM = M;
     	prev_gU = gU;
     	prev_gV = gV;
     	prev_gM = gM;

	if dist == 'n'
		handler1 = @(U)NormalU(X, U, V, M, L);
	elseif dist == 'p' 
		handler1 = @(U)PoissonU(X, U, V, M, L);
	elseif dist == 'b'
		handler1 = @(U)BernoulliU(X, U, V, M, L);
        end	
	%keep M and V fixed and update U
	[U, fval, eflag, out, gr] = fminunc(handler1, prevU, options);
	gU = norm(gr,'fro');% / (n*k);
	fprintf('\n[xpca:U]iter %d; #iters U %d; opt cond %.4f; norm(gr) %.4f\n', ii, out.iterations, out.firstorderopt, gU); 			
	
	%plot crt objective
	iter = iter +1; plotCrtObj(crt_axis, iter, fval); 

	if dist == 'n'
		handler2 = @(V)NormalV(X, U, V, M, L);
	elseif dist == 'p'	
		handler2 = @(V)PoissonV(X, U, V, M, L);
	elseif dist == 'b'
		handler2 = @(V)BernoulliV(X, U, V, M, L);
	end
	%keep U and M fixed and update V
	[V, fval, eflag, out, gr] = fminunc(handler2, prevV, options);
        gV = norm(gr,'fro');% / (d*k);
	fprintf('\n[xpca:V]iter %d; #iters V %d; opt cond %.4f; norm(gr) %.4f\n', ii, out.iterations, out.firstorderopt, gV); 			
	
	%plot crt objective
	iter = iter + 1; plotCrtObj(crt_axis, iter, fval); 
        
	if dist == 'n'
		handler3 = @(M)NormalM(X, U, V, M, L);
	elseif dist == 'p'
		handler3 = @(M)PoissonM(X, U, V, M, L);
	elseif dist == 'b'
		handler3 = @(M)BernoulliM(X, U, V, M, L);
	end
	%keep U and V fixed and update M
	[M, fval, eflag, out, gr] = fminunc(handler3, prevM, options);
        gM = norm(gr,'fro'); % / d;
	fprintf('\n[xpca:M]iter %d; #iters M %d; opt cond %.4f; norm(gr) %.4f\n', ii, out.iterations, out.firstorderopt, gM); 			
	
	%plot crt objective
	iter = iter + 1; plotCrtObj(crt_axis, iter, fval); 

	%check stopping criteria for U and V
        diff_U = norm(U - prevU, 'fro') / (n*k);
        diff_V = norm(V - prevV, 'fro') / (d*k);
        diff_M = norm(M - prevM, 'fro') / d;
        if diff_U<tol & diff_V<tol & diff_M<tol
       		fprintf('\t diff_U %.6f\n\n', diff_U);
       		fprintf('\t diff_V %.6f\n\n', diff_V);
       		fprintf('\t diff_M %.6f\n\n', diff_M);
       		break;
        end
	
	%if gU<=tol & gV<=tol & gM<=tol
       	%	fprintf('\t gU %.6f\n\n', gU);
       	%	fprintf('\t gV %.6f\n\n', gV);
       	%	fprintf('\t gM %.6f\n\n', gM);
       	%	break;
        %end

end
	iter = round(iter/3);
	
end
