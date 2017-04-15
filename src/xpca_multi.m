function [U, V1, M1, V2, M2, V3, M3, iter] = xpca_multi(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3, tol, crt_axis)

n = size(X1, 1);
d1 = size(X1, 2);
d2 = size(X2, 2);
d3 = size(X3, 2);
k = size(U, 1) / n;

options = optimset('LargeScale', 'off' ,'HessUpdate', 'bfgs','GradObj', 'on');
%options = optimoptions(@fmincon,'Algorithm','quasi-newton', 'HessUpdate', 'bfgs', 'GradObj', 'on')

iter=0;
while (true) 
    	
	prevU = U;
	prevV1 = V1;
	prevM1 = M1;
	prevV2 = V2;
	prevM2 = M2;
	prevV3 = V3;
	prevM3 = M3;
       	
	iter = iter + 1;
	
	%keep Us and Vs fixed and update M1
	handler3 = @(M1)HandlerMMulti(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3, 1);
	%solve for M
	[M1, fval] = fminunc(handler3, prevM1, options);
        %plot
        plotCrtObj(crt_axis, iter, fval); 
       
	[f, g1]=handler3(prevM1);
	[f, g2]=handler3(M1);
	gM1=norm(g1, 'fro') * 1/(d1);
	gM2=norm(g2, 'fro') * 1/(d1);
        fprintf('\n\n[xpca:M1]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gM1, gM2); 			
	iter = iter + 1;
	
	%keep U and M fixed and update V1
	handler2 = @(V1)HandlerVMulti(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3, 1);
	%solve for V
	[V1, fval] = fminunc(handler2, prevV1, options);
        %plot
        plotCrtObj(crt_axis, iter, fval); 
                     	
	[f, g1]=handler2(prevV1);
	[f, g2]=handler2(V1);
	gV1=norm(g1, 'fro') * 1/(d1*k);
	gV2=norm(g2, 'fro') * 1/(d1*k);
        fprintf('\n\n[xpca:V1]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gV1, gV2);

	iter = iter + 1;
	
	%keep Us and Vs fixed and update M2
	handler3 = @(M2)HandlerMMulti(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3, 2);
	%solve for M
	[M2, fval] = fminunc(handler3, prevM2, options);
        %plot
        plotCrtObj(crt_axis, iter, fval); 
       
	[f, g1]=handler3(prevM2);
	[f, g2]=handler3(M2);
	gM1=norm(g1, 'fro') * 1/(d2);
	gM2=norm(g2, 'fro') * 1/(d2);
        fprintf('\n\n[xpca:M2]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gM1, gM2); 			
	iter = iter + 1;

	%keep U and M fixed and update V2
	handler2 = @(V2)HandlerVMulti(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3, 2);
	%solve for V
	[V2, fval] = fminunc(handler2, prevV2, options);
        %plot
        plotCrtObj(crt_axis, iter, fval); 
                     	
	[f, g1]=handler2(prevV2);
	[f, g2]=handler2(V2);
	gV1=norm(g1, 'fro') * 1/(d2*k);
	gV2=norm(g2, 'fro') * 1/(d2*k);
        fprintf('\n\n[xpca:V2]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gV1, gV2); 			
	iter = iter + 1;
	
	%keep U and V fixed and update M3
	handler3 = @(M3)HandlerMMulti(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3, 3);
	%solve for M
	[M3, fval] = fminunc(handler3, prevM3, options);
        %plot
        plotCrtObj(crt_axis, iter, fval); 
       
	[f, g1]=handler3(prevM3);
	[f, g2]=handler3(M3);
	gM1=norm(g1, 'fro') * 1/(d3);
	gM2=norm(g2, 'fro') * 1/(d3);
        fprintf('\n\n[xpca:M3]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gM1, gM2); 			
	iter = iter + 1;	
	
	%keep U and M fixed and update V3
	handler2 = @(V3)HandlerVMulti(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3, 3);
	%solve for V
	[V3, fval] = fminunc(handler2, prevV3, options);
        %plot
        plotCrtObj(crt_axis, iter, fval); 
                     	
	[f, g1]=handler2(prevV3);
	[f, g2]=handler2(V3);
	gV1=norm(g1, 'fro') * 1/(d3*k);
	gV2=norm(g2, 'fro') * 1/(d3*k);
        fprintf('\n\n[xpca:V3]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gV1, gV2); 		
	
	iter = iter + 1;
        
	%keep M and V fixed and update U
	handler1 = @(U)HandlerUMulti(X1, X2, X3, U, V1, M1, V2, M2, V3, M3, c, a1, a2, a3);
	%solve for U
	[U, fval] = fminunc(handler1, prevU, options);
	%plot
	plotCrtObj(crt_axis, iter, fval); 
                     	
	[f, g1]=handler1(prevU);
	[f, g2]=handler1(U);
	gU1=norm(g1, 'fro') * 1/(n*k);
	gU2=norm(g2, 'fro') * 1/(n*k);
	fprintf('\n\n[xpca:U]iter %d; g1 %.6f; g2 %.6f\n\n', iter, gU1, gU2); 			
        %check convergence criteria for U and V
        diffU = norm(U - prevU, 'fro') * 1/(n*k);
        diffV1 = norm(V1 - prevV1, 'fro') * 1/(d1*k);
        diffV2 = norm(V2 - prevV2, 'fro') * 1/(d2*k);
        diffV3 = norm(V3 - prevV3, 'fro') * 1/(d3*k);
        diffM1 = norm(M1 - prevM1, 'fro') * 1/d1;
        diffM2 = norm(M2 - prevM2, 'fro') * 1/d2;
        diffM3 = norm(M3 - prevM3, 'fro') * 1/d3;
        if diffU < tol & diffV1 < tol & diffV2 < tol & diffV3 < tol & diffM1 < tol & diffM2 < tol & diffM3 < tol
       		fprintf('\t diffU %.4f\n\n', diffU);
       		fprintf('\t diffV1 %.4f\n\n', diffV1);
       		fprintf('\t diffV2 %.4f\n\n', diffV2);
       		fprintf('\t diffV3 %.4f\n\n', diffV3);
       		fprintf('\t diffM1 %.4f\n\n', diffM1);
       		fprintf('\t diffM2 %.4f\n\n', diffM2);
       		fprintf('\t diffM3 %.4f\n\n', diffM3);
       		break;
        end

end
	iter = round(iter/7);
	
end
