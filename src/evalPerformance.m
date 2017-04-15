function [rss, srU, srV, rM, aU, paU, aV, paV] = evalPerformance(dist, X, TrueU, TrueV, TrueM, U, S, V, M) 

   [n, d] = size(X);
   k = size(U, 2);
   
   %compute Pearson RSS
    Theta = U * S *  V' + ones(n, 1) * M';    
    if dist == 'n'
    	Xmu = Theta;
    	%xmu(xmu<0)=0;
    elseif dist == 'p' 
    	Xmu = exp(Theta);
    	%xmu(xmu<0)=0;
    elseif dist == 'b' 
    	Xmu = sigmoid(Theta);
    	%xmu(xmu<0)=0; 
    	%xmu(xmu>1)=0.99;
    end
    %rss = norm(X - Xmu, 'fro') / std(std(xmu));
    rss = norm(X - Xmu, 'fro') * 1 / (n*d);
    fprintf('\n\n[simPerformance]Pearson RSS %.6f\n\n', rss);
	
    %compute subspace recovery measure
    UN = orth(U);
    VN = orth(V);
    srU = norm(UN * UN' - TrueU * TrueU', 'fro') * 1/(n*k);	
    srV = norm(VN * VN' - TrueV * TrueV', 'fro') * 1/(d*k);	
    rM = norm(M - TrueM, 'fro') * 1/d;	
    fprintf('\n\n[simPerformance]SR U %.4f\n', srU);
    fprintf('[simPerformance]SR V %.4f\n', srV);
    fprintf('[simPerformance]R M %.4f\n', rM);
    fprintf('[simPerformance]norm(TrueM) %.4f, norm(M) %.4f\n\n', norm(TrueM, 'fro'), norm(M, 'fro'));
    
    %display('U U');
    %U' * U
    %display('V V');
    %V' * V
    
    %display('UN UN');
    %UN' * UN
    %display('VN VN');
    %VN' * VN

    %compute principal angle for U
    %%[UN, R] = mygs(U);
    UU = UN' * TrueU;
    [US, S, VS] = svd(UU);
    paU = acos(S(end, end)) * 180/pi;
    aU = subspace(UN, TrueU); 
    fprintf('\n\n[simPerformance]Angle U %.6f\n', aU);
    fprintf('[simPerformance]Principal angle U %.6f\n\n', paU);
    
    %compute principal angle for V
    %[VN, R] = mygs(V);
    VV = VN' * TrueV;
    [US, S, VS] = svd(VV);
    paV = acos(S(end, end)) * 180/pi;
    aV = subspace(VN, TrueV); 
    fprintf('\n\n[simPerformance]Angle V %.6f\n', aV);
    fprintf('[simPerformance]Principal angle V %.6f\n\n', paV);
    
    %paM = subspace(M, TrueM); 
    %fprintf('[simPerformance]Principal angle M %.6f\n\n', paM);

    %display('UO UO');
    %UO' * UO
    %display('VO VO');
    %VO' * VO
end
