%[U, V, M] = function exppca(algo, dist, X, types, k, initU, initS, initV)
function [U, V, M] = main(algo, dist, indir, varargin)
 	
	tol = 1e-4;
	maxit = 200;
   
	fprintf('Total number of input arguments = %d\n',nargin);
    
	%output mat file to save important vars
	outFile = strcat(algo, '_', dist, '_saved.mat');
    
	if nargin == 3
    		n = 100; d = 60; k = 3;
    		%make d a multiple of k
    		multiple = 1;
    		while (mod(d, multiple * k)~=0) 
    			d = d - 1;
    		end
    		%generate simulated data
		[X, TrueU, TrueV, TrueM, initU, initS, initV, initM] = generateData(dist, outFile, n, d, k, multiple);
    
	elseif nargin == 9
    		X = varargin{1};
		[n, d] = size(X);
		k = varargin{2};
    		initU = varargin{3};
    		initS = varargin{4};
    		initV = varargin{5};
    		L = varargin{6};
		initM = zeros(d, 1);
		TrueU = inf; 
		TrueV = inf; 
		TrueM = inf; 
	end

    %plot objective function value at each iteration
    figure('Name','Objective function','NumberTitle','off');
    crt_axis=gca;
    set(crt_axis,'XLim',[0 200]);
    set(crt_axis,'XTick',[0:20:200]);
    pause(0.03); hold on;
    
    %measure elapsed time
    tStart = tic;
    
    %call XPCA
    if strcmp(algo, 'fminunc')
    	[U, V, M, iter] = xpca(dist, X, initU(:), initV(:), initM, L, tol, maxit, crt_axis);
	U = reshape(U,n,k);	
    	V = reshape(V,d,k);
    elseif strcmp(algo, 'newton') 
    	[U, V, M, iter] = xpca_newton(dist, X, initU, initV, initM, L, tol, maxit, crt_axis);
    end
    
    %measure elapsed time
    tElapsed = toc(tStart);
    
    assert(rank(U)==k, 'Rank of U is not k!');
    assert(rank(V)==k, 'Rank of V is not k!');
   
    if nargin == 3
    	%plot output V
    	VN = normc(V);
    	figure('Name','Output V','NumberTitle','off');
    	plot(VN, 'LineWidth', 3); 
    	pause(1);
    end	
    
    %compute performance measures for PCA and XPCA
    if nargin == 3
    	disp('**SVD**');
	evalPerformance('n', X, TrueU, TrueV, TrueM, initU, initS, initV, zeros(d, 1));
    	disp('**XPCA**')
	evalPerformance(dist, X, TrueU, TrueV, TrueM, U, eye(k), V, M); 	
    end
    
    %number of iterations and elapsed time	
    fprintf('[main]Number of iterations %d\n', iter);
    fprintf('[main]Elapsed time %.2f\n', tElapsed);
    
    %save vars to mat file
    if nargin == 3
    	    	save(outFile, 'U', 'V', 'M', 'iter', 'tElapsed', '-append');
   elseif nargin == 9
    	save(outFile, 'initU', 'initS', 'initV', 'U', 'V', 'M', 'iter', 'tElapsed');
   end

end
