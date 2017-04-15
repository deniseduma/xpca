%function [U, V1, M1, V2, M2, V3, M3] = exppca_multi(method, X1, X2, X3, ctypes, k, initU, initS1, initV1, initS2, initV2, initS3, initV3);
function [U, V1, M1, V2, M2, V3, M3] = main_multi(method, varargin)
    
    fprintf('Total number of input arguments = %d\n',nargin);
    
    %output mat file to save important vars
    outFile = strcat(method, '_', 'multi', '_', 'saved.mat');
    
    if nargin == 1 %simulations
    	n = 100; d1 = 60; d2= 60; d3 = 60; k = 3;
    	
	%make n a multiple of k
    	multiple = 1;
    	while (mod(n, multiple * k)~=0) 
    		n = n - 1;
    	end
    	
	%generate simulated data
	[X1, X2, X3, TrueU, TrueV1, TrueM1, TrueV2, TrueM2, TrueV3, TrueM3, initU, initS1, initV1, initM1, initS2, initV2, initM2, initS3, initV3, initM3] = generateDataMulti(outFile, n, d1, d2, d3, k, multiple);
    
    elseif nargin == 13 %real data
    	X1 = varargin{1};
    	X2 = varargin{2};
    	X3 = varargin{3};
        types = varargin{4};
	k = varargin{5};
    	initU = varargin{6};
    	initS1 = varargin{7};
    	initV1 = varargin{8};
    	initS2 = varargin{9};
    	initV2 = varargin{10};
    	initS3 = varargin{11};
    	initV3 = varargin{12};
	
	n = size(X1, 1);
	d1 = size(X1, 2);
	d2 = size(X2, 2);
	d3 = size(X3, 2);
	
	TrueU = inf; 
	TrueV1 = inf; TrueV2 = inf; TrueV3 = inf; 
	TrueM1 = inf; TrueM2 = inf; TrueM3 = inf;
	initM1 = zeros(d1, 1); initM2 = zeros(d2, 1); initM3 = zeros(d3, 1);
    
    end

    tol = 1e-3; 		
    c=1; a1=1; a2=1; a3=0;
    
    %adding smoothness constraints
    %A = diag(ones(n-1,1), 1) + diag(ones(n-1,1), -1);
    %D = diag(A * ones(n,1));
    %L = D - A; 
   
    %rearrange U and V into vector form
    initU = reshape(initU', n*k, 1);
    initV1 = reshape(initV1', d1*k, 1);
    initV2 = reshape(initV2', d2*k, 1);
    initV3 = reshape(initV3', d3*k, 1);
    
    %plot objective function value at each iteration
    figure('Name','Objective function','NumberTitle','off');
    crt_axis=gca;
    set(crt_axis,'XLim',[0 100]);
    set(crt_axis,'XTick',[0:10:100]);
    pause(0.03); hold on;

    %measure elapsed time
    tStart = tic;
    
    %call XPCA
    if strcmp(method, 'fminunc')
    	[U, V1, M1, V2, M2, V3, M3, iter] = xpca_multi(X1, X2, X3, initU, initV1, initM1, initV2, initM2, initV3, initM3, c, a1,a2, a3, tol, crt_axis);
    elseif strcmp(method, 'newton') 
    	[U, V1, M1, V2, M2, V3, M3, iter] = xpca_newton_multi(X1, X2, X3,  initU, initV1, initM1, initV2, initM2, initV3, initM3, c, a1, a2, a3, tol, crt_axis);
    end
    
    %measure elapsed time
    tElapsed = toc(tStart);

    %rearrange U and V back into original form
    U = (reshape(U, k, n))';
    V1 = (reshape(V1, k, d1))';
    V2 = (reshape(V2, k, d2))';
    V3 = (reshape(V3, k, d3))';

    assert(rank(U)==k, 'Rank of U is not k!');
    assert(rank(V1)==k, 'Rank of V1 is not k!');
    assert(rank(V2)==k, 'Rank of V2 is not k!');
    assert(rank(V3)==k, 'Rank of V3 is not k!');
   
   if nargin == 1
    	%plot output U
    	figure('Name','Output U','NumberTitle','off');
    	plot(U, 'LineWidth', 3); 
    	UN = normc(U);
    	figure('Name','Output U (normalized)','NumberTitle','off');
    	plot(UN, 'LineWidth', 3); 
    	pause(1);
    end	

    %compute performance measures for PCA and IPCA
    if nargin == 2
        %eval perf of HO GSVD
    	simPerformanceMulti(X1, X2, X3, TrueU, TrueV1, TrueM1, TrueV2, TrueM2, TrueV3, TrueM3, initU, initS1, initV1, zeros(d1, 1), initS2, initV2, zeros(d2, 1), initS3, initV3, zeros(d3, 1));
    	%eval perf of IPCA
	simPerformanceMulti(X1, X2, X3, TrueU, TrueV1, TrueM1, TrueV2, TrueM2, TrueV3, TrueM3, U, eye(k), V1, M1, eye(k), V2, M2, eye(k), V3, M3); 	
    end

        
    %number of iterations and elapsed time	
    fprintf('[exppca_multi]Number of iterations %d\n', iter);
    fprintf('[exppca_multi]Time elapsed %.2f\n', tElapsed);
    
    if nargin == 1
    	%save vars to mat file
    	save(outFile, 'U', 'V1', 'M1', 'V2', 'M2', 'V3', 'M3', 'iter', 'tElapsed', '-append');
   %elseif nargin == 13
   % 	save(outFile, 'initU', 'initS', 'initV', 'U', 'V', 'M', 'iter', 'tElapsed');
   end

end
