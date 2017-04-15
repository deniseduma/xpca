function [X, TrueU, TrueV, TrueM, initU, initS, initV, initM] = generateData(dist, outFile, n, d, k, multiple)

    %generate U with sin and cos signals
    TrueV=[]; 
    B = linspace(-4*pi, 4*pi, d);
    s = -1;
    for iter=1:k
	%%amp_sin = 4*(1 - iter/10);
	%%amp_cos = 1.5*iter;
     	amp_sin = iter;
    	amp_cos = iter;
    	if mod(iter, 2) ~= 0
    		s = -s;
    		newCol = s*amp_sin*sin(1/iter*B);
    	   	TrueV = [TrueV newCol'];
    	else
    		newCol = amp_cos*cos(iter*B);
    		TrueV = [TrueV newCol'];
    	end	
   end

    %ortgonalize V
    region = d / (multiple * k);
    no_regions = d / region;
    for iter=1:no_regions
    	if (mod(iter,k)==0)
    		comp = k;
    	else 	
    		comp = mod(iter, k);
    	end
	%all components except for comp are set to zero in region iter
    	TrueV((1 + (iter-1) * region):iter * region, [1:(comp-1) (comp+1):k]) = zeros(region, k-1);
    end

    %generate U at random
    %TrueU = randn(n, k);
    
    %generate V with morlet 
    %TrueV = zeros(d,k);
    %TrueV(:,1)=morlet(-8,0,d);
    %TrueV(:,2)=gauswavf(0,8,d);
    %TrueV(:,3)=mexihat(-8,8,d);
    
    %generate V at random
    TrueU = randn(n, k);
    
    %generate V with heaviside function
    %TrueV=[]; 
    %BV = linspace(-2*pi, 2*pi, d);
    %s=-1;
    %for iter=1:k
    %	s = (-1) * s; 
    %	newColV=1/iter * heaviside(s*BV);
    %	TrueV = [TrueV newColV'];
    %end
    
    %generate V with step function
    %TrueV=zeros(d, k);
    %for iter=1:k
    %	s = 1;
    %	step = 5;
    %	amp = 1 - (iter / 10);
    %	for iter1=1:step:floor(d / step) * step
    %		s = -s;
    %		for iter2=iter1:iter1 + step - 1
    % 			TrueV(iter2, iter) = amp * s;
    %		end	
    %	end
    %	s = -s;
    %	for iter1=(floor(d / step) * step):d
    %	   	TrueV(iter2, iter) = amp * s;
    %	end
    %end	

    display(['size of TrueU ' num2str(size(TrueU)) ' rank ' num2str(rank(TrueU))]);
    assert(rank(TrueU)==k, 'Rank of U is not k!');
    display(['size of TrueV ' num2str(size(TrueV)) ' rank ' num2str(rank(TrueV))]);
    assert(rank(TrueV)==k, 'Rank of V is not k!');
    
    disp('TrueV before col normalization');
    TrueV;	
    TrueV'*TrueV
    
    TrueU=orth(TrueU);
    TrueV = normc(TrueV);	
    
    disp('TrueV after col normalization');
    TrueV;	     
    TrueV'*TrueV
    
    disp('TrueU after orthonormalization');
    TrueU'*TrueU	    
    
    %plot true U and V
    %figure('Name','True U','NumberTitle','off');
    %plot(TrueU,'LineWidth', 3);
    figure('Name','True V','NumberTitle','off');
    plot(TrueV,'LineWidth', 3);
    pause(1);
    
    %generate intercept
    TrueM = zeros(d, 1);
    %TrueM(1 : round(d/5)) = -1;
    %TrueM(round(d/5) + 1 : round(2*d/5)) = -0.5;
    %TrueM(round(2*d/5) + 1 : round(3*d/5)) = 0;
    %TrueM(round(3*d/5) + 1 : round(4*d/5)) = 0.5;
    %TrueM(round(4*d/5) + 1 : end) = 1;
    
    %scaling matrix D
    m = min(n, d);
    lm = log(m);
    if dist == 'n'
    	D = diag([ 2*m 1.5*m 1*m])
    elseif dist == 'p'
    	D = diag([ 1*m 0.75*m 0.5*m])
    elseif dist == 'b'
    	D = diag([ 3*m 2*m 1*m])
    end
    D = diag([ 1*m 1*m 1*m])

    %compute theta 
    O = TrueU * D * TrueV' + ones(n, 1) * TrueM';
    
    %generate X according to various exponential family distributions 
    if dist == 'n' %Normal
    	disp('Normal dist');
	mu = O;
    	X = mu + randn(n, d);
    elseif dist == 'p' %Poisson
    	disp('Poisson dist');
    	mu = exp(O);
    	X = poissrnd(mu);
    elseif dist == 'b' %Bernoulli
    	disp('Bernoulli dist');
    	mu = sigmoid(O);
    	X = binornd(1, mu);
    end

    figure('Name','Distribution of X','NumberTitle','off');
    hist(X, 20);
    pause(1);
    %if dist == 'p'
    %	figure('Name','Distribution of log(X)','NumberTitle','off');
    %	hist(log(1 + X), 20);
    %	pause(1);
    %end
    figure('Name','Distribution of O','NumberTitle','off');
    hist(O, 20);
    pause(1);
    figure('Name','g(O) vs O','NumberTitle','off');
    plot(O, mu);
    pause(1);
    
    MO=mean(mean(O))
    MX=mean(mean(X))

    %init U and V by SVD
    %if dist == 'p'
    %	[initU, initS, initV] = svds(log(1 + X), k);
    %else 
    	[initU, initS, initV] = svds(X, k);
    %end
    assert(rank(initU)==k);
    assert(rank(initV)==k);
    %init M
    initM = zeros(d, 1);
   
    %plot init U and V
    %figure('Name','Init U','NumberTitle','off');
    %plot(initU,'LineWidth', 3);
    figure('Name','Init V','NumberTitle','off');
    plot(initV,'LineWidth', 3);
    pause(1);

    %save vars to mat file
    save(outFile, 'X', 'TrueU', 'D', 'TrueV', 'TrueM', 'initU', 'initS', 'initV', 'initM'); 	
end
