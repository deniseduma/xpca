function [X1, X2, X3, TrueU, TrueV1, TrueM1, TrueV2, TrueM2, TrueV3, TrueM3, initU, initS1, initV1, initM1, initS2, initV2, initM2, initS3, initV3, initM3] = generateDataMulti(outFile, n, d1, d2, d3, k, multiple)

    %make n a multiple of k
    %while (mod(n, multiple * k)~=0) 
    %	n = n - 1;
    %end

    %generate U with sin and cos signals
    TrueU=[]; 
    B = linspace(-4*pi, 4*pi, n);
    s = -1;
    for iter=1:k
	%%amp_sin = 4*(1 - iter/10);
	%%amp_cos = 1.5*iter;
     	amp_sin = iter;
    	amp_cos = iter;
    	if mod(iter, 2) ~= 0
    		s = -s;
    		newCol = s*amp_sin*sin(1/iter*B);
    	   	TrueU = [TrueU newCol'];
    	else
    		newCol = amp_cos*cos(iter*B);
    		TrueU = [TrueU newCol'];
    	end	
   end

    %ortgonalize U
    region = n / (multiple * k);
    no_regions = n / region;
    for iter=1:no_regions
    	if (mod(iter,k)==0)
    		comp = k;
    	else 	
    		comp = mod(iter, k);
    	end
	%all components except for comp are set to zero in region iter
    	TrueU((1 + (iter - 1) * region):iter * region, [1:(comp-1) (comp+1):k]) = zeros(region, k-1);
    end

    %generate V  
    TrueV1 = zeros(d1, k);
    TrueV2 = zeros(d2, k);
    TrueV3 = zeros(d3, k);

    TrueV1(:,1)=morlet(-8,0,d1);
    TrueV1(:,2)=morlet(-16,0,d1);
    TrueV1(:,3)=morlet(-24,0,d1);
    
    TrueV2(:,1)=mexihat(-8,8,d2);
    TrueV2(:,2)=mexihat(-16,16,d2);
    TrueV2(:,3)=mexihat(-24,24,d2);
    
    TrueV3(:,1)=gauswavf(0,8,d3);
    TrueV3(:,2)=gauswavf(0,16,d3);
    TrueV3(:,3)=gauswavf(0,24,d3);

    
    display(['size of TrueU ' num2str(size(TrueU)) ' rank ' num2str(rank(TrueU))]);
    assert(rank(TrueU)==k, 'Rank of U is not k!');
    %disp('TrueU * TrueU');
    %TrueU' * TrueU	     
    
    TrueU=normc(TrueU);
    TrueV1 = orth(TrueV1);	
    TrueV2 = orth(TrueV2);	
    TrueV3 = orth(TrueV3);	

    display(['size of TrueV1 ' num2str(size(TrueV1)) ' rank ' num2str(rank(TrueV1))]);
    assert(rank(TrueV1)==k, 'Rank of V1 is not k!');
    %disp('TrueV1 * TrueV1');
    %TrueV1' * TrueV1	     
    
    display(['size of TrueV2 ' num2str(size(TrueV2)) ' rank ' num2str(rank(TrueV2))]);
    assert(rank(TrueV2)==k, 'Rank of V2 is not k!');
    %disp('TrueV2 * TrueV2');
    %TrueV2' * TrueV2	     
    
    display(['size of TrueV3 ' num2str(size(TrueV3)) ' rank ' num2str(rank(TrueV3))]);
    assert(rank(TrueV3)==k, 'Rank of V3 is not k!');
    %disp('TrueV3 * TrueV3');
    %TrueV3' * TrueV3	     
    
    %DEBUG
    disp(['rank(TrueU)=' num2str(rank(TrueU))]);
    disp(['rank(TrueV1)=' num2str(rank(TrueV1))]);
    disp(['rank(TrueV2)=' num2str(rank(TrueV2))]);
    disp(['rank(TrueV3)=' num2str(rank(TrueV3))]);
    fprintf('\n');
    
    %plot true U and V
    figure('Name','True U','NumberTitle','off');
    plot(TrueU,'LineWidth', 3);
    %V1
    %figure('Name','True V1','NumberTitle','off');
    %plot(TrueV1,'LineWidth', 3);
    %V2
    %figure('Name','True V2','NumberTitle','off');
    %plot(TrueV2,'LineWidth', 3);
    %V3
    %figure('Name','True V3','NumberTitle','off');
    %plot(TrueV3,'LineWidth', 3);
    pause(1);
    
    %generate intercepts
    TrueM1 = zeros(d1, 1);
    TrueM2 = zeros(d2, 1);
    TrueM3 = zeros(d3, 1);
    
    %TrueM1
    TrueM1(1 : round(d1/5)) = 0.3;
    TrueM1(round(d1/5) + 1 : round(2*d1/5)) = 0.7;
    TrueM1(round(2*d1/5) + 1 : round(3*d1/5)) = 0.2;
    TrueM1(round(3*d1/5) + 1 : round(4*d1/5)) = 0.9;
    TrueM1(round(4*d1/5) + 1 : end) = 0.5;
    %TrueM1 = 5 * rand(d1, 1);
    
    %TrueM2
    TrueM2(1 : round(d2/5)) = -2;
    TrueM2(round(d2/5) + 1 : round(2*d2/5)) = -1;
    TrueM2(round(2*d2/5) + 1 : round(3*d2/5)) = 0;
    TrueM2(round(3*d2/5) + 1 : round(4*d2/5)) = 1;
    TrueM2(round(4*d2/5) + 1 : end) = 2;
    %TrueM2 = 10 * rand(d2, 1);
    
    %TrueM3
    TrueM3(1 : round(d3/5)) = -7;
    TrueM3(round(d3/5) + 1 : round(2*d3/5)) = -5;
    TrueM3(round(2*d3/5) + 1 : round(3*d3/5)) = 0;
    TrueM3(round(3*d3/5) + 1 : round(4*d3/5)) = 5;
    TrueM3(round(4*d3/5) + 1 : end) = 7;
    %TrueM3 = 15 * rand(d3, 1);%DEBUG
    
    %scaling matrix D
    m1 = min(n, d1);
    m2 = min(n, d2);
    m3 = min(n, d3);
    %Poisson: C=0.5, D = diag([ 1/2m 1/3m 1/4m]), a1=a2=0.2
    %Poisson break SVD: D = diag([ 1m 0.75m 0.5m]), a1=a2=0.5
    %Poisson break SVD: C=0.5, D = diag([ 2m 1.5m 1m]), a1=a2=0.2
    %if dist == 'n'
    	D1 = diag([ 2*m1 1.5*m1 1*m1]);
    %elseif dist == 'p'
    	D2 = diag([ 1*m2 0.75*m2 0.5*m2]);
    %elseif dist == 'b'
    	D3 = diag([ 3*m3 2*m3 1*m3]);
    %end

    %compute theta params
    Theta1 = TrueU * D1 * TrueV1';
    Theta2 = TrueU * D2 * TrueV2';
    Theta3 = TrueU * D3 * TrueV3';
    
    Theta1 = Theta1 + ones(n, 1) * TrueM1';
    Theta2 = Theta2 + ones(n, 1) * TrueM2';
    Theta3 = Theta3 + ones(n, 1) * TrueM3';
    
    %DEBUG
    disp(['rank(Theta1)=' num2str(rank(Theta1))]);
    disp(['rank(Theta2)=' num2str(rank(Theta2))]);
    disp(['rank(Theta3)=' num2str(rank(Theta3))]);
    fprintf('\n');
    
    %generate X's according to various exponential family distributions 
    %if dist == 'n' %Normal
    	mu1 = Theta1;
    	X1 = mu1 + randn(n, d1);
    %elseif dist == 'p' %Poisson
    	mu2 = exp(Theta2);
    	X2 = poissrnd(mu2);
    %elseif dist == 'b' %Bernoulli
    	mu3 = sigmoid(Theta3);
    	X3 = binornd(1, mu3);
    %end
    
    %DEBUG
    disp(['rank(X1)=' num2str(rank(X1))]);
    disp(['rank(X2)=' num2str(rank(X2))]);
    disp(['rank(X3)=' num2str(rank(X3))]);
    %assert(rank(X1)==d1, 'Rank of X1 is not d1!');
    %assert(rank(X2)==d2, 'Rank of X2 is not d2!');
    %assert(rank(X3)==d3, 'Rank of X3 is not d3!');
    fprintf('\n');

    %figure('Name','Distribution of X','NumberTitle','off');
    %hist(X, 20);
    %pause(1);
    %if dist == 'p'
    %	figure('Name','Distribution of log(X)','NumberTitle','off');
    %	hist(log(1 + X), 20);
    %	pause(1);
    %end
    
    %figure('Name','Distribution of Theta','NumberTitle','off');
    %hist(Theta, 20);
    %pause(1);
    figure('Name','g(Theta1) vs Theta1','NumberTitle','off');
    plot(Theta1, mu1);
    figure('Name','g(Theta2) vs Theta2','NumberTitle','off');
    plot(Theta2, mu2);
    figure('Name','g(Theta3) vs Theta3','NumberTitle','off');
    plot(Theta3, mu3);
    pause(1);
    
    MTheta1=mean(mean(Theta1))
    MX1=mean(mean(X1))
    MTheta2=mean(mean(Theta2))
    MX2=mean(mean(X2))
    MTheta3=mean(mean(Theta3))
    MX3=mean(mean(X3))

    %init U and V by HO GSVD
    %[initU, initS, initV] = svds(X, k);
    %[U, S, V] = svds(log(1 + X), k); %HO GSVD init
    
    %init U, S and V by HO GSVD
    X2 = log(1 + X2);
   
    A1 = X1 * X1'; 
    A2 = X2 * X2'; 
    A3 = X3 * X3';
        
    I1 = inv(A1);
    I2 = inv(A2);
    I3 = inv(A3);
    
    %DEBUG
    disp(['rank(A1)=' num2str(rank(A1))]);
    disp(['rank(A2)=' num2str(rank(A2))]);
    disp(['rank(A3)=' num2str(rank(A3))]);
    fprintf('\n');

    disp(['rank(I1)=' num2str(rank(I1))]);
    disp(['rank(I2)=' num2str(rank(I2))]);
    disp(['rank(I3)=' num2str(rank(I3))]);
    fprintf('\n');
   
    S = 1/6 * (A1 * I2 + A2 * I1 + A1 * I3 + A3 * I1 + A2 * I3 + A3 * I2);
    [initU, D] = eig(S);

    initU = initU(:, 1:k); 
   
    B1 = X1' * initU;
    B2 = X2' * initU;
    B3 = X3' * initU;

    initS1 = diag(sqrt(diag(B1' * B1)));
    initS2 = diag(sqrt(diag(B2' * B2)));
    initS3 = diag(sqrt(diag(B3' * B3)));
   
    initV1 = B1 * diag(1 ./ diag(initS1));
    initV2 = B2 * diag(1 ./ diag(initS2));
    initV3 = B3 * diag(1 ./ diag(initS3));

    %initS1 = eye(k); initS2 = eye(k); initS3 = eye(k);
    %initV1 = B1; initV2 = B2; initV3 = B3;
    
    %DEBUG
    disp(['rank(initU)=' num2str(rank(initU))]);
    disp(['rank(initV1)=' num2str(rank(initV1))]);
    disp(['rank(initV2)=' num2str(rank(initV2))]);
    disp(['rank(initV3)=' num2str(rank(initV3))]);
    fprintf('\n');

    %init intercepts
    initM1 = zeros(d1, 1);
    initM2 = zeros(d2, 1);
    initM3 = zeros(d3, 1);
   
    %plot init U and V
    figure('Name','Init U','NumberTitle','off');
    plot(initU,'LineWidth', 3);
    %figure('Name','Init V','NumberTitle','off');
    %plot(V,'LineWidth', 3);
    pause(1);

    %save vars to mat file
    save(outFile, 'X1', 'X2', 'X3', 'TrueU', 'D1', 'TrueV1', 'TrueM1', 'D2', 'TrueV2', 'TrueM2', 'D3', 'TrueV3', 'TrueM3', 'initU', 'initS1', 'initV1', 'initM1', 'initS2', 'initV2', 'initM2', 'initS3', 'initV3', 'initM3'); 	


end
