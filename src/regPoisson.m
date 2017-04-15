function [regU, regV, regV2] = regPoisson(n, d, k)
	
	regU = 1;%1/sqrt(n);
	regV = 1;%1/sqrt(d);
	%regV2 = sqrt(d); %when svd(log(1 + X))
	regV2 = 1;%1/sqrt(d);
end
