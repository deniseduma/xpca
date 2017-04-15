function netsmooth(inputdir)

[X,vname,cname]=tblread(['/home/duma/TCGA/data/', inputdir, '/', inputdir, '.txt'], ' ');
whos X
whos vname
whos cname
[x1,x2]=size(X);

afile = ['/home/duma/TCGA/data/', inputdir, '/', inputdir, '_adj.txt']
A = dlmread(afile);
whos A
[a1,a2]=size(A);

outFile = [inputdir, '_smooth.txt'];

alpha = 0.7;

iters=0;
oldF = X;
while(1) 
	newF = alpha*oldF*A + (1-alpha)*X;
	diff = norm(newF-oldF, 'fro');
	if (diff<1e-6) 
		disp(['final norm ', num2str(diff)]);
		break;
	end
	%disp(['iter ', num2str(iters), ', diff ', num2str(diff)]);
	oldF = newF;
	iters = iters + 1;
end
disp(['num iters to convergence ', num2str(iters)]);
whos newF

%quantile normalize the rows of newF
normF = quantilenorm(newF', 'display',1);
whos normF
F = normF';

%Save smoothed matrix to file
writeTXT(F, cname, char('ID', vname), inputdir, outFile);

end
