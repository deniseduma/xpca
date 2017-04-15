  function [X, U, S, V] = ov_wrapper(method, dist, algo, k)
  
  %% load mRNA data 
  inFile = ['../data/ov_mRNA_tum_sorted.txt']
  [X, vname, cname] = tblread(inFile, ' ');
  %%read cancer types
  auxFile = ['../data/ov_mRNA_types_types.txt']
  fileID = fopen(auxFile);
  ctypes = textscan(fileID, '%s');
  fclose(fileID);
  
  colors = ctypes{1, 1};
  cmap = {'1' [0 0 1]; '2' [0 1 1]; '3' [1 0 0]; '4' [0 1 0]; '5' [0 0 0]};
  [~, index] = ismember(colors, cmap(:,1));
  colors = cmap(index, 2);
  colors = cell2mat(colors);
 
  %%load mutation data
  %[X, vname, cname]= tblread('../data/bc_mut_cna_sorted_small.txt', ' ');
  %%read cancer types
  %fileID = fopen('../data/bc_mut_cna_types_types.txt');
  %ctypes = textscan(fileID, '%s');
  %fclose(fileID);
  
  %colors = ctypes{1, 1};
  %cmap = {'1' [0 0 1]; '2' [0 1 1]; '3' [1 0 0]; '4' [0 1 0]; '5' [0 0 0]};
  %[~, index] = ismember(colors, cmap(:,1));
  %colors = cmap(index, 2);
  %colors = cell2mat(colors);
  
  %% load methylation data
  %[X,vname,cname]= tblread('../data/bc_methylO_sorted_small.txt',' ');
  %%read cancer types
  %fileID = fopen('../data/bc_methylO_types_types.txt');
  %ctypes = textscan(fileID, '%s');
  %fclose(fileID);
  
  %colors = ctypes{1, 1};
  %cmap = {'1' [0 0 1]; '2' [0 1 1]; '3' [1 0 0]; '4' [0 1 0]; '5' [0 0 0]};
  %[~, index] = ismember(colors, cmap(:,1));
  %colors = cmap(index, 2);
  %colors = cell2mat(colors);
  
  whos X
  maxX = max(max(X))
  %whos vname
  %whos cname
  %whos colors
  
  [n, d] = size(X)
  
  if (strcmp(dist,'p'))
  	disp('(SVD)Poisson data');
	%for poisson data
	[initU, initS, initV] = svds(log(1 + X), k);
  else	
	[initU, initS, initV] = svds(X, k);
  end

  %return values
  U = initU; S = initS; V = initV;

  assert(rank(initU)==k, 'rank(initU) is not k!');
  assert(rank(initV)==k, 'rank(initV) is not k!');
  
  comps = [];
  for i=1:k
	comps = [comps; int2str(i)];
  end
  vname2 = strvcat('ID', comps);

  if (strcmp(method, 'svd')) %SVD
  	disp('***SVD method***');

  	initU = initU * initS;
	
	%plotScatter(initU, k, colors, 'SVD(X)')  
  	%initU_T = plotDendrogram(initU, k, 'initU');
  	
	%k-means clustering
	%kmeansopt = statset('kmeans');
	%kmeansopt.replicates = 10;
	%initU_T = kmeans(initU, k, 'options', kmeansopt);
	initU_T = kmeans(initU, k, 'Replicates', 1000);

	%write matrix to file
  	writeCSV(initU, cname, vname2, 'initU');
  	writeCSV(initU_T, cname, strvcat('ID', 'Cluster'), 'initU_T');
	      
   elseif (strcmp(method, 'mds')) %MDS
  	disp('***MDS method***');
   	
	if (strcmp(dist,'p'))
  		disp('(MDS)Poisson data');
   		%for Poisson data
   		D = pdist(log(1 + X), 'euclidean');
	else
		D = pdist(X, 'euclidean');
        end
  	[Y, eigvals] = cmdscale(D);
	
	%[ny, dy] = size(Y)
	%eigvals = eigvals(1:dy);
	%Y = Y * diag(eigvals);
	
	Y = Y(:, 1:k);
	%disp('Y * Y');
	%Y' * Y

	%plotScatter(Y, k, colors, 'MDS(X)');
   	%Y_T = plotDendrogram(Y, k, 'Y');
   
	%k-means clustering
	Y_T = kmeans(Y, k, 'Replicates', 1000);
   	
	%write matrix to file
   	writeCSV(Y, cname, vname2, 'Y');
   	writeCSV(Y_T, cname, strvcat('ID', 'Cluster'), 'Y_T');
   	
 elseif (strcmp(method, 'nmf')) %NMF
  	disp('***NMF method***');
	
	%read NMF factorization from file (done in R)
  	inFile2 = ['../data/ov_mRNA_tum_nmf.txt']
	[Z, vnameZ, cnameZ] = tblread(inFile2,' ');
	
	[nZ, dZ] = size(Z)
	
	%DEBUG
	for i=1:length(cname)
		assert(strcmp(cnameZ(i), cname(i)), 'Patient ids in different order!');
	end

	%plotScatter(Z, k, colors, 'NMF(X)')  
  	%Z_T = plotDendrogram(Z, k, 'Z');
  	
	Z_T = kmeans(Z, k, 'Replicates', 1000);

	%write matrix to file
  	writeCSV(Z, cname, vname2, 'Z');
  	writeCSV(Z_T, cname, strvcat('ID', 'Cluster'), 'Z_T');

   elseif (strcmp(method, 'xpca')) %XPCA
  	disp('***XPCA method***');
   	
  	%initU = initU * initS;
	
	%plot the histogram of X values
   	figure('Name','Histogram of X','NumberTitle','off');
   	hist(X, 20);
   	pause(1);

   	%call exppca
   	%first param is method, either 'fminunc' or 'newton'
   	[U, V, M] = exppca(algo, dist, X, ctypes, k, initU, initS, initV);

   	% norm U
	UN = normc(U);
	normColU = sqrt(diag(U' * U));
	D_U = diag(normColU);
	% norm V
	VN = normc(V);
	normColV = sqrt(diag(V' * V));
	D_V = diag(normColV);
	% new U
	U = UN * D_U * D_V;
	
	%principal components scatter plot	
   	%plotScatter(U, k, colors, 'XPCA(X)');
   	%U_T = plotDendrogram(U, k, 'U');
   
	%k-means clustering
	U_T = kmeans(U, k, 'Replicates', 1000);
   	
	%write matrix to file
   	writeCSV(U, cname, vname2, 'U');
   	writeCSV(U_T, cname, strvcat('ID', 'Cluster'), 'U_T');
   
   end
   	%write all matrices to file
   	%writeCSVs(X, initU, Y, U, vname, cname, types);
   
end


