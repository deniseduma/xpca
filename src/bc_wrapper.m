function [X, U, S, V] = bc_wrapper(method, dist, algo, inputdir, k)
  
  classifier = 'K'; 
  if (strcmp(classifier, 'C'))
  	distance = 'complete';
  elseif (strcmp(classifier, 'W'))
  	distance = 'ward';
  end

  %Read in survival data
  if (strcmp(inputdir, 'ucec_mrna') | strcmp(inputdir, 'ucec_mut') | strcmp(inputdir, 'ucec_met')) 
  	srvFile = ['../data/' , inputdir , '/', inputdir, '_histology.txt']
  else 
  	srvFile = ['../data/' , inputdir , '/', inputdir, '_survival.txt']
  end
  [SRV,vname_srv,cname_srv] = tblread(srvFile, ' ');
  whos SRV
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%
  
  if (strcmp(inputdir, 'bc_mrna') | strcmp(inputdir, 'luad_mrna') | strcmp(inputdir, 'ov_mrna') | strcmp(inputdir, 'ucec_mrna')) 
  	dta = 'mrna';
  elseif (strcmp(inputdir, 'bc_mut') | strcmp(inputdir, 'luad_mut') | strcmp(inputdir, 'ov_mut') | strcmp(inputdir, 'ucec_mut')) 
  	dta = 'mut';
  elseif (strcmp(inputdir, 'bc_met') | strcmp(inputdir, 'luad_met') | strcmp(inputdir, 'ov_met') | strcmp(inputdir, 'ucec_met')) 
  	dta = 'met';
  end 
  
  % Read input data
  inFile = ['../data/' , inputdir , '/', inputdir, '.txt']
  [X,vname,cname] = tblread(inFile, ' ');
  if (strcmp(dist,'p') | strcmp(dist,'n'))
  	X = log10(X+1);
  end
  [n,d] = size(X);
  % Center columns of X
  X = X - ones(n,1)*mean(X,1);
  whos X
  disp(['rank(X) ', num2str(rank(X))]);  

  %adding smoothness constraints
  % DUMMY Laplacian
  %A = diag(ones(d-1,1), 1) + diag(ones(d-1,1), -1);
  %D = diag(A * ones(d,1));
  %L = D - A;

  %L = laplacian('STRING.txt', indir);
  L = dlmread(['../data/' inputdir '/' inputdir '_laplacian.txt'], '\t');

  %%read cancer types
  typesFile = ['../data/' , inputdir , '/', inputdir, '_types.txt']
  fileID = fopen(typesFile);
  ctypes = textscan(fileID, '%s');
  fclose(fileID);
  ctypes = ctypes{1, 1};
  %%blue cyan yellow green black
  %%cmap = {'LuminalA' 'blue'; 'LuminalB' 'cyan'; 'HER2' 'red'; 'Basal' 'green'; 'Normal' 'black'}
  cmap = {'LuminalA' [0 0 1]; 'LuminalB' [0 1 1]; 'HER2' [1 0 0]; 'Basal' [0 1 0]; 'Normal' [0 0 0]};
  [~, index] = ismember(ctypes, cmap(:,1));
  ctypes = cmap(index, 2);
  ctypes = cell2mat(ctypes);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load mutation data %%%%%%%%%%%%%%%%%%%
  %[X, vname, cname]= tblread('../data/bc_mut_cna_sorted_small.txt', ' ');
  %%read cancer types
  %fileID = fopen('../data/bc_mut_cna_types_types.txt');
  %ctypes = textscan(fileID, '%s');
  %fclose(fileID);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load methylation data %%%%%%%%%%%%%%%%%%%
  %[X,vname,cname]= tblread('../data/bc_methylO_sorted_small.txt',' ');
  %%read cancer types
  %fileID = fopen('../data/bc_methylO_types_types.txt');
  %ctypes = textscan(fileID, '%s');
  %fclose(fileID);
  %ctypes = ctypes{1, 1};
  %cmap = {'1' [0 0 1]; '2' [0 1 1]; '3' [1 0 0]; '4' [0 1 0]; '5' [0 0 0]};
  %[~, index] = ismember(ctypes, cmap(:,1));
  %ctypes = cmap(index, 2);
  %ctypes = cell2mat(ctypes);
    
  % Intersect data and survival observations 
  interXSRV = intersect(cellstr(cname_srv), cellstr(cname));   
  SRV = SRV( ismember( cellstr(cname_srv), interXSRV ), : );
  X = X( ismember( cellstr(cname), interXSRV ), : );
  cname = char(interXSRV);
  disp('inter(X, SRV)')
  whos X
  whos SRV
  
  % Take SVD of X to init algo 
  [initU, initS, initV] = svds(X, k);
  U = initU; S = initS; V = initV;

  comps = [];
  for i=1:k
	comps = [comps; 'PC', int2str(i)];
  end
  vname2 = strvcat('SampleID', comps);

  if (strcmp(method, 'svd')) %SVD
  	disp('***SVD method***');

  	initU = initU * initS;
	
  	if (strcmp(classifier, 'C') | strcmp(classifier, 'W')) 
		%hierarchical clustering
		initU_T = plotDendrogram(initU, k, ['initU_', classifier, '_k', num2str(k)], distance);
	elseif (strcmp(classifier, 'K')) 
		%k-means clustering
		initU_T = kmeans(initU, k, 'Replicates', 1000);
	end
	initU_T = [SRV initU_T];
		
	%write matrix to file
  	writeCSV(initU, cname, vname2, [inputdir,'/',method], ['U', '_k', num2str(k), '.csv']);
  	%write survival data and cluster memberships to file
	if (strcmp(inputdir, 'ucec_mrna') | strcmp(inputdir, 'ucec_mut') | strcmp(inputdir, 'ucec_met')) 
  		writeCSV(initU_T, cname, strvcat('ID', 'Histology', 'Cluster'), [inputdir,'/',method], ['U_', classifier, '_k', num2str(k), '.csv']);
	else 
  		writeCSV(initU_T, cname, strvcat('ID', 'OS_Event', 'OS_Time', 'Cluster'), [inputdir,'/',method], ['U_', classifier, '_k', num2str(k), '.csv']);
	end
	plotScatter(initU, k, ctypes, method, ['../data/',inputdir,'/',method,'/']);
   
   elseif (strcmp(method, 'mds')) %MDS
  	disp('***MDS method***');
   	
	D = pdist(X, 'euclidean');
  	[Y, eigvals] = cmdscale(D);
	Y = Y(:, 1:k);
	whos Y

	if (strcmp(classifier, 'C') | strcmp(classifier, 'W')) 
		%hierarchical clustering
		Y_T = plotDendrogram(Y, k, ['Y_', classifier, '_k', num2str(k)], distance);
	elseif (strcmp(classifier, 'K')) 
		%k-means clustering
		Y_T = kmeans(Y, k, 'Replicates', 1000);
	end
	Y_T = [SRV Y_T];
   	
	%write matrix to file
   	writeCSV(Y, cname, vname2, [inputdir,'/',method], ['Y', '_k', num2str(k),'.csv']);
  	%write survival data and cluster memberships to file
  	if (strcmp(inputdir, 'ucec_mrna') | strcmp(inputdir, 'ucec_mut') | strcmp(inputdir, 'ucec_met')) 
  		writeCSV(Y_T, cname, strvcat('ID', 'Histology', 'Cluster'), [inputdir,'/',method], ['Y_', classifier, '_k', num2str(k),'.csv']);
	else 
  		writeCSV(Y_T, cname, strvcat('ID', 'OS_Event', 'OS_Time', 'Cluster'), [inputdir,'/',method], ['Y_', classifier, '_k', num2str(k),'.csv']);
	end
	plotScatter(Y, k, ctypes, method, ['../data/',inputdir,'/',method,'/']);
   
   elseif (strcmp(method, 'nmf')) %NMF
  	disp('***NMF method***');
	
	%read NMF factorization from file (done in R)
  	inFile2 = ['../data/', inputdir, '/', method, '/', 'nmf_k', num2str(k), '.txt']
	[Z, vnameZ, cnameZ] = tblread(inFile2,' ');
	Z = Z( ismember( cellstr(cnameZ), interXSRV ), : );
  	cnameZ = char(interXSRV);
	
	%DEBUG
	for i=1:length(cname)
		assert(strcmp(cnameZ(i), cname(i)), 'Patient ids in different order!');
	end

	if (strcmp(classifier, 'C') | strcmp(classifier, 'W')) 
		%hierarchical clustering
		Z_T = plotDendrogram(Z, k, ['Z_', classifier, '_k', num2str(k)], distance);
	elseif (strcmp(classifier, 'K')) 
		%k-means clustering
		Z_T = kmeans(Z, k, 'Replicates', 1000);
	end
	Z_T = [SRV Z_T];
	
	%write Z matrix to file
  	writeCSV(Z, cname, vname2, [inputdir,'/',method], ['Z', '_k', num2str(k), '.csv']);
  	%write survival data and cluster memberships to file
  	if (strcmp(inputdir, 'ucec_mrna') | strcmp(inputdir, 'ucec_mut') | strcmp(inputdir, 'ucec_met')) 
  		writeCSV(Z_T, cname, strvcat('ID', 'Histology', 'Cluster'), [inputdir,'/',method], ['Z_', classifier, '_k', num2str(k), '.csv']);
	else 
  		writeCSV(Z_T, cname, strvcat('ID', 'OS_Event', 'OS_Time', 'Cluster'), [inputdir,'/',method], ['Z_', classifier, '_k', num2str(k), '.csv']);
	end
	plotScatter(Z, k, ctypes, method, ['../data/',inputdir,'/',method,'/']);
   
   elseif (strcmp(method, 'xpca')) %XPCA
  	disp('***XPCA method***');
  	%plot the histogram of X values
   	figure('Name','Histogram of X','NumberTitle','off');
   	histogram(X);
   	pause(1);
   	%call xpca
	%initU = initU * initS;
   	%first param is algo, either 'fminunc' or 'newton'
   	[U, V, M] = main(algo, dist, inputdir, X, k, initU, initS, initV, L);
   	% normalize U
	UN = normc(U);
	normColU = sqrt(diag(U' * U));
	D_U = diag(normColU);
	% normalize V
	VN = normc(V);
	normColV = sqrt(diag(V' * V));
	D_V = diag(normColV);
	% new U
	U = UN * D_U * D_V;
	
   	if (strcmp(classifier, 'C') | strcmp(classifier, 'W')) 
		%hierarchical clustering
		U_T = plotDendrogram(U, k, ['U_', classifier, '_k', num2str(k)], distance);
	elseif (strcmp(classifier, 'K')) 
		%k-means clustering
		U_T = kmeans(U, k, 'Replicates', 1000);
	end
	U_T = [SRV U_T];
	
	% Write U matrix to file
  	writeCSV(U, cname, vname2, [inputdir,'/',method], ['U', '_k', num2str(k), '.csv']);
  	% Write survival data and cluster memberships to file
	if (strcmp(inputdir, 'ucec_mrna') | strcmp(inputdir, 'ucec_mut') | strcmp(inputdir, 'ucec_met')) 
  		writeCSV(U_T, cname, strvcat('ID', 'Histology', 'Cluster'), [inputdir,'/',method], ['U_', classifier, '_k', num2str(k), '.csv']);
	else 
		writeCSV(U_T, cname, strvcat('ID', 'OS_Event', 'OS_Time', 'Cluster'), [inputdir,'/',method], ['U_', classifier, '_k', num2str(k), '.csv']);
	end	
	plotScatter(U, k, ctypes, method, ['../data/',inputdir,'/',method,'/']);
   end

end
