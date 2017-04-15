  function [X1, X2, X3, U, V1, M1, V2, M2, V3, M3] = bc_wrapper_multi(method)
  
  %% load methylation data
  [X1, vname1, cname1]= tblread('../data/bc_methyl2_sorted_small_common.txt',' ');

  %% load mRNA data
  [X2, vname2, cname2]= tblread('../data/bc_mRNA_sorted_common.txt',' ');
  
  %% load mutation data
  %[X, vname, cname]= tblread('../data/bc_mutation_sorted_small.txt',' ');
  [X3, vname3, cname3]= tblread('../data/bc_mutation_mr10_sorted_small_common.txt', ' ');
  
    
  %%read cancer types
  fileID = fopen('../data/bc_mRNA_methyl2_mutation_types_types.txt');
  ctypes = textscan(fileID, '%s');
  fclose(fileID);
  
  ctypes = ctypes{1, 1};
  
  colors = ctypes;
  %blue cyan yellow green black
  %cmap = {'LuminalA' 'blue'; 'LuminalB' 'cyan'; 'HER2' 'red'; 'Basal' 'green'; 'Normal' 'black'}
  cmap = {'LuminalA' [0 0 1]; 'LuminalB' [0 1 1]; 'HER2' [1 0 0]; 'Basal' [0 1 0]; 'Normal' [0 0 0]};
  [~, index] = ismember(colors, cmap(:,1));
  colors = cmap(index, 2);
  
  colors = cell2mat(colors);
  
  whos X1
  whos X2
  whos X3
  r1 = rank(X1)
  r2 = rank(X2)
  r3 = rank(X3)
  %whos vname
  %whos cname
  whos ctypes
  %types(1:5)
  %whos color
  %color(1:5)
  %exit('error');
  
  k=5;
  n = size(X1, 1); 
  d1 = size(X1, 2);
  d2 = size(X2, 2);
  d3 = size(X3, 2);

   %HO GSVD init
   X2 = log(1 + X2);
   
   A1 = X1 * X1'; A2 = X2 * X2'; A3 = X3 * X3';
   I1 = inv(A1);
   I2 = inv(A2);
   I3 = inv(A3);
   
   S = 1/6 * (A1 * I2 + A2 * I1 + A1 * I3 + A3 * I1 + A2 * I3 + A3 * I2);
   [U, D] = eig(S);

   U = U(:, 1:5); 
   
   B1 = X1' * U;
   B2 = X2' * U;
   B3 = X3' * U;

   S1 = diag(sqrt(diag(B1' * B1)));
   S2 = diag(sqrt(diag(B2' * B2)));
   S3 = diag(sqrt(diag(B3' * B3)));
   
   V1 = B1 * diag(1 ./ diag(S1));
   V2 = B2 * diag(1 ./ diag(S2));
   V3 = B3 * diag(1 ./ diag(S3));
   
   M1 = zeros(d1, 1); M2 = zeros(d2, 1); M3 = zeros(d3, 1); 
   %plotScatter(U, 'HO GSVD(X)', colors);
   %return;
   
   %plot the histogram of X values
   %figure('Name','Histogram of X','NumberTitle','off');
   %hist(X, 20);
   %pause(1);

   %call exppca
   %first param is method, either 'fminunc' or 'newton'
   %[U, V1, M1, V2, M2, V3, M3] = exppca_multi(method, X1, X2, X3, ctypes, k, U, eye(size(B1, 2)), B1, eye(size(B2, 2)), B2, eye(size(B3, 2)), B3);
   [U, V1, M1, V2, M2, V3, M3] = exppca_multi(method, X1, X2, X3, ctypes, k, U, S1, V1, S2, V2, S3, V3);

   %principal components scatter plot	
   plotScatter(U, 'IPCA(X)', colors);

   %write csv to file
   %writeCSV(X, initU, Y, U, vname, cname, types);
    
