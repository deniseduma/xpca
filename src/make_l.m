function L = laplacian(inputdir)
%inputdir is one of bc/ov/ucec/luad/gbm 

fileID = fopen('/home/duma/xpca/data/STRING.txt')
S = textscan(fileID, '%s\t%s');
fclose(fileID);
N = [S{1, 1} S{1, 2}];
B = [S{1, 1}; S{1, 2}];
[n1, n2] = size(N);
[b1, b2] = size(B);
whos N
whos B

%inFile = ['/home/duma/TCGA/data/' inputdir '/' inputdir '_genes.txt']
%fileID = fopen(inFile);
%G = textscan(fileID, '%s');
%fclose(fileID);
%G = G{1, 1};
%g = length(G);
%whos G

matFile = [inputdir, '.txt'];
[X, vname, cname] = tblread(['../data/', inputdir, '/', matFile], ' ');
whos X
whos vname
whos cname

outFile = ['/home/duma/xpca/data/' inputdir '/' inputdir '_laplacian.txt']
outFile2 = ['/home/duma/xpca/data/' inputdir '/' inputdir '_adj.txt']

% Only keep the genes which are in STRING
idx = ismember(cellstr(vname), B);
X = X(: , idx);
vname = vname(idx ,:);
whos X
whos vname

G = cellstr(vname);
g = length(G);

D = zeros(g, g);
A = zeros(g, g);

num_edges = 0;
for index=1:n1 
	[i1, g1] = ismember(N{index, 1}, G);	
	[i2, g2] = ismember(N{index, 2}, G);
	if (i1==1 & i2==1 & g1~=g2)
		A(g1, g2) = 1;
		A(g2, g1) = 1;
		num_edges = num_edges + 1;
	end
end	

D = diag(sum(A, 2));

%% Remove 0 lines from D
idx = all(D==0, 2);

[d1, d2] = size(D)
D(idx, :) = []; D(:, idx) = [];
A(idx, :) = []; A(:, idx) = [];
[d1, d2] = size(D)

[x1, x2] = size(X)
X(:, idx) = [];
[x1, x2] = size(X)
[v1, v2] = size(vname)
vname(idx, :) = [];
[v1, v2] = size(vname)

G = G(~idx);
g = length(G);
whos G

%% Check top degrees
disp('diagonal of D');
[val, index] = sort(diag(D), 'descend');
val(1:min(g, 20))

%% Compute Laplacian
L = D - A;
DD = D^(-1/2);
L = DD * L * DD;

%% Compute degree normalized adjacency matrix
normA = A*D^(-1); 

disp(['num edges 1 ', num2str(num_edges)]);
%% Write network edges to file 
index = 0;
E1 = ''; E2 = '';
%DEBUG
for i=1:g
	for j=1:i-1
		if (A(i, j) == 1)
			index = index + 1;
			E1 = char(E1, G{i, 1});
			E2 = char(E2, G{j, 1});
		end
	end
end	
disp(['num edges 2 ', num2str(index)]);

E1 = E1(2:end, :); E2 = E2(2:end, :);
E = [cellstr(E1) cellstr(E2)]; 
whos E
cell2csv(['/home/duma/xpca/data/' inputdir '/' inputdir '_genes_net_edges.txt'], E, '\t');

%% Write laplacian to file
dlmwrite(outFile, L, 'delimiter', '\t', 'precision', 3);
%% Write adj mat to file
dlmwrite(outFile2, normA, 'delimiter', '\t', 'precision', 3);

%% Write data matrix back to file
%tblwrite(X, vname, cname, matFile,' ');
writeTXT(X, cname, char('ID', vname), inputdir, matFile);

%% Write new gene list back to file
fileID = fopen(['/home/duma/xpca/data/' inputdir '/' inputdir '_genes_net.txt'],'w');
fprintf(fileID,'%s\n', G{:});
fclose(fileID);

end
