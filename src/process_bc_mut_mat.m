%[X, vname, cname] = tblread('../data/Mutations.CNAabs1.cat.noY.m05.txt','\t');
[X, vname, cname] = tblread('../data/CNA.part.matrix.thabs1.cat.noY.txt','\t');

[n, d] = size(X)

whos X
whos cname
whos vname

% check the presence of somatic mutations
numMut = length(find(X==2))

%replace -1 and 2 by 1 
X(X==-1) = 1;
X(X==2) = 1;

%only keep genes that appear in at least 10% of the patients
sumX = sum(X);
indexX = find(sumX >= 0.1 * n);
X = X(:, indexX);
vname = vname(indexX, :);

%now remove patients which are left with no mutation
sumX = sum(X');
indexX = find(sumX ~= 0);
X = X(indexX, :);
cname = cname(indexX, :);

[n1, d1] = size(X)

whos X
whos vname
whos vname

vname = strvcat('PID',  vname);

cnamearr = cellstr(cname);
vnamearr = cellstr(vname)';

xcellarr=mat2cell(X,ones(n1,1), ones(d1,1));
xcellarr=[cnamearr xcellarr];
xcellarr = [vnamearr; xcellarr];
whos xcellarr

%write X cell array to file
cell2csv('../data/bc_mut_cna.txt', xcellarr, ' ');
