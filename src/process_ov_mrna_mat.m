[X, cname, vname] = tblread('../data/ov_mRNA_original_tum.txt','\t');

X = X';
[n, d] = size(X)

whos X
whos cname
whos vname

%trim patients ids
cnamearr = cellstr(cname);

for i = 1:n
	remain = cnamearr{i};
   	[str1, remain] = strtok(remain, '-');
   	[str2, remain] = strtok(remain, '-');
   	[str3, remain] = strtok(remain, '-');
	cnamearr(i) = cellstr([str1 '-' str2 '-' str3]);
end

%retain only top 1000 most variable genes
stdX = std(X);
[val, index] = sort(stdX, 'descend');
size_index1 = size(index)
index = index(1:1000);
size_index2 = size(index)
X = X(:, index);
vname = vname(index,  :);

[n1, d1] = size(X)

%write X cell array to file
vname = strvcat('PID',  vname);
vnamearr = cellstr(vname)';

xcellarr=mat2cell(X, ones(n1,1), ones(d1,1));
xcellarr=[cnamearr xcellarr];
xcellarr = [vnamearr; xcellarr];

whos xcellarr

cell2csv('../data/ov_mRNA_tum.txt', xcellarr, ' ');
