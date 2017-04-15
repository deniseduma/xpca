function process_mrna_mat(wd)

%[X, cname, vname] = tblread('../data/bc_mRNA_original_tum.txt','\t');
[X, cname, vname] = tblread(['../data/', wd, '/', wd, '_mRNA_original_firehose.txt'],'\t');

X = X';
[n, d] = size(X)

whos X
whos cname
whos vname

cnamearr = cellstr(cname);
vname = strvcat('PID',  vname);
vnamearr = cellstr(vname);

whos cnamearr
whos vnamearr

%shorten patient ids and only keep primary tumor patients
%for i = 1:n
%	tokens = strsplit(cname(i, :), '-');
%	cname(i, :) = strjoin(tokens(1:3), '-');
%end
toKeep = [];
for i = 1:n
	remain = cnamearr{i};
	if (strfind(remain, '-01A-')) 
		toKeep = [toKeep i];
	end
	[str1, remain] = strtok(remain, '-');
   	[str2, remain] = strtok(remain, '-');
   	[str3, remain] = strtok(remain, '-');
	cnamearr(i) = cellstr([str1 '-' str2 '-' str3]);
end

X = X(toKeep, :);
cnamearr = cnamearr(toKeep);

[n1, d1] = size(X)

%retain only top 2000 most variable genes
stdX = std(X);
[val, index] = sort(stdX, 'descend');
[i1, i2] = size(index)
index = index(1:2000);
[i1, i2] = size(index)
X = X(:, index);
vnamearr = vnamearr([1 index + 1]);

[n2, d2] = size(X)

whos vnamearr
vnamearr(1:5)

%remove | from gene name 
vnamearr = regexprep(vnamearr, '\|(\d+)', '');

whos vnamearr
vnamearr(1:5)

xcellarr=mat2cell(X, ones(n2,1), ones(d2,1));
xcellarr=[cnamearr xcellarr];
xcellarr = [vnamearr'; xcellarr];

%write X cell array to file
%cell2csv('../data/bc_mRNA_tum.txt', xcellarr, ' ');
cell2csv(['../data/', wd, '/', wd, '_mRNA_tum_firehose.txt'], xcellarr, ' ');

end
