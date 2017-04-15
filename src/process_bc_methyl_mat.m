[X, cname, vname]= tblread('../data/bc_methyl.txt',' ');

X = X';
[n, d] = size(X)

whos X
whos cname
whos vname

%%%replace any 0 in X by the min in that col%%%
%replace each 0 in X by 100
X(X==0) = 100;
%find min in X col by col
minX = min(X);
%map columns of X to their min
mymap = [(1:d)' minX'];
assert(size(mymap, 1) == d, 'size(mymap, 1) != d');
assert(size(mymap, 2) == 2, 'size(mymap, 2) != 2');
%find indices of 100 in X
[r100 c100] = find(X==100);
%map c100 to corresponding min per col
[~, index] = ismember(c100, mymap(:, 1));
%replace 100 by corresponding min per col
X(X==100) = mymap(index, 2);

%apply the probit function to the data
%X = norminv(X , 0, 1);
%apply the logit function to the data
%X = log(X ./ (1 - X));

vname = strvcat('PID',  vname);

cnamearr = cellstr(cname);
vnamearr = cellstr(vname)';

xcellarr=mat2cell(X,ones(n,1), ones(d,1));
xcellarr=[cnamearr xcellarr];
xcellarr = [vnamearr; xcellarr];
whos xcellarr

%write X cell array to file
%cell2csv('../data/bc_methyl2.txt', xcellarr, '\t');
%write data without probit/logit transformation
cell2csv('../data/bc_methylO.txt', xcellarr, '\t');
