function writeTXT(X, cname, vname, inputdir, name)
  
  [n d] = size(X);
  
  %whos X
  %whos cname
  %whos vname

  cnamearr = cellstr(cname);
  vnamearr = cellstr(vname);
  
  xcellarr=mat2cell(X, ones(n,1), ones(d,1));
  xcellarr=[cnamearr xcellarr];
  xcellarr = [vnamearr'; xcellarr];
  %whos xcellarr
  
  %write X cell array to file
  cell2csv(['../data/', inputdir, '/', name], xcellarr, ' ');

