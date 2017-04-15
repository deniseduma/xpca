function writeCSV(X, cname, vname, inputdir, name)
  
  [n, d] = size(X);

  cnamearr = cellstr(cname);
  vnamearr = cellstr(vname);
  
  xcellarr=mat2cell(X, ones(n,1), ones(d,1));
  xcellarr=[cnamearr xcellarr];
  xcellarr = [vnamearr'; xcellarr];
  
  %write X cell array to file
  cell2csv(['../data/', inputdir, '/', name], xcellarr, ',');
