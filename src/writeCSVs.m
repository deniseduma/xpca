function writeCSV(X, initU, Y, U, vname, cname, type)

  [n d] = size(X)
  [ny dy] = size(Y)
  k = size(U, 2);
  
  comps = [];
  for i=1:k
  	comps = [comps; int2str(i)];
  end

  vname = strvcat('ID',  vname, 'TYPE');
  vname2 = strvcat('ID', comps, 'TYPE');
  
  %write table to file (requires newer Matlab version)
  %T = cell2table(cacx, 'VariableNames', mat2cell(vname), 'RowNames', mat2cell(cname));
  cnamearr = cellstr(cname);
  vnamearr = cellstr(vname)';
  vnamearr2 = cellstr(vname2)';
  
  xcellarr=[mat2cell(X,ones(n,1), ones(d,1)) type];
  xcellarr=[cnamearr xcellarr];
  xcellarr = [vnamearr; xcellarr];
  whos xcellarr
  %write X cell array to file
  cell2csv('X.csv', xcellarr, ',');
  
  ycellarr = [mat2cell(Y,ones(ny,1), ones(dy,1)) type];
  ycellarr = [cnamearr ycellarr];
  if (dy == d)
  	ycellarr = [vnamearr; ycellarr];
  else	
	ycellarr = [[vnamearr(1) vnamearr(2:dy + 1) vnamearr(end)]; ycellarr];
  end
  whos ycellarr
  %write Y2 cell array to file
  cell2csv('Y.csv', ycellarr, ',');
  
  uinitcellarr = [mat2cell(initU,ones(n,1), ones(k,1)) type];
  uinitcellarr = [cnamearr, uinitcellarr];
  uinitcellarr = [vnamearr2; uinitcellarr];
  whos uinitcellarr
  %write U2 cell array to file
  cell2csv('initU.csv', uinitcellarr, ',');
  
  ucellarr = [mat2cell(U,ones(n,1), ones(k,1)) type];
  ucellarr = [cnamearr, ucellarr];
  ucellarr = [vnamearr2; ucellarr];
  whos ucellarr
  %write U cell array to file
  cell2csv('U.csv', ucellarr, ',');

end
