function plotSurvival(k)

[X, vname, cname]= tblread('../data/bc_mRNA_clusters_U.csv',',');
[n, d] = size(X)
%X = X(find(X(:, 2)~=0), :);
%[n1, d1] = size(X)


%colors for the plots
cc = ['r' 'g' 'b' 'm' 'k' 'y' 'c'];
mycc = cc(1:k)


% plot cdf for each cluster
mylegend = [];
figure('Name','ecdf','NumberTitle','off');

for index=1:k
	%%% M1
	%X3 = X(:, 3);
	%rows = find(X3==i);
	%Y = X(rows, :)
	%%% M2
	%X3 = X(:, 3);
	%Y = X(X3==index, :)
	% M3
	Y = [X(X(:,3)==index, 2), X(X(:,3)==index, 1)];
	[ny, dy] = size(Y)
	if (ny >= 5)
		% estimate empirical survivor function
		[ff, xx] = ecdf(Y(:,1), 'censoring', ~Y(:,2), 'function','survivor');
		%plot estimated survivor functions
		stairs(xx, ff, mycc(index), 'LineWidth', 2);
		hold on
		mylegend = strvcat(mylegend, ['Cluster ' num2str(index) ' ']);
		
		% fit Weibull distributions to OS times 
		%pd = fitdist(Y(:,2),'wbl','censoring', Y(:,1))
		%whos pd;
		%plot(0:1:max(x), 1 - cdf('wbl', 0:1:max(x), pd.a, pd.b), ['--' mycc(i)])
		%hold on
		%mylegend = strvcat(mylegend , ['WBL ' num2str(index) ' ']);
	end
end
hold off

% add legend
legend(mylegend, 'Location', 'NorthEast', 'Orientation', 'Vertical')
xlabel('Days to death/last contact');
ylabel('Survival probability');

% fit Cox prop hazards regression where clusters are the explanatory variable
[b, logL, H, stats] = coxphfit(X(:,3), X(:,2), 'censoring', ~X(:,1));
stats

%figure('Name','coxphfit','NumberTitle','off');

%sx = size(X)
%sh = size(H)
%stairs(H(:, 1), exp(-H(:, 2)), 'LineWidth', 2)

% plot the Cox estimate of the baseline survivor function
%for i=1:k
	%if (ny >= 5)
%		index = find(X(:, 3) == i);
%		si = size(index)
%		stairs(H(index', 1), exp(-H(index', 2)), 'LineWidth', 2)
	%end	
%end

% add legend
%legend(mylegend, 'Location', 'NorthEast', 'Orientation', 'Vertical')
%xlabel('Days to death/last contact');
%ylabel('Survival probability');

end

