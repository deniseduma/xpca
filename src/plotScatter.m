function plotScatter(X, k, colors, title, file) 
for i = 1:k-1
	for j = i+1:k
  		figure('Name',title,'NumberTitle','off');
  		crt_axis=gca;
  		scatter(X(:,i),X(:,j),50,colors,'filled');
		legend(crt_axis, [num2str(i),' vs ',num2str(j)]);
  		legend(crt_axis,'show');
		saveas(gcf,[file,num2str(i),'vs',num2str(j)],'epsc');
	end
end	
%pause(2);
end
