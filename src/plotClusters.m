function plotClusters(X, IDX, C, k, title)

figure('Name', title ,'NumberTitle','off');

% plot the clusters and the cluster centroids

for i =1:k 

	plot(X(IDX==i, :), X,'r.','MarkerSize',12);
	hold on
	plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12);
	plot(ctrs(:,1),ctrs(:,2),'kx', 'MarkerSize',12,'LineWidth',2);
     	plot(ctrs(:,1),ctrs(:,2),'ko', 'MarkerSize',12,'LineWidth',2);
     	legend('Cluster 1','Cluster 2','Centroids', 'Location','NW');
	hold off;
end
