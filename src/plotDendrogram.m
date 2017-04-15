function T = plotDendrogram(X, p, title, dist)

[n, d] = size(X)
Z = linkage(X, dist);
[nz, dz] = size(Z)
%colors = get(gca,'ColorOrder');
%figure('Name', title ,'NumberTitle','off');
[H, T, P] = dendrogram(Z, p);
%imagesc((colors(Perm, :) * [4 2 1]')')
display('size(T)');
size(T)
display('size(P)');
size(P)

end

 
