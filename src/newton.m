function [x, fval, niter] = newton(fun, x0, citer, crt_axis)

% termination tolerance
tol = 1e-3;

% maximum number of allowed iterations
maxiter = 1000;

% minimum allowed perturbation
dxmin = 1e-6;

% step size ( 0.33 causes instability, 0.2 quite accurate)
%alpha = 1e-9;
%alpha and beta for backtracking line search
alpha = 0.01;
beta = 0.8;

% initialize gradient norm, optimization vector, iteration counter, perturbation
x = x0;
ss = 1 / (size(x, 1) * size(x, 2));
niter = 0; gnorm = inf; dx = inf;

% calculate f, grad
[fval, g, dir] = fun(x);
gnorm = norm(g, 'fro') * ss;

% plot objective function contours for visualization
%ezcontour(f); axis equal; hold on
%iter_vals = [citer + niter];
%f_vals = [f];
plotCrtObj(crt_axis, citer + niter, fval); 

%DEBUG
display(['[newton]niter ' num2str(niter) ',  fval ' num2str(fval) ',  gnorm ' num2str(gnorm)]);

if ~isfinite(dir)
	dir
	error('[newton]dir is infinite or NaN ');
end

%while and(gnorm>tol, and(niter < maxiter, dx > tol))
while and(gnorm>tol, dx>dxmin)
%while (gnorm>tol)
%while (abs(gg) > tol)
    
    t = 1;
    %WRONG WAY OF COMPUTING DIR
    %dir = (-1) * h \ g;
    %dir = (-1) * h_inv * g;
    %%grad descent
    %%dir = (-1) * g;
    %backtracking line search
    xnew = x + t * dir;
    [fval2, g2] = fun(xnew);
    
    %DEBUG
    %%display(['	[newton2]f2 ' num2str(f2) ', o ' num2str(f + alpha * t * g' * dir) ', g2*dir2 ' num2str(abs(g2' * dir2)) ', b*g*dir ' num2str(beta*abs(g' * dir)) ', g2n ' num2str(norm(g2, 'fro')) ', dx ' num2str(norm((xnew - x),'fro')) ', t ' num2str(t)]);
    display(['	[newtonB]fval2 ' num2str(fval2) ',  other ' num2str(fval + sum(sum(alpha * t * g' * dir))) ',  g2norm ' num2str(norm(g2, 'fro'))]);
    
    while (fval2 > (fval + alpha * t * g' * dir)) % & abs(g2' * dir2) > beta * abs(g' * dir))
    	t = beta * t;
    	xnew = x + t * dir;
    	[fval2, g2] = fun(xnew);
    	
	%DEBUG
	%%display(['	[newton2]f2 ' num2str(f2) ', o ' num2str(f + alpha * t * g' * dir) ', g2*dir2 ' num2str(abs(g2' * dir2)) ', b*g*dir ' num2str(beta*abs(g' * dir)) ', g2n ' num2str(norm(g2, 'fro')) ', dx ' num2str(norm((xnew - x),'fro')), ', t ' num2str(t)]);
   	%display(['	[newton2]fval2 ' num2str(fval2) ',  other ' num2str(fval + sum(sum(alpha * t * g' * dir))) ',  g2norm ' num2str(norm(g2, 'fro')) ',  t ' num2str(t)]);
    
   end	
    
    % check step
    if ~isfinite(xnew)
        error(['[newton]niter ' num2str(niter) ', x is inf or NaN ' num2str(xnew)])
    end
    
    %update termination metrics
    %fval = fval2; g = g2; dir2 = dir;
    [fval, g, dir] = fun(xnew);
    gnorm = norm(g, 'fro') * ss;
    dx = norm(xnew - x, 'fro');
    
    niter = niter + 1;
    
    %plot current point
    %plot(crt_axis, citer + niter, f, 'o','LineWidth',5); hold on;
    %iter_vals = [iter_vals; citer + niter];
    %f_vals = [f_vals; f];
    plotCrtObj(crt_axis, citer + niter, fval); 
    
    %DEBUG
    display(['[newtonF]niter ' num2str(niter) ', fval ' num2str(fval) ', gnorm ' num2str(gnorm) ', dx ' num2str(dx) ', t ' num2str(t)]);
    
    x = xnew;
end

%for i=1:length(iter_vals)
%    colRow = rem(i,size(colOrd, 1));
%    if colRow == 0
%    	colRow = size(colOrd, 1);
%    end
%    % Get the color
%    col = colOrd(colRow,:);
%    %plot
%    plot(crt_axis, iter_vals(i), f_vals(i), 'o-', 'LineWidth', 3, 'Color', col); 
%    hold on;
%end
