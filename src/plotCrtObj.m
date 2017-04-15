function plotCrtObj(crt_axis, iter, fval)

colOrd = get(crt_axis,'ColorOrder');

m = size(colOrd, 1);
colRow = rem(iter,m);
if colRow == 0
	colRow = m;
end
% Get the color
col = colOrd(colRow,:);
plot(crt_axis, iter, fval, 'o-', 'LineWidth', 3, 'Color', col); 
pause(0.03); hold on;

end
