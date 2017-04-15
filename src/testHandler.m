function testHandler

a = 1; b = 1; c = 1;

%handler = @(a)funHandler(a, b, c);

for i=1:10
	a = a + 1;
	b = b * 2;
	c = c * 3;
	handler = @(a)funHandler(a, b, c);
	fprintf('i=%d\n', i);
	handler(a);
end

end
