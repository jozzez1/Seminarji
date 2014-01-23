function [x, p] = phaseD (x1, p1, D, K)
	x(1) = x1;
	p(1) = p1;

	N = length (D);
	for i = 1:N-1
		p(i+1) = p(i) + D(i)*K(i)*sin(x(i));
		x(i+1) = mod(x(i) + D(i)*p(i+1), 2*pi);
	end

	return;
endfunction
