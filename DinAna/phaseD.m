function [x, p] = phaseD (x1, p1, D, K)
	x(1) = x1;
	p(1) = p1;

	N = length (D);
	for i = 1:N-1
		p(i+1) = p(i) + K(i)*sin(x(i));
		x(i+1) = x(i) + D(i) * p(i+1);
	end

	return;
endfunction
