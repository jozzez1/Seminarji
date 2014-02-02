function [X, P, t] = phases (N, tmax, tau, theta, k1, k2)
	% first of all, we calculate D, K and t -- we need those
	[D, K, t] = Delta (tmax, tau, theta, k1, k2);

	% now calculate X and with different initial positions
	for j = 1:N
		x1 = rand(1)*2*pi;
		p1 = rand(1)*28 - 14;
		[x, p] = phaseD (x1, p1, D, K);

		% phaseD returnes line-vectors, but we want 
		% row-vectors
		X(:,j) = x.';
		P(:,j) = p.';
	end

	% for the sake of consistency, let's also make 't'
	% a row-vector
	t = t.';

	% let's also arrange the results for plotting
	[nr, nc] = size(X);
	n = nr * nc;
	A(:,1) = reshape (X, n, 1);
	A(:,2) = reshape (P, n, 1);

	save -ascii prelim.txt A;

	return;
endfunction
