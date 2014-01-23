function [D, K, t] = Delta (tmax, tau, theta, k1, k2)
	M = (tmax - theta) / tau;
	N = floor(tmax);

	T1 = 0:N;
	for j = 0:M
		T2 (j+1) = j*tau + theta;
	end

	T (1,1:N+1) = T1;
	T (1,N+2:N+M+2) = T2;

	% we sort all those time indices
	[T, i] = sort (T);

	% we have to filter out all those
	% time entries that repeat themselves
	k = 1;
	for j = 1:N+M+1
		if (T(j+1) != T(j))
			t(k) = T(j);
			k++;
		endif
	end

	N = length (t);

	% now we create vectors D and K
	for j = 1:N-1
		D(j) = t(j+1) - t(j);
		if (mod(t(j), 1) == 0)
			if (mod(t(j), tau) - theta == 0)
				K(j) = k1 + k2;
			else
				K(j) = k1;
			endif
		else
			K(j) = k2;
		endif
	end

	return;
endfunction
