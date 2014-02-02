function [b,lam,cas] = ljapunov1 (tmax, T, tau, theta, k1, k2)
	% first let's get some operational parameters
	d0 = 1e-6;
	d = rand(1,2);
	d = d0 * (d ./ norm(d));

	x(1) = rand(1)*2*pi;
	x(2) = rand(1)*18 - 14;

	% a small perturbation
	y(1) = x(1) + d0;
	y(2) = x(2);

	[D, K, t] = Delta (tmax, tau, theta, k1, k2);
	N = length(t);

	% every time when mod(t, 10) == 0 we will freeze the lyapunov stuff
	L = 0;
	k = 1;

	% now we have to propagate those two bastards
	for i = 1:N-1
		x(2) += K(i)*sin(x(1));
		x(1) += D(i)*x(2);

		y(2) += K(i)*sin(y(1));
		y(1) += D(i)*y(2);

		if mod(t(i),T) == 0
			dt   = norm(x(1) - y(1));
			if dt > 0
				L   += log(dt/d0);
	
				lam(k) = L;
				cas(k) = t(i);

				k++;
	
				% and now we have to rescale the 'y' vector
				y(1) = x(1) + d0;
			endif
		endif
	end

%	L = L/(cas(length(cas)));
	% we will L through linear fit
	dol = length(cas);
	kappa = 1;
	for k = 2:dol
		if (cas(k) >= 2000) && (cas(k-1) < 2000)
			kappa = k;
		endif
	end

	Y      = lam.'(kappa:dol);
	X(:,1) = cas.'(kappa:dol);
	X(:,2) = ones(length(cas)-kappa+1, 1); 

	[b, sigma, r] = ols (Y, X);


	return;
endfunction
