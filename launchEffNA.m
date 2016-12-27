function dir = launchEffNA(NA, nt)
    % The direction cosines of the initial photons are modelled to have even
    % spread angle between 0 and NA
    theta = asin(NA/nt);  
    cost = rand()*(cos(theta)-1)+1;
    sint = sqrt(1-cost^2);
	psi = 2.0 * pi * rand();
	cosp = cos(psi);
    sinp = sin(psi);
    dir = [sint*cosp sint*sinp cost];    
end