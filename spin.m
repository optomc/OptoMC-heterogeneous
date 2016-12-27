function dir = spin(g, photonDir)
% It choose a new direction for photon by sampling the polar deflection
% angle theta and azimuthal angle psi
if (g==0)
    costheta=2*rand()-1;
else
    temp=(1-g*g)/(1-g+2*g*rand());
    costheta=(1+g*g-temp*temp)/(2*g);
end
sintheta=sqrt(1-costheta*costheta);

% psi: [0 - 2*pi]
psi=2*pi*rand();
cospsi = cos(psi);
if (psi<pi)
    sinpsi=sqrt(1-cospsi*cospsi);
else
    sinpsi=-sqrt(1-cospsi*cospsi);
end

% New trajectory. 
if (1-abs(photonDir(3)) <= 1E-12)
    tmpDir= [sintheta*cospsi sintheta*sinpsi costheta*sign(photonDir(3))];
else
    tmp=sqrt(1-photonDir(3).^2); 
    tmpDir=[sintheta*(photonDir(1)*photonDir(3)*cospsi-photonDir(2)*sinpsi)/tmp+photonDir(1)*costheta ...
                 sintheta*(photonDir(2)*photonDir(3)*cospsi+photonDir(1)*sinpsi)/tmp+photonDir(2)*costheta ...
                 -sintheta*cospsi*tmp+photonDir(3)*costheta];
end
dir=tmpDir;
end