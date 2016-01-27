% vprop.m   Milad Sharif  02-19-13
function Efilt = vprop(Ein,M)
% Polarization-mode dispersion
	Epmd1 = ifft(fft(Ein(1,:)).*reshape(M(1,1,:),1,size(M,3)) + fft(Ein(2,:)).*reshape(M(1,2,:),1,size(M,3)));
    Epmd2 = ifft(fft(Ein(1,:)).*reshape(M(2,1,:),1,size(M,3)) + fft(Ein(2,:)).*reshape(M(2,2,:),1,size(M,3)));
    Efilt = [Epmd1; Epmd2];    
end