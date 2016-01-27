function M = jonesMatrix(taumeanps,Nsect,omega)
    dtau = taumeanps*1e-12/sqrt(Nsect);

    phi = rand(1, 3)*2*pi;

    U1 = [exp(-j*phi(1)/2), 0; 0 exp(j*phi(1)/2)];
    U2 = [cos(phi(2)/2) -j*sin(phi(2)/2); -j*sin(phi(2)/2) cos(phi(2)/2)];
    U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];

    U = U1*U2*U3;
   
    M = repmat(U,[1,1,length(omega)]);

    for k = 1:Nsect
        phi = rand(1, 3)*2*pi;
        U1 = [exp(-j*phi(1)/2), 0; 0 exp(j*phi(1)/2)];
        U2 = [cos(phi(2)/2) -j*sin(phi(2)/2); -j*sin(phi(2)/2) cos(phi(2)/2)];
        U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];

        U = U1*U2*U3;

        for m = 1:length(omega)
            D = [exp(-j*dtau*omega(m)/2), 0; 0, exp(j*dtau*omega(m)/2)];
            M(:,:,m) = M(:,:,m)*U'*D*U;
        end
    end
end