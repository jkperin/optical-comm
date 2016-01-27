function     [Ys, tsamp] = Sampler(t, Y, Nstep, Nsymb,oversamp)
    
    % split into four real signals
    Y1i = real(Y(1,:));
    Y1q = imag(Y(1,:));
    Y2i = real(Y(2,:));
    Y2q = imag(Y(2,:));

    %sample
    firstsamp = 1 + floor(0.5*Nstep);                                     % index of first sample
    sampincr = floor(Nstep/oversamp);
    lastsamp = floor(0.5*Nstep) + oversamp*Nsymb*floor(Nstep/oversamp);   % index of first sample

    sampind = mod(firstsamp:sampincr:lastsamp,Nsymb*Nstep);          % indices of samples in x pol. (modulo wraps around)
    tsamp = t(sampind);           % times at which x pol. sampled

    Y1is = Y1i(sampind);           % samples of inphase in x pol.
    Y1qs = Y1q(sampind);           % samples of quad in x pol.
    Y2is = Y2i(sampind);           % samples of inphase in y pol.
    Y2qs = Y2q(sampind);           % samples of quad in y pol.


    Ys = [Y1is + 1i*Y1qs; Y2is + 1i*Y2qs];

end
