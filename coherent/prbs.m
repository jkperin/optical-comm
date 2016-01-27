function b = prbs(n,taps,N)
% n    : Number of bits in the shifte-register
% taps : Bit positions that affect the next stat of SR.
% N    : Number of bits.
    if (nargin < 1)
        n=7;
        taps = [6,7];
    end
    
    if (nargin < 2)
        switch (n)
            case 7
                taps = [6 7];
            case 8
                taps = [4 5 6 8];
            case 9
                taps = [5 9];
            case 10
                taps = [7 10];
            case 11
                taps = [9 11];
            case 12
                taps = [4 10 11 12];
            case 13
                taps = [8 11 12 13];
            case 14
                taps = [2 12 13 14];
            case 15
                taps = [14 15];
            otherwise   
                error('Not enough input arguments');
        end
    end
       
    if (nargin < 3)
        N = 2^n-1;
    end
    
    while(1)
        RandomSeed = randi(2,1,n)-1;
        if sum(RandomSeed)>0, break; end;
    end

    b = zeros(1,N);
    SR = RandomSeed;                % Shift Register
    Taps = zeros(1,n);
    Taps(taps) = ones(size(taps));
    
    for i = 1:N
        b(i) = SR(n);
        SR_1 = mod(SR*Taps',2);
        SR(2:n)=SR(1:n-1);
        SR(1) = SR_1;   
    end
        
