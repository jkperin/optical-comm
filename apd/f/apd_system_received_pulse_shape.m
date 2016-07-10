function [HrxPshape, H] = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim)

f = sim.f/sim.fs;
% Channel frequency response. Hch does not include pulse shape
Nhold = sim.Mct/Tx.DAC.ros;
Hpshape = mpam.Hpshape(sim.f);
Hdac = Tx.DAC.filt.H(f).*freqz(1/Nhold*ones(1, Nhold), 1, 2*pi*f).*exp(1j*2*pi*f*(Nhold-1)/2); % ZOH and DAC response
Hmod = Tx.Mod.H(sim.f);
Hfiber = Fiber.Himdd(sim.f, Tx.Laser.wavelength, Tx.Mod.alpha, 'small signal');
Hapd = Apd.H(sim.f);
Hch = Hdac.*Hmod.*Hfiber.*Hapd; % tx pulse shape is not included in channel frequency response

if sim.WhiteningFilter
    Hw = Apd.Hwhitening(sim.f, Tx.Ptx*Fiber.link_attenuation(Tx.Laser.wavelength), Rx.N0);
else
    Hw = ones(size(sim.f));
end

Hch = Hw.*Hch;

Hadc = Rx.ADC.filt.H(f);
switch lower(Rx.filtering)
    case 'antialiasing' % receiver filter is specified in ADC.filt
        Hrx = Hadc;
        HrxPshape = Hpshape.*Hch.*Hrx; % received pulse shape frequency response
    case 'matched' % receiver filter is matched filter
        HrxPshape = Hpshape.*Hch; % received pulse shape frequency response
        Hmatched = conj(HrxPshape); % matched filter matched to received signal pulse shape
        Hrx = Hmatched;
    otherwise
        error('ber_preamp_sys_montecarlo: Rx.filtering must be either antialiasing or matched')
end     

H.pshape = Hpshape;
H.dac = Hdac;
H.mod = Hmod;
H.fiber = Hfiber;
H.apd = Hapd;
H.ch = Hch;
H.w = Hw;
H.adc = Hadc;
H.rx = Hrx;
