function ptail = calc_ptail_fix_thresh(Ptx, Pth, soa, rx, OptFilt, EleFilt, sim)

[px, x] = px_soa_klse_freq(Ptx, soa, rx, OptFilt, EleFilt, sim);

ptail = trapz(x(x < Pth), px(x < Pth));
