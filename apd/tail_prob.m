function Ptail = tail_prob(x, b, apdd, ther_var, dt)

level_pdf = apdd.levels_pdf(b + x, dt);

I = level_pdf.I - level_pdf.mean;

dI = abs(I(1) - I(2));

xtilde = x*apdd.Gain;
if ther_var ~= 0
    thermal_pdf = pdf('normal', I, 0, sqrt(ther_var));

    noise_pdf = dI*conv(level_pdf.p, thermal_pdf, 'same');
else
    noise_pdf = level_pdf.p;
end

% plot(I, noise_pdf, I(I <= -xtilde), noise_pdf(I <= -xtilde))

if sum((I <= -xtilde)) > 1
    Ptail = trapz(I(I <= -xtilde), noise_pdf(I <= -xtilde));
elseif sum((I <= -xtilde)) == 1
    Ptail = dI*noise_pdf(I <= -xtilde);
else
    Ptail = 0;
end

1;

