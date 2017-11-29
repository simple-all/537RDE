x = [10 20 30 40 50];
nmi = [280 427 805 1009 1243];
fmf = [0.19 0.28 0.32 0.34 0.38];

pp_nmi = polyfit(x, nmi, 2);
pp_fmf = polyfit(x, fmf, 2);

xs = [1:1:50];
for i = 1:numel(xs)
    nms(i) = polyval(pp_nmi, xs(i));
    fms(i) = polyval(pp_fmf, xs(i));
end

% Vehicle params
lbm_s = 3.35 * 2.20462; % lbm/s
range = 183.5; % 0.1 fmf
range = 354; % 0.17 fmf
range = 305; % 0.15 fmf

fmf_v = polyval(pp_fmf, lbm_s);

figure;
hold on;
plot(x, fmf, 'x');
plot(xs, fms);
xlabel('Air Mass Flow (lbm/s)');
ylabel('Fuel Mass Fraction');


figure;
hold on;
plot(x, nmi, 'x');
plot(xs, nms);
xlabel('Air Mass Flow (lbm/s)');
ylabel('Range (Nautical Miles)');
% Operating point
plot(lbm_s, range, 'kx');