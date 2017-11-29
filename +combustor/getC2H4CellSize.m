function a = getC2H4CellSize(P)
    size = [0.065 0.061 0.06 0.05 0.04 0.04]';
    pressure = 101325 * [0.3 0.4 0.45 0.5 0.6 0.9]';
    pp = fit(pressure, size, 'exp1');
    a = pp(P);
    P = 101325 * [0.3:0.1:2];
end