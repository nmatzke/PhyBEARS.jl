

# Specify the Qij matrix
i = [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6];
j = [3, 4, 3, 5, 1, 2, 6, 1, 6, 2, 6, 3, 4, 5];
M0 = [d, t12, d, t12, e, e, t12, t12, d, t12, d, t12, e, e];
M1 = [d, t12, d, t12, e, e, t12, t21, d, t21, d, t21, e, e];
M2 = [d, t12, d, m3t12, e, e, m3t12, t12, d, m3t12, d, m3t12, e, e];
M3 = [d, t12, d, m3t12, e, e, m3t12, t21, d, t21, d, t21, e, e];
M4 = [d, t12, d, t12, e, e, t12, t21, d, m4t21, d, m4t21, e, e];
M5 = [d, t12, d, m3t12, e, e, m3t12, t21, d, m4t21, d, m4t21, e, e];


