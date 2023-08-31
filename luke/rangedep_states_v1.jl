



rownames = ["t12", "t21", "m1", "m2", "m3", "m4"]
type_vec = ["fixed", "fixed", "fixed", "fixed", "fixed", "fixed"]
init_vec = [0.01, 0.01, 0.0, 0.0, 0.0, 0.0]
min_vec = [1.0e-12, 1.0e-12, 0.0, 0.0, 0.0, 0.0]
max_vec = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
est_vec = [0.01, 0.01, 0.0, 0.0, 0.0, 0.0]
note_vec = ["devel", "devel", "devel", "devel", "devel", "devel"]
desc_vec = ["base transition rate from H (hermaphroditic/monoecious) to D (dioecious)", "base transition rate from D (dioecious) to H (hermaphroditic/monoecious)", "multiplier on t12 when in range L", "multiplier on t21 when in range L", "multiplier on t12 when in range M", "multiplier on t21 when in range M"]

bmo_add = DataFrames.DataFrame(rownames=rownames, type=type_vec, init=init_vec, min=min_vec, max=max_vec, est=est_vec, note=note_vec, desc=desc_vec)

# Specify the Qij matrix
i = [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6];
j = [3, 4, 3, 5, 1, 2, 6, 1, 6, 2, 6, 3, 4, 5];
M0 = [d, t12, d, t12, e, e, t12, t12, d, t12, d, t12, e, e];
M1 = [d, t12, d, t12, e, e, t12, t21, d, t21, d, t21, e, e];
M2 = [d, t12, d, m3t12, e, e, m3t12, t12, d, m3t12, d, m3t12, e, e];
M3 = [d, t12, d, m3t12, e, e, m3t12, t21, d, t21, d, t21, e, e];
M4 = [d, t12, d, t12, e, e, t12, t21, d, m4t21, d, m4t21, e, e];
M5 = [d, t12, d, m3t12, e, e, m3t12, t21, d, m4t21, d, m4t21, e, e];


