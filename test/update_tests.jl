# Save results
results_fn = "/GitHub/PhyBEARS.jl/test/test_results.txt" # lnLs and times

txtvec = ["test: ", "apes_M0_DEC_v1.jl", "apes", "areas:2", "states:4", "DEC+Yule", "1 like", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, R_bgb_lnL]

append_txtvec(results_fn, txtvec)

write_trdf_ancstates(results_fn, ["botstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")

write_trdf_ancstates(results_fn, ["topstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")



# Save results
txtvec = ["test: ", "apes_M0_DEC_v1.jl", "apes", "areas:2", "states:4", "DEC+Yule", "ML", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, R_bgb_lnL]

append_txtvec(results_fn, txtvec)

write_trdf_ancstates(results_fn, ["botstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")

write_trdf_ancstates(results_fn, ["topstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")




# Save results
txtvec = ["test, v7 ancstates: ", "runtests_ClaSSE_tree_n8_DEC.jl", "Psychotria", "areas:4", "states:16", "DEC+Yule", "1 like", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, DEC_R_result_total_LnLs1]

append_txtvec(results_fn, txtvec)

df1 = df1bot = bgb_ancstates_AT_branchBots_df = R_ancstates_corners
df2 = df2bot = Julia_ancstates_corners_v5
write_trdf_ancstates(results_fn, ["botstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")

df1 = df1top = bgb_ancstates_AT_branchTops_df = R_ancstates_nodes
df2 = df2top = Julia_ancstates_nodes_v5
write_trdf_ancstates(results_fn, ["topstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")

