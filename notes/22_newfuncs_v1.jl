

# Translate SSE results to tables
states_list_txt = ["A", "B", "AB"]

Rorder = sort(trdf, "Rnodenums").nodeIndex;
trdf[Rorder,:]
downpass_normlikes_by_Rnode = vec_of_vecs_to_df(res.normlikes_at_each_nodeIndex_branchTop[Rorder]; colnames=states_list_txt)

downpass_branchBots_by_Rnode = vec_of_vecs_to_df(res.normlikes_at_each_nodeIndex_branchBot[Rorder]; colnames=states_list_txt)

downpass_normlikes_by_Rnode
downpass_branchBots_by_Rnode

trdf[Rorder,:][149,:]
downpass_normlikes_by_Rnode[149,:]
downpass_branchBots_by_Rnode[149,:]

outfn = "downpass_normlikes_by_Rnode.txt"
CSV.write(outfns, downpass_normlikes_by_Rnode; delim="\t")
outfn = "downpass_branchBots_by_Rnode.txt"
CSV.write(outfns, downpass_branchBots_by_Rnode; delim="\t")


"""
vec = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
states_list_txt = ["A", "B", "AB"]
vec_of_vecs_to_df(vec; colnames=:auto)
vec_of_vecs_to_df(vec; colnames=states_list_txt)
"""
function vec_of_vecs_to_df(vec; colnames=:auto)
	ncols = length(vec[1])
	df = DataFrame([getindex.(vec, i) for i in 1:ncols], colnames, copycols=false)
	return(df)
end # END function vec_of_vecs_to_df(vec; colnames=:auto)


