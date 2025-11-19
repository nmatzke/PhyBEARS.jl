# Nearest-neighbor hidden state prediction of discrete or continuous traits for tips on a tree.
hsp_nearest_neighbor = function(tree, 
								tip_states, 	# vector of size Ntips. May include values of any type, including NA (for unknown states)
								check_input	= TRUE){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	
	# basic error checking
	if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if(check_input){
		if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
	}
	
	# find known & unknown tips and calculate nearest distances between them
	known_tips = which(!is.na(tip_states));
	if(length(known_tips)==0) stop("ERROR: All tip states are hidden");
	unknown_tips = which(is.na(tip_states));
	if(length(unknown_tips)==0){
		nearest_distances	= rep(0, times=Ntips)
		nearest_neighbors 	= c(1:Ntips)
		states 				= tip_states
	}else{
		nearest 			= find_nearest_tips(tree, only_descending_tips=FALSE, target_tips=known_tips, check_input=FALSE)
		nearest_distances 	= nearest$nearest_distance_per_tip
		nearest_neighbors 	= nearest$nearest_tip_per_tip
		states				= tip_states
		states[unknown_tips]= tip_states[nearest_neighbors[unknown_tips]]
	}
		
	return(list(success				= TRUE,
				states				= states,
				nearest_neighbors	= nearest_neighbors,
				nearest_distances	= nearest_distances))
}


