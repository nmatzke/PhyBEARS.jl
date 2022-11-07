"""
Update the Qarray and Carray vectors at time t
"""
function update_QC_mats_time_t!(p, t)
 	# The row of bmo that refers to "u", the effect of area on extinction rate
 	# u_row = (1:Rnrow(bmo))[bmo.rownames .== "u"][]
 	max_extinction_rate = p.setup.max_extinction_rate

 	# Get the area of areas at time t
	p.setup.area_of_areas .= p.area_of_areas_interpolator(t)
 	
  # Possibly varying parameters
  n = p.n
  #mu = p.params.mu_vals # base extinction rate of each range
  #mu_t = p.params.mu_t_vals # mu_t = mu at time t
  
  # Populate changing mus with time
  @inbounds for i in 1:n
  	# total_area = get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator)
  	p.params.mu_t_vals[i] = p.params.mu_vals[i] * get_area_of_range(t, p.states_as_areas_lists[i], p.setup.area_of_areas)^p.bmo.est[p.setup.bmo_rows.u_e]
  end
  # Correct "Inf" max_extinction_rates
  p.params.mu_t_vals[p.params.mu_t_vals .> max_extinction_rate] .= max_extinction_rate
  
  
  # Get the e_vals for the Qij matrix, at time t
  update_Qij_e_vals!(p)
  # (updates p.params.Qij_vals)

 # Get the d_vals for the Qij matrix, at time t
  # 1. Update the distance matrices etc.
  p.setup.distmat .= p.distances_interpolator(t)
  
  
  # Update the vicariance minimum distance 
  p.setup.vicdist_t .= p.vicariance_mindists_interpolator(t)
    
  # Using the current t's distmat, etc. update the dmat_t, then 
  # propagate through the Q matrix
  update_Qij_d_vals!(p)
  
  # Using the current t's distmat, etc. update the jmat_t, then
  # propagate through the C matrix
  # Replaces: update_Cijk_j_rates!(p)
	# updates j events
  # updates vicariance also
  update_Cijk_rates!(p)

end # END function update_QC_mats_time_t!(p, t)