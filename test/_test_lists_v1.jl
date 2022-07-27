cd("/GitHub/PhyBEARS.jl/test/")
include("runtests_trad_vs_Gflow_v3.jl")
include("runtests_trad_vs_Gflow_v3_linreg.jl")
include("runtests_trad_vs_Gflow_v4.jl")
include("runtests_trad_vs_Gflow_v4_linreg.jl")


include("runtests_trad_vs_Gflow_v5_Psychotria.jl")
include("runtests_trad_vs_Gflow_v5_Psychotria_linreg.jl")

# Show that correlation is high between trad SSE and Gflow (and somewhat to Gflow_arrays, 
# although more issues there.)

# Check likelihood calculations with extinction and +J -- works!
include("runtests_trad_vs_Gflow_v5_Psychotria_wExtinction+J.jl")

# Check likelihood calculations across a variety of parameter values
include("/GitHub/PhyBEARS.jl/test/runtests_trad_vs_Gflow_v5_Psychotria_linreg_NoCrap.jl")

# Speed tests
include("/GitHub/PhyBEARS.jl/test/speedtests_Cyrtandra_wExtinction+J.jl")

include("")
include("")
include("")
include("")
