#######################################################
# (Archived version of Wallis's code; she has further updates in a diff repository)
#######################################################

module pygplatesfunctions
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process  skip the precompile  caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

print("\n\nStarting module 'pygplatesfunctions'...loading dependencies...\n")
using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA           # for lsoda()
using Sundials        # for CVODE_BDF(linear_solver=:GMRES)
using DifferentialEquations
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks
#using Plots						# for plot
using DataFrames          # for DataFrame()
using Query               # for Querying dataframes!
using PhyBEARS.TrUtils # for flat2() (similar to unlist)
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs

# (1) List all function names here:
export sayhello4, txtdf_read, csvdf_read, location_used, land_last_touch, distance_given, distance_interp, oldest_polygon


"""

can really just use:

using Dataframe
df = readtable("Desktop/School/Rhacophoridae/Gplates/output.csv")
df = readtable("Desktop/School/Rhacophoridae/Gplates/output.txt", separator = \t)

"""


#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/PhyBEARS.jl/notes/")
include("pygplatesfunctions.jl")
"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

say_hello4() = println("Hello dude4!")
say_hello4()


function test_error(x)
    #blahblahblah
    y = x - 1

    if x == 0
        txt = "hmm"
        error("testing error ", x, " didn't work")
    end
    return y
end

# df = readtable("Desktop/School/Rhacophoridae/Gplates/output.csv")
# time = 55
# land1 = "India"
# land2 = "Sundaland"


function location_used(df, time, land1, land2)
    #blah blah blah

    
end


"""
To test queries:

df = DataFrame(name=["John", "Sally", "Roger", "Alice"], age=[54., 34., 79., 41.], children=[0, 2, 4, 5], adult=["Alice", "Bob", "Charlie", "John"])

q = @from i in df begin
    @where i.age > 40 && i.children > 0
    @select i.name
    @collect
end


q2 = @from i in df begin
    @where i.age > 40 && ((i.name == "John" && i.adult == "Alice") || (i.name == "Alice" && i.adult=="John"))
    @select i.children
    @collect
end


"""

"""
    distance_given(df, time, land1, land2)

Provides the distance at a given time point, provided that the time point requested is
along the interval originally outputted from pygplates into the dataframe. Allows the 
user to quickly access distance already recorded within their dataframe. 

For interpolated distances at time points not previously recorded, please use 'distance_interp'

* `df` - dataframe created by pygplates output created by Wallis Bland
contains variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

* `time` - timepoint for distance requested. Within distance_given, this timepoint should be
within the original dataframe.

* `land1` - first landmass chosen for distance comparison. Ensure this is within "string" format.

* `land2` - second landmass chosed for distance comparison. Ensure this is within "string" format

NOTE: Lands 1 and 2 do not need to line up with Land1 or 2 within the dataframe in this function, as it will search both!

# Examples
```julia-repl
julia> test_df = DataFrame(Reconstruction_Time_Ma=[0, 10, 20, 30, 40, 50, 60, 0, 10, 20, 30, 40, 50, 60], Land1=["India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India"], Land2=["Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar","Madagascar"], Closest_Distance_km=[0, 1000, 2000, 3500, 4100, 4900, "NA", 0, 1000, 2000, 3500, 4100, 4900, "NA"], Point1_Lat=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point1_Lon=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point2_Lat=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"], Point2_Lon=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"])

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1   Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Int64                   String  String      Any                  Any         Any         Any         Any        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │                      0  India   Asia        0                    1           1           1           1
   2 │                     10  India   Asia        1000                 2           2           2           2
   3 │                     20  India   Asia        2000                 3           3           3           3
   4 │                     30  India   Asia        3500                 4           4           4           4
   5 │                     40  India   Asia        4100                 5           5           5           5
   6 │                     50  India   Asia        4900                 6           6           6           6
   7 │                     60  India   Asia        NA                   NA          NA          7           7
   8 │                      0  India   Madagascar  0                    1           1           1           1
   9 │                     10  India   Madagascar  1000                 2           2           2           2
  10 │                     20  India   Madagascar  2000                 3           3           3           3
  11 │                     30  India   Madagascar  3500                 4           4           4           4
  12 │                     40  India   Madagascar  4100                 5           5           5           5
  13 │                     50  India   Madagascar  4900                 6           6           6           6
  14 │                     60  India   Madagascar  NA                   NA          NA          NA          NA


julia> distance_given(test_df, 40, "India", "Asia")
DataValue{Any}(4100)

julia> distance_given(test_df, 40, "Asia", "India")
DataValue{Any}(4100)


# POTENTIAL ERRORS:

* TTIMESTAMP NOT WITHIN TABLE
    Variables potentially identified, but time stamp not within the written intervals of the dataframe. Use function pygplatesfunctions.distance_interp() instead.

    julia> distance_given(test_df, 44, "Asia", "India")
    ERROR: STOP ERROR in pygplatesfunctions.distance_given().
     Pairing: 'Asia' and 'India' has not been found within given dataframe at this timestamp.
     Please check spelling of compared landmasses.
     If using function distance_given(), please check dataframe timestamps.
 	    if time intervals taken from original pygplates do not have your needed timestamp, please us function pygplatesfunctions.distance_interp() 
     If using function distance_interp(), please ensure time intervals within dataframe go PAST the requested timestamp. Thank you! 



* INCORRECT VARIABLES / MISSPELLINGS
    In this case the variable was either spelt wrong or for some other reason Variable was never printed within the given dataframe.

    julia> distance_given(test_df, 40, "Asia", "Indai")
    ERROR: STOP ERROR in pygplatesfunctions.distance_given().
     Pairing: 'Asia' and 'Indai' has not been found within given dataframe at this timestamp.
     Please check spelling of compared landmasses.
     If using function distance_given(), please check dataframe timestamps.
        if time intervals taken from original pygplates do not have your needed timestamp, please us function pygplatesfunctions.distance_interp() 
     If using function distance_interp(), please ensure time intervals within dataframe go PAST the requested timestamp. Thank you! 



* INCORRECT PAIRING
    Lands exist within the dataframe, but were never compared distance wise 
    If original pygplates code was run correctly, this SHOULD NOT OCCUR. More likely need to check spelling

    julia> distance_given(test_df, 40, "Asia", "Madagascar")
    ERROR: STOP ERROR in pygplatesfunctions.distance_given().
     Pairing: 'Asia' and 'Madagascar' has not been found within given dataframe at this timestamp.
     Please check spelling of compared landmasses.
     If using function distance_given(), please check dataframe timestamps.
        if time intervals taken from original pygplates do not have your needed timestamp, please us function pygplatesfunctions.distance_interp() 
     If using function distance_interp(), please ensure time intervals within dataframe go PAST the requested timestamp. Thank you! 



* BOTH LANDS MISSING POLYGONS
    Variables identified, but polygons for original gplates gpml file for BOTH landmasses were nvever identified at that timestamp

    julia> distance_given(test_df, 60, "India", "Madagascar")
    ERROR: STOP ERROR in pygplatesfunctions.distance_given(). Both landmasses, 'India' and 'Madagascar', 
     do not have polygons at this timestamp (time=60).
     Please use land_begin to find the earliest timestamp with polygons present. If land_begin returns an earlier timestamp, please check gplates file.



* SINGLE LAND MISSING POLYGONS
    Variables identified, but polygons from original gplates gpml file were never identified at that timestamp.
    Note below, whether India is land 1 or 2, it is noted that that is the one missing polygons


    julia> distance_given(test_df, 60, "India", "Asia") 
    ERROR: STOP ERROR in pygplatesfunctions.distance_given(). The landmass, 'India', 
     does not have a polygon at this timestamp (time=60).
     Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.

    
    julia> distance_given(test_df, 60, "Asia", "India")
    ERROR: STOP ERROR in pygplatesfunctions.distance_given(). The landmass, 'India', 
     does not have a polygon at this timestamp (time=60).
     Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.
```
"""

# Hey Nick! A question here! mshould land1  land2 be input as strings when someone uses the fuction?
# i.e: distance_given(df, 40, "India", "Asia") ??


# or should I enter them as "land1", "land2" within this definition of the function?



function distance_given(df, time, land1, land2)
    
    q = @from i in df begin
        @where i.Reconstruction_Time_Ma == time && ((i.Land1 == land1 && i.Land2 == land2) || (i.Land2 == land1 && i.Land1 == land2))
        @select i.Closest_Distance_km
        @collect
    end

    if length(q) == 0
        error("STOP ERROR in pygplatesfunctions.distance_given().\n Pairing: '", land1, "' and '", land2, "' has not been found within given dataframe at this timestamp.\n Please check spelling of compared landmasses.\n If using function distance_given(), please check dataframe timestamps.\n \t if time intervals taken from original pygplates do not have your needed timestamp, please us function pygplatesfunctions.distance_interp() \n If using function distance_interp(), please ensure time intervals within dataframe go PAST the requested timestamp. Thank you! ")
    end


    if q[1] == "NA"
        q1 = @from i in df begin
            @where i.Reconstruction_Time_Ma == time && i.Land1 == land1 && i.Land2 == land2
            @select i.Point1_Lat
            @collect
        end

        q2 = @from i in df begin
            @where i.Reconstruction_Time_Ma == time && i.Land1 == land1 && i.Land2 == land2
            @select i.Point2_Lat
            @collect
        end
        
        if length(q1) == 0
            q3 = @from i in df begin
                @where i.Reconstruction_Time_Ma == time && i.Land1 == land2 && i.Land2 == land1
                @select i.Point1_Lat
                @collect
            end

            q4 = @from i in df begin
                @where i.Reconstruction_Time_Ma == time && i.Land2 == land1 && i.Land1 == land2 
                @select i.Point2_Lat
                @collect
            end

            if q4[1]=="NA" && q3[1]!="NA"
                error("STOP ERROR in pygplatesfunctions.distance_given(). The landmass, '", land1, "',\n does not have a polygon at this timestamp (time=", time, ").\n Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.")
            end # end error check

            if q3[1]=="NA" && q4[1]!="NA"
                error("STOP ERROR in pygplatesfunctions.distance_given(). The landmass, '", land2, "',\n does not have a polygon at this timestamp (time=", time, ").\n Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.")
            end # end error check

            if q3[1]=="NA" && q4[1]=="NA"
                error("STOP ERROR in pygplatesfunctions.distance_given(). Both landmasses, '", land1, "' and '", land2, "',\n do not have polygons at this timestamp\n (time=", time, ").\n Please use land_begin to find the earliest timestamp with polygons present. If land_begin returns an earlier timestamp, please check gplates file.")
            end
        end


        if q1[1]=="NA" && q2[1]!="NA"
            error("STOP ERROR in pygplatesfunctions.distance_given(). The landmass, '", land1, "',\n does not have a polygon at this timestamp (time=", time, ").\n Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.")
        end # end error check

        if q2[1]=="NA" && q1[1]!="NA"
            error("STOP ERROR in pygplatesfunctions.distance_given(). The landmass, '", land2, "',\n does not have a polygon at this timestamp (time=", time, ").\n Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.")
        end # end error check

        if q1[1]=="NA" && q2[1]=="NA"
            error("STOP ERROR in pygplatesfunctions.distance_given(). Both landmasses, '", land1, "' and '", land2, "',\n do not have polygons at this timestamp \n (time=", time, ").\n Please use land_begin to find the earliest timestamp with polygons present. If land_begin returns an earlier timestamp, please check gplates file.")
        end

    end 


    distance_km = q[1]
    return distance_km

end


"""
    distance_interp(df, time, land1, land2)

Provides the distance at a given time point, to be used when timepoint requested is
not within the dataframe.

To pull up distance timestamps already recorded within the dataframe, a faster function 
would be distance_given

* `df` - dataframe created by pygplates output created by Wallis Bland
contains variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

* `time` - timepoint for distance requested. May be any given time as long as the recorded
timespan within the dataframe is extends beyond requested timepoint.

* `land1` - first landmass chosen for distance comparison. Ensure this is within 'string' format

* `land2` - second landmass chosed for distance comparison. Ensure this is within 'string' format

NOTE: Lands 1 and 2 do not need to line up with Land1 or 2 within the dataframe in this function, as it will search both!

Process: Takes the two closest timestamps recorded in dataframe  produces a weighed average

# Examples
```julia-repl


julia> test_df = DataFrame(Reconstruction_Time_Ma=[0, 10, 20, 30, 40, 50, 60, 0, 10, 20, 30, 40, 50, 60], Land1=["India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India"], Land2=["Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar","Madagascar"], Closest_Distance_km=[0, 1000, 2000, 3500, 4100, 4900, "NA", 0, 1000, 2000, 3500, 4100, 4900, "NA"], Point1_Lat=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point1_Lon=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point2_Lat=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"], Point2_Lon=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"])

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1   Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Int64                   String  String      Any                  Any         Any         Any         Any        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │                      0  India   Asia        0                    1           1           1           1
   2 │                     10  India   Asia        1000                 2           2           2           2
   3 │                     20  India   Asia        2000                 3           3           3           3
   4 │                     30  India   Asia        3500                 4           4           4           4
   5 │                     40  India   Asia        4100                 5           5           5           5
   6 │                     50  India   Asia        4900                 6           6           6           6
   7 │                     60  India   Asia        NA                   NA          NA          7           7
   8 │                      0  India   Madagascar  0                    1           1           1           1
   9 │                     10  India   Madagascar  1000                 2           2           2           2
  10 │                     20  India   Madagascar  2000                 3           3           3           3
  11 │                     30  India   Madagascar  3500                 4           4           4           4
  12 │                     40  India   Madagascar  4100                 5           5           5           5
  13 │                     50  India   Madagascar  4900                 6           6           6           6
  14 │                     60  India   Madagascar  NA                   NA          NA          NA          NA

julia> distance_km = distance_interp(test_df, 44, "India", "Asia")

DataValue{Float64}(4420.0)


# POTENTIAL ERRORS

* INCORRECT VARIABLES / MISSPELLINGS
    In this case the variable was either spelt wrong or for some other reason Variable was never printed within the given dataframe.

    julia> distance_interp(test_df, 44, "Asia", "Indai")
    ERROR: STOP ERROR in pygplatesfunctions.distance_interp().
     Pairing: 'Asia' and 'Indai' has not been found within given dataframe at this timestamp.



* INCORRECT PAIRING
    Lands exist within the dataframe, but were never compared distance wise 
    If original pygplates code was run correctly, this SHOULD NOT OCCUR. More likely need to check spelling

    julia> distance_interp(test_df, 44, "Asia", "Madagascar")
    ERROR: STOP ERROR in pygplatesfunctions.distance_interp().
     Pairing: 'Asia' and 'Madagascar' has not been found within given dataframe.



* TIME INVALID
    If the timestamp requested is beyond the reconstructed years (in example our further back year is 60, we requested 64)
    there is no higher time for the interpolator to use as a weighted average. 
    Return to pygplates python code to change the output years to cover the correct range.

    julia> distance_km = distance_interp(test_df, 64, "India", "Asia")
    ERROR: STOP ERROR in pygplatesfunctions.distance_interp().
     Timestamp: 64 is beyond the reconstructed time in the given dataframe



* NAs within Outputs
    These errors will appear as errors within the pygrplatesfunction.distance_given() code.
    When the selected higher and lower time are sent through the distance_given function to retrieve the distances for
    the weighted average, they will return an error stating which land has failed at which time stamp.

    This is most likely to occur if there are missing polygons within outputs, especially when polygons appear and disappear.
    Polygon location and 'closest distance' can be set to 0s before they appear. However that is up to the user to decide.

    julia> distance_interp(test_df, 56, "India", "Asia") 
    ERROR: STOP ERROR in pygplatesfunctions.distance_given(). The landmass, 'India', 
     does not have a polygon at this timestamp (time=60).
     Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.
    

```
"""


# HEY SHOULD THIS SKIP ONES THAT HAVE NA DISTANCES?

function distance_interp(df, time, land1, land2)
    
    List_of_times = @from i in df begin
        @where ((i.Land1 == land1 && i.Land2 == land2) || (i.Land2 == land1 && i.Land1 == land2))
        @select i.Reconstruction_Time_Ma
        @collect
    end


    # should this actually be
    """
    List_of_times = @from i in df begin
        @where i.Closest_Distance_km != "NA" && ((i.Land1 == land1 && i.Land2 == land2) || (i.Land2 == land1 && i.Land1 == land2))
        @select i.Reconstruction_Time_Ma
        @collect
    end

    """

    if (length(List_of_times) == 0)
        error("STOP ERROR in pygplatesfunctions.distance_interp().\n Pairing '", land1, "' and '", land2, "' has not been found within given dataframe")
    end

    if time - last(List_of_times) > 0
        error("STOP ERROR in pygplatesfunctions.distance_interp().\n Timestamp: ", time, " (mya) is beyond the reconstructed time in the given dataframe")
    end

    i = 1
    global time_low, time_high
    for timecheck in List_of_times

        if timecheck == "NA"
            continue
        end

        if (time - timecheck <= 0)
            time_low = List_of_times[i-1]
            time_high = List_of_times[i]

            break
        end
        
        i = i+1
    end
    
    dist_high = distance_given(df, time_high, land1, land2)
    dist_low = distance_given(df, time_low, land1, land2)
    
    t2b_high = abs(time - time_high)
    t2b_low = abs(time - time_low)
    t2bsum = t2b_high + t2b_low
    
    weight_high = (t2bsum - t2b_high)/t2bsum
    weight_low = (t2bsum - t2b_low)/t2bsum
    
    weight_dist_high = weight_high * dist_high
    weight_dist_low = weight_low * dist_low
    
    distance_km =  weight_dist_high + weight_dist_low
    return distance_km

end


"""
take the two closest timestamps

by tens, you want 44
you would take 40 & 50 (array of those two distances)

if time_low > time_low_save && time_low < time
    time_low_save = time_low
end

if time_high < time_high_save && time_high > time
    time_high_save = time_high
end

weighted average?
borders times: 40 & 50 in this case 
distances at border timestamps: 1000 & 1300
time to borders: 4 & 6 (time_high_save - time & time - time_low_save)

weights: ((sum of d2b) - d2b)/(sum of d2b)
weights by original border distance
 then sum those

we want 44

q1 = @from i in df begin
    @where i.Land1 == land1 && i.Land2 == land2
    @select i.time
    @collect
end


time_dists = abs.(q1 .- time)
t2b_low = min(time_dists)
t2b_high = 

TF = q1 .== time_dists
rownums =  collect(1:length(q1))
time_low = q1[TF][1]
time_high = q1[rownums[TF][1] + 1]



q1 = @from i in df begin
    @where i.Land1 == land1 && i.Land2 == land2
    @select i.time
    @collect
end

for i in q1
    time_low = i
    time_high = i
    
    if time_low > time_low_save && time_low <= time
        time_low_save = time_low
    end
    
    if time_high < time_high_save && time_high >= time
        time_high_save = time_high
    end
end

List_of_times = @from i in df begin
    @where i.Land1 == land1 && i.Land2 == land2
    @select i.time
    @collect
end

for timecheck in List_of_times
    if time - timecheck <= 0
        time_low = q1[timecheck-1]
        time_high = q1[timecheck]
        break
    end
end

distance_high = distance(df, time_high, land1, land2)
distance_low = distance(df, time_low, land1, land2)

t2b_high = abs(time - time_high)
t2b_low = abs(time - time_low)
t2bsum = t2b_high + t2b_low

weight_high = (t2bsum - t2b_high)/t2bsum
weight_low = (t2bsum - t2b_low)/t2bsum

weight_dist_high = weight_high * distance_high
weight_dist_low = weight_low * distance_low

distance_km =  weight_dist_high + weight_dist_low

"""

"""
    oldest_polygon(df, land)

Provides the last timestamp in which a polygon can be found. Some gpml files have polygons that disappear and reappear!

* `df` - dataframe created by pygplates output created by Wallis Bland
contains variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

* `time` - timepoint for distance requested. May be any given time as long as the recorded
timespan within the dataframe is extends beyond requested timepoint.

* `land1` - first landmass chosen for distance comparison. Ensure this is within 'string' format

* `land2` - second landmass chosed for distance comparison. Ensure this is within 'string' format


# Examples
```julia-repl


julia> test_df = DataFrame(Reconstruction_Time_Ma=[0, 10, 20, 30, 40, 50, 60, 0, 10, 20, 30, 40, 50, 60], Land1=["India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India"], Land2=["Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar","Madagascar"], Closest_Distance_km=[0, 1000, 2000, 3500, 4100, 4900, "NA", 0, 1000, 2000, 3500, 4100, 4900, "NA"], Point1_Lat=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point1_Lon=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point2_Lat=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"], Point2_Lon=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"])

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1   Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Int64                   String  String      Any                  Any         Any         Any         Any        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │                      0  India   Asia        0                    1           1           1           1
   2 │                     10  India   Asia        1000                 2           2           2           2
   3 │                     20  India   Asia        2000                 3           3           3           3
   4 │                     30  India   Asia        3500                 4           4           4           4
   5 │                     40  India   Asia        4100                 5           5           5           5
   6 │                     50  India   Asia        4900                 6           6           6           6
   7 │                     60  India   Asia        NA                   NA          NA          7           7
   8 │                      0  India   Madagascar  0                    1           1           1           1
   9 │                     10  India   Madagascar  1000                 2           2           2           2
  10 │                     20  India   Madagascar  2000                 3           3           3           3
  11 │                     30  India   Madagascar  3500                 4           4           4           4
  12 │                     40  India   Madagascar  4100                 5           5           5           5
  13 │                     50  India   Madagascar  4900                 6           6           6           6
  14 │                     60  India   Madagascar  NA                   NA          NA          NA          NA

  julia> oldest_polygon(test_df, "Asia")
  oldest_time_stamp_with_polygon = 60

  julia> oldest_polygon(test_df, "India")
  oldest_time_stamp_with_polygon = 50

  julia> oldest_polygon(test_df, "Indai")
  ERROR: STOP ERROR in pygplatesfunctions.oldest_polygon.
   Landmass 'Indai' either does not exist in dataframe or does not contain any polygons.

"""


function oldest_polygon(df, land)
    
    List_of_times = @from i in df begin
        @where i.Point1_Lat != "NA" && i.Land1 == land
        @select i.Reconstruction_Time_Ma
        @collect
    end

    if length(List_of_times) == 0
        List_of_times = @from i in df begin
            @where i.Point2_Lat != "NA" && i.Land2 == land
            @select i.Reconstruction_Time_Ma
            @collect
        end
    end

    if length(List_of_times) == 0
        error("STOP ERROR in pygplatesfunctions.oldest_polygon.\n Landmass '", land, "' either does not exist in dataframe or does not contain any polygons.")
    end

    global oldest_time_stamp_with_polygon
    oldest_time_stamp_with_polygon = last(List_of_times)
    print("oldest_time_stamp_with_polygon = ", oldest_time_stamp_with_polygon)

end




function land_last_touch(df, land1, land2)
    distance = 0    
    
    q = @from i in df begin
        @where i.Closest_Distance_km == 0 && i.Land1 == land1 && i.Land2 == land2
        @select i.Reconstruction_Time_Ma
        
        @collect
   end

    
    if (length(areas_list) != numareas)
        txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). In your input file,\n'", lgdata_fn, "',\nthe given number of areas in line 1 (numareas=", numareas, ") does not equal the length of the areas_list (", parts2, ").\nPlease correct your input geography file  retry. Have a nice day."])
        error(txt)
    end # end error check


end



"""
# Goal: Read in various text files,
# starting with the output from wallis' pygplates code
# 
# looks like:

# =================
# Reconstruction_Time_Ma	Land1	Land2	Closest_Distance_km	Point1_Lat	Point1_Lon	Point2_Lat	Point2_Lon
# 0	   India	Africa	1995.3783881556299	23.58040000000074	66.3544	12.660000510000089	51.33756291700017
# 1	   India	Africa	1984.115611081836	23.163004686839308	66.22018058892301	12.473779939253923	51.1800340122845
# 2    India	Africa	1972.858728890917	22.746898758675545	66.08224988394925	12.287658928675025	51.02252712595236
# ...
# =================
#
# This is basically a normal dataframe?
# Should look into dataframes.jl for ideas?


include("/GitHub/PhyBEARS.jl/notes/Parsers.jl")
import .Parsers
pygp_data = "Desktop/School/Rhacophoridae/Gplates/output.txt"
pygp_df = Pygplates_parser.txtdf_read
"""


"""
Function to read in Lagrange/BioGeoBEARS-type PHYLIP-formatted
 geography input files.

pygp_data = lagrange-style geography input file
block_allQs = give error if a row has all-question-marks
"""



function txtdf_read(pygp_data; block_allQs=true)

	# Obtain a file handle
	fhandle = open(pygp_data)

	# Read the lines to an array
	lines = readlines(fhandle)
	close(fhandle)
	
	# Count the species read in
	compare = 0
	
	# We have to assume that the first line is a header
	for i in 1:length(lines)
		# If the line is blank, skip it
		line = lines[i]
		if (length(line) == 0)
			continue
		end
		#println(i, ": ", line)

		# Parse the header line
		if (i == 1)
			# Split on "("
			parts = split(line, "(")
			# remove whitespace, convert to Integer using parse.()
			parts1 = strip(parts[1])
			global (numtaxa, numareas) = parse.(Int, split(strip(parts1)))

			# Part 2: remove trailing ")", parse the rest to areas list
			parts2 = strip(strip(strip(parts[2]), [')']))
			# Split on whitespace to get the list of areas
			global areas_list = split(parts2)

			# Check that the given number matches the actual number of areas
			if (length(areas_list) != numareas)
				txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). In your input file,\n'", lgdata_fn, "',\nthe given number of areas in line 1 (numareas=", numareas, ") does not equal the length of the areas_list (", parts2, ").\nPlease correct your input geography file  retry. Have a nice day."])
				error(txt)
			end # end error check
	 
			# Set up the output matrix
			global sp_names = collect(repeat([""], numtaxa))
			tmp_areas = collect(repeat(["0"], numareas))
			global sp_areas = collect(repeat([tmp_areas], numtaxa))

			continue # Go to next round of for-loop without finishing this one
		end # END if (i == 1)

		# All of the other lines, fill in sp_names  sp_areas
		words = strip.(split(line))
		compare = compare + 1
	
		# Error check
		if (spnum > numtaxa)
			txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). While reading in species #", spnum, ", you exceeded the limit declared in the first line of your input geography file, where numtaxa=", numtaxa, " species. Please correct your input geography file  retry. Have a nice day."])
			error(txt)
		end

		
		sp_names[spnum] = words[1]
		# Split areanums to "1", "0", etc., parse to Integers
		#areanums_for_this_species = parse.(Int, split(words[2], ""))
		areanums_for_this_species = split(words[2], "")
		#print(areanums_for_this_species)
		# [:] avoids creating a linked reference
		sp_areas[spnum] = areanums_for_this_species[:]

	end

	# Error check
	if (spnum != numtaxa)
		txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). While reading in species, only ", spnum, " taxa were read in, however the first line of your input geography file declared there would be numtaxa=", numtaxa, " species. Please correct your input geography file  retry. Have a nice day."])
		error(txt)
	end

	# DataFrame of the area presence/absences
	# Flatten to vector, reshape to array, transpose to get
	# back original inputs in matrix form
	sp_areas2 = permutedims(reshape(flat2(sp_areas), numareas, numtaxa))
	sp_areas2_Matrix = convert(Matrix, sp_areas2)
	sp_areas3_Matrix = Rcbind(sp_names, sp_areas2_Matrix)
	geog_df = convert(DataFrame, sp_areas3_Matrix)
	new_names = Rrbind(["tipnames"], areas_list)
	geog_df = rename(geog_df, new_names)

	geog_df[!,:K]
	geog_df[!,:tipnames]

	return geog_df
end # END getranges_from_LagrangePHYLIP


#######################################################
# Put the tip ranges into the likelihoods
#######################################################
function tipranges_to_tiplikes(inputs, geog_df)
	dfnames = names(geog_df)
	area_column_nums = 2:length(dfnames)
	areas_txt_list = dfnames[area_column_nums]
	numareas = length(inputs.setup.areas_list)
	
	# Check if the number of areas in the geography file matches the number in the geog_df
	if (inputs.setup.numareas != numareas)
		txt = paste0(["STOP ERROR in tipranges_to_tiplikes(): inputs.setup.numareas=", numareas, ", but the number of areas in geog_df is ", numareas, ". Please fix  re-run."])
		error(txt)
	end
	
# 	statenums = collect(1:length(states_list))
# 	observed_statenums = collect(repeat([0], nrow(geog_df)))
	trdf_nodenums = collect(1:nrow(inputs.trdf))

	# Convert the observed geography into a list of states comparable to states_list
	geog_observed_list = []  # empty array

	for i in 1:nrow(geog_df)
		tmprow = geog_df[i,area_column_nums]
		tmpnums = parse.(Int, flat2(tmprow))
		range_as_areanums = inputs.setup.areas_list[tmpnums .== 1]
		# Compare observed range_as_areanums to full states_list
		TF = [range_as_areanums] .== inputs.setup.states_list
		if (sum(TF) != 1)
			txt = paste0(["STOP ERROR: An observed range in your geography file, from tipname '", geog_df[i,:tipnames], "', is not found in the list of states in inputs.setup.states_list. Printing range_as_areanums, then inputs.setup.states_list"])
			print(txt)
			print("\nrange_as_areanums (this is the observed range that was not found):\n")
			print(range_as_areanums)
			print("\nstates_list:\n")
			print(inputs.setup.states_list)
			error(txt)
		end
	
		# Yay, a single match to the states_list was found!
		# Convert to a state index number
		inputs.setup.observed_statenums[i] = inputs.setup.statenums[TF][1]
	end

	# Go through the geography file tips, match each one to the correct node of trdf,
	# then update the tip likelihoods at that node.

	for i in 1:nrow(geog_df)
		spname = geog_df[i,:tipnames]
		TF = spname .== inputs.trdf[!,:nodeName]
		nodeNum = trdf_nodenums[TF][1]
	
		# Input likelihoods of 1 for the observed state, 0 otherwise
		inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0         # zero out
		inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0     # zero out
		inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
		inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
		inputs.res.sumLikes_at_node_at_branchTop[nodeNum] = 1.0
	end

	return inputs
end # END function tipranges_to_tiplikes ()


end # ENDING Parsers
