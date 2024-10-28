using DataFrames
using Graphs

include("./modules/modules_scr.jl")

N, border, filt_sharp, keep_gaps, input_SDs, gaps_file, outdir = Script_funcs.parse_arguments()

if !isfile(input_SDs)
    println(stderr, "ERROR: No input SDs file")
    exit(-1)
end

# Create interval object
a = SDmoduleInit.Initialize_SD_object(input_SDs);
b = SDmoduleInit.Initialize_pyr_object(a.int_DF);
c, ind_arr, ind_arr2 = SDmoduleInit.Initialize_SD_object(input_SDs, comp3="filter")

# Create graph object
nex = deepcopy(b)
graph_obj = NetworkAnalysis.Initialize_SD_network(b, a.int_DF)
println("\n")
graph_obj2 = NetworkAnalysis.Initialize_SD_network(nex, c.int_DF)
println("\nGraphs are built")

# Identify suspicious edges
suspicious_edges = Script_funcs.handle_suspicious_scr(b, graph_obj, graph_obj2, filt_sharp, keep_gaps, border; gaps_filename=gaps_file)

# Predict MST
edg_lis = Script_funcs.Run_MST_scr!(b, graph_obj, suspicious_edges, N)

# Save
Script_funcs.save_all_dataframes(b, edg_lis, graph_obj, outdir)

println("Done")


