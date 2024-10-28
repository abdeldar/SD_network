
module SDmoduleInit
    
    using IntervalTrees, DataFrames

    export Initialize_SD_object, Initialize_pyr_object, Interval_object, Itree_from_DF, Itree_from_DF_ape_lineage
    
    mutable struct Interval_object
        int_DF::DataFrame
        int_tree::Dict{String, IntervalTrees.IntervalBTree}
    end
    
    function Read_SDs_file_to_DF(filename::String, comp::Int = 1)
        SD_df = DataFrame(chr = String[], coor_s = Int64[], coor_e = Int64[],
            chr2 = String[], coor_s2 = Int64[], coor_e2 = Int64[], identity=Float64[])
        open(filename) do sdfile
            lines = readlines(sdfile);
            if comp == 1
                coors_str = map(x -> strip(x) |> split |> x -> x[[1,2,3,7,8,9,26]], lines)
            elseif comp == 2
                coors_str = map(x -> strip(x) |> split |> x -> x[[1,2,3,4,5,6,22]], lines)
            elseif comp == 3
                coors_str = map(x -> strip(x) |> split |> x -> x[[1,2,3,5,6,7,12]], lines)
            end
            map(x -> push!(SD_df ,[replace(x[1], r"chr" => "") [parse.(Int64, x[2:3])...]' replace(x[4], r"chr" => "") [parse.(Int64, x[5:6])...]' parse(Float64, x[7])]), coors_str)
        end
        return SD_df
    end

    function Filter_SDs_dataframe(SD_df::DataFrame, comp::String = "left")
        SD_df_filt = similar(SD_df, 0)
        for line in eachrow(SD_df)
            chr1 = line[:chr]
            chr2 = line[:chr2]
            if comp == "left"
                if chr1 < chr2
                    push!(SD_df_filt, line)
                elseif (chr1 == chr2) & (line[:coor_s] < line[:coor_s2])
                    push!(SD_df_filt, line)
                end
            elseif comp == "both"
                push!(SD_df_filt, line)
            end
        end
        SD_df_filt[!, :ids] = map(x -> "id"*string(x), 1:size(SD_df_filt)[1])
        return SD_df_filt
    end
    
    function Itree_from_DF(DF::DataFrame, comp::String="SD_IntTree", gap_value::Int=5)
        Itree::Dict{String, IntervalTrees.IntervalBTree} = Dict()
        if comp == "SD_IntTree"
            for line in eachrow(DF)
                Itree[line[:chr]] = get(Itree, line[:chr], IntervalTree{Int, IntervalValue{Int, String}}())
                push!(Itree[line[:chr]], IntervalValue(line[:coor_s], line[:coor_e], line[:ids]))
                Itree[line[:chr2]] = get(Itree, line[:chr2], IntervalTree{Int, IntervalValue{Int, String}}())
                push!(Itree[line[:chr2]], IntervalValue(line[:coor_s2], line[:coor_e2], line[:ids]))
            end
        elseif comp == "Pyrs_IntTree"
            for line in eachrow(DF)
                Itree = Itree_merge_overlaps(Itree, line[:chr], line[:coor_s], line[:coor_e], gap_value)
                Itree = Itree_merge_overlaps(Itree, line[:chr2], line[:coor_s2], line[:coor_e2], gap_value)
            end        
        end
        return Itree
    end
                
    function Check_sharp_borders_to_ids(SD_df::DataFrame, SD_tree::Dict, comp::String, gap::Int64, gap2::Int64)
        ind_arr = ones(size(SD_df)[1])
        for (i, line) in enumerate(eachrow(SD_df))
            if comp == "first"
                chr, coor_s, coor_e = line[:chr], line[:coor_s], line[:coor_e]
            elseif comp == "second"
                chr, coor_s, coor_e = line[:chr2], line[:coor_s2], line[:coor_e2]
            end
            over = collect(intersect(SD_tree[chr], (coor_s - gap, coor_e + gap)))
            for inter in over
                if ((inter.first - gap <= coor_s) && (coor_s <= inter.first + gap)) || ((inter.last - gap <= coor_e) && (coor_e <= inter.last + gap))
                    if coor_e - coor_s + gap2 < inter.last - inter.first   # don't have to check the id with gap3
                        ind_arr[i] = 0
                    end
                end
                if ((inter.last - gap <= coor_s) && (coor_s <= inter.last + gap)) || ((inter.first - gap <= coor_e) && (coor_e <= inter.first + gap))
                    ind_arr[i] = 1
                    break
                end
            end
        end
        return ind_arr
    end
    
    function Filt_sharp_borders_SDs(SD_obj::Interval_object, gap::Int64=15, gap2::Int64=100)
        SD_df, SD_tree = SD_obj.int_DF, SD_obj.int_tree
        ind_arr = Check_sharp_borders_to_ids(SD_df, SD_tree, "first", gap, gap2)
        ind_arr2 = Check_sharp_borders_to_ids(SD_df, SD_tree, "second", gap, gap2)
        SD_df2 = SD_df[(ind_arr .== 1) .| (ind_arr2 .== 1), :]                         # in question &
        SD_tree2 = Itree_from_DF(SD_df2, "SD_IntTree")
        return Interval_object(SD_df2, SD_tree2), ind_arr, ind_arr2
    end
    
                
    function Itree_merge_overlaps(Itree::Dict, chr::String, coor_s::Int64, coor_e::Int64, gap_value::Int=5)
        Itree[chr] = get(Itree, chr, IntervalTree{Int, IntervalValue{Int, String}}())
        over = collect(intersect(Itree[chr], (coor_s - gap_value, coor_e + gap_value)))
        if length(over) == 0
            push!(Itree[chr], IntervalValue(coor_s, coor_e, ""))
        else
            left_b = minimum(push!(map(x -> x.first, over), coor_s))
            right_b = maximum(push!(map(x -> x.last, over), coor_e))
            map(x -> delete!(Itree[chr], (x.first, x.last)), over)
            push!(Itree[chr], IntervalValue(left_b, right_b, ""))
        end
        return Itree
    end
                    
    function DF_from_Itree(Itree::Dict)
        Pyrs_df = DataFrame(chr = String[], coor_s = Int64[], coor_e = Int64[], ids = String[])
        Itree2::Dict{String, IntervalTrees.IntervalBTree} = Dict()
        id_num = 0
        for chr in sort(collect(keys(Itree)))
            for inter in Itree[chr]
                id_num += 1
                push!(Pyrs_df, [chr inter.first inter.last "id"*string(id_num)])
                Itree2[chr] = get(Itree2, chr, IntervalTree{Int, IntervalValue{Int, String}}())
                push!(Itree2[chr], IntervalValue(inter.first, inter.last, "id"*string(id_num)))
            end
        end
        return Itree2, Pyrs_df
    end


    function Initialize_SD_object(filename::String, comp::Int=1, comp2::String="left"; comp3::String="all")
        DF_prel = Read_SDs_file_to_DF(filename, comp)
        DF_filt = Filter_SDs_dataframe(DF_prel, comp2)
        Itree_SDs = Itree_from_DF(DF_filt, "SD_IntTree")
        SD_obj = Interval_object(DF_filt, Itree_SDs)
        if comp3 == "filter"
            SD_obj = Filt_sharp_borders_SDs(SD_obj)
        end
        return SD_obj
    end
                    
    function Initialize_pyr_object(DF_filt::DataFrame; gap_value::Int=5)
        Itree_prel = Itree_from_DF(DF_filt, "Pyrs_IntTree", gap_value)
        Itree_pyrs, DF_pyrs = DF_from_Itree(Itree_prel)
        Pyr_obj = Interval_object(DF_pyrs, Itree_pyrs)                
        return Pyr_obj
    end


    function Itree_from_DF_ape_lineage(DF::DataFrame, column)
        Itree::Dict{String, IntervalTrees.IntervalBTree} = Dict()
        for line in eachrow(DF)
            Itree[line[:chr]] = get(Itree, line[:chr], IntervalTree{Int, IntervalValue{Int, Float64}}())
            if column != nothing
                push!(Itree[line[:chr]], IntervalValue(line[:coor_s], line[:coor_e], line[column]))
            else
                push!(Itree[line[:chr]], IntervalValue(line[:coor_s], line[:coor_e], 0.0))
            end
        end         
        return Itree
    end
    
end
                                                                    
###############################################################################################

module NetworkAnalysis
    
    using Graphs, MetaGraphs, IntervalTrees, DataFrames, Main.SDmoduleInit
    
    export Graph_pars_object, Initialize_SD_network
                
    mutable struct Graph_pars_object
        mg::MetaGraph
        self_loops::Dict
        edges_double::Dict
        edges_type::Dict
        edges_tandem::Dict
        edges_ident::Dict
    end
    
    function Extract_fields_from_Graph_object(graph_obj::Graph_pars_object)
        out_array = [getfield(graph_obj, field) for field in [:mg, :self_loops, :edges_double, :edges_type, :edges_tandem, :edges_ident]]
        return out_array
    end
                
    function Graph_from_DF(int_obj::Interval_object, SD_df::DataFrame, gap_value::Int)
        self_loops::Dict = Dict()
        edges_double::Dict = Dict()
        edges_type::Dict = Dict()
        edges_tandem::Dict = Dict()
        edges_ident::Dict = Dict()
        mg = MetaGraph(size(int_obj.int_DF)[1], 1.0)     # 1.0 - default weight

        map(x -> set_prop!(mg, x, :id, int_obj.int_DF[!, :ids][x]), 1:size(int_obj.int_DF)[1])
        set_indexing_prop!(mg, :id)
        for line in eachrow(SD_df)
            over1 = collect(intersect(int_obj.int_tree[line[:chr]], (line[:coor_s] - gap_value, line[:coor_e] + gap_value)))
            over2 = collect(intersect(int_obj.int_tree[line[:chr2]], (line[:coor_s2] - gap_value, line[:coor_e2] + gap_value)))
            if length(over1) != 1 || length(over2) != 1
                display("Initialize_graph_from_DF warning")
            end
            vert1 = mg[over1[1].value, :id]
            vert2 = mg[over2[1].value, :id]
            
            if (Edge(vert1, vert2) in edges(mg)) && (vert1 != vert2)
                edges_double[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))] = get(edges_double, (get_prop(mg, vert1, :id), get_prop(mg, vert2, :id)), 0) + 1
                #edges_ident[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))] = max(edges_ident[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))], line[:identity])
            elseif !(Edge(vert1, vert2) in edges(mg)) && (vert1 != vert2)
                add_edge!(mg, vert1, vert2)
                edges_ident[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))] = line[:identity]
            end
            
            if vert1 == vert2
                self_loops[get_prop(mg, vert1, :id)] = get(self_loops, get_prop(mg, vert1, :id), 0) + 1
            else
                if line[:chr] == line[:chr2]
                    edges_type[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))] = get(edges_type, (get_prop(mg, vert1, :id), get_prop(mg, vert2, :id)), 1)
                    if min(abs(line[:coor_s] - line[:coor_e2]), abs(line[:coor_s2] - line[:coor_e])) < 5e5
                       edges_tandem[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))] = get(edges_tandem, (get_prop(mg, vert1, :id), get_prop(mg, vert2, :id)), 1)
                    end                     
                else
                    edges_type[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))] = get(edges_type, (get_prop(mg, vert1, :id), get_prop(mg, vert2, :id)), 0)
                end
            end
        end
        graph_obj = Graph_pars_object(mg, self_loops, edges_double, edges_type, edges_tandem, edges_ident)
        return graph_obj 
    end
    
    function Filter_single_nodes_everywhere!(int_obj::Interval_object, graph_obj::Graph_pars_object)
        mg, self_loops, edges_double, edges_type, edges_tandem, edges_ident = Extract_fields_from_Graph_object(graph_obj)
        node_single = [e[1] for e in filter(x -> length(x) == 1, connected_components(mg))]
        id_single = map(x -> get_prop(mg, x, :id), node_single)
        map(x -> rem_prop!(mg, mg[x, :id], :id), id_single)
        mg = mg[filter_vertices(mg, :id)]
        for id in id_single
            line = int_obj.int_DF[int_obj.int_DF[!, :ids] .== id, :]
            over = collect(intersect(int_obj.int_tree[line[!, :chr][1]], (line[!, :coor_s][1], line[!, :coor_e][1])))
            if length(over) == 1
                delete!(int_obj.int_tree[line[!, :chr][1]], (over[1].first, over[1].last))
            elseif length(over) > 1
                display("Filter_single_nodes_everywhere warning")
            end
            int_obj.int_DF = int_obj.int_DF[int_obj.int_DF[!, :ids] .!= id, :]
            delete!(self_loops, id)
        end
        set_indexing_prop!(mg, :id)
        graph_obj2 = Graph_pars_object(mg, self_loops, edges_double, edges_type, edges_tandem, edges_ident)
        return graph_obj2
    end
                            
    function Initialize_SD_network(int_obj::Interval_object, SD_df::DataFrame; gap_value::Int=5)
        graph_obj = Graph_from_DF(int_obj, SD_df, gap_value);
        graph_obj = Filter_single_nodes_everywhere!(int_obj, graph_obj);
        mg, self_loops, edges_double, edges_type, edges_tandem = Extract_fields_from_Graph_object(graph_obj);
        intra_num = reduce(+, [edges_type[i] for i in keys(edges_type)]);
        inter_num = length(edges_type) - intra_num;

        display("Double edges = $(length(edges_double)) = $(length(edges_double)/ne(mg)) of all edges")
        display("Self loops = $(length(self_loops)) = $(length(self_loops)/nv(mg)) of all nodes (or $(404/(nv(mg)+404-259)) without filtering singleton nodes)")
        display("Tandem edges = $(length(edges_tandem)) = $(length(edges_tandem)/ne(mg)) of all edges and $(length(edges_tandem)/intra_num) of all intra")
        display("Intra = $intra_num and Inter = $inter_num edges")
                                    
        return graph_obj;
    end
    
end
                                                            
###########################################################################################################
                                                                                            
module MST_from_SD

    using DataFrames, Main.SDmoduleInit, Graphs, MetaGraphs, SparseArrays, Distributions, StatsBase
    
    export graph_nosec_from_graph, assign_weights_edges!


    function assign_weights_edges!(graph::MetaGraph; comp="normal", suspicious_edges=nothing)
        cc = connected_components(graph)
        cl = map(length, cc)
        bws = betweenness_centrality(graph)
        if (comp == "sharp_out") & (suspicious_edges == nothing)
            display("suspicious_edges file is required. No output")
            return
        end
        for j in edges(graph)
            set_prop!(graph, j, :weight, 1.0)
        end
        for i in vertices(graph)
            neis = neighbors(graph, i)
            for n in neis
                sum_bws = sum(bws[neis])
                if sum_bws == 0
                   sum_bws = 1 
                end
                neis2 = neighbors(graph, n)
                mat_inds = findall(x -> x in neis, neis2)
                if Graphs.weights(graph)[i, n] == 1.0
                    if comp == "normal"
                        wei = 1.001 - (length(mat_inds) + 1)/length(neis)
                    elseif (comp == "normal_noise") | (comp == "sharp_out")
                        wei = 1.001 - (length(mat_inds) + 1)/length(neis) + rand(Uniform(-0.00005, 0.00005))
                    elseif comp == "betweenness_weight"
                        wei = 1.001 - bws[n]/sum_bws
                    elseif comp == "shuffle"
                        wei = 1.001 + rand()
                    end
                    set_prop!(graph, i, n, :weight, wei)
                else
                    if comp == "normal"
                        wei = minimum([1.001 - (length(mat_inds) + 1)/length(neis), Graphs.weights(graph)[i, n]])
                    elseif comp == "normal_noise"
                        wei = minimum([1.001 - (length(mat_inds) + 1)/length(neis), Graphs.weights(graph)[i, n]]) + rand(Uniform(-0.00005, 0.00005))
                    elseif comp == "sharp_out"
                        if (Edge(i, n) in suspicious_edges) | (Edge(n, i) in suspicious_edges)
                            wei = 100.0
                        else
                            wei = minimum([1.001 - (length(mat_inds) + 1)/length(neis), Graphs.weights(graph)[i, n]]) + rand(Uniform(-0.00005, 0.00005))
                        end
                    elseif comp == "betweenness_weight"
                        wei = minimum([1.001 - bws[n]/sum_bws, Graphs.weights(graph)[i, n]])
                    elseif comp == "shuffle"
                        wei = minimum([1.001 + rand(), Graphs.weights(graph)[i, n]])
                    end
                    set_prop!(graph, i, n, :weight, wei)
                end
            end
        end
    end


    function custom_cmp(x::String)
        number_idx = findfirst(isdigit, x)
        str, num = SubString(x, 1, number_idx-1), SubString(x, number_idx, length(x))
        return str, parse(Int, num)
    end


    function graph_nosec_from_graph(graph::MetaGraph; comp="mst_krus", comp_vis=false)
        if comp == "mst_krus"
            edg_list = kruskal_mst(graph)
            gr1 = induced_subgraph(graph, edg_list)[1]
        elseif comp == "mst_prim"
            ids, n_deg_int = [], []
            ccs = connected_components(graph)
            for cc in ccs
                cc_gr = induced_subgraph(graph, cc)[1]
                edg_list = prim_mst(cc_gr)
                gr1 = induced_subgraph(cc_gr, edg_list)[1]
                emap = vertices(gr1)
                append!(ids, map(x -> get_prop(gr1, x, :id), emap))
                append!(n_deg_int, map(x -> length(neighbors(gr1, x)), emap))
            end
            n_deg_int = n_deg_int[sortperm(ids, rev=false, by=custom_cmp)]
            ids = sort(ids, by=custom_cmp)
            return n_deg_int, ids                   # only 2 outputs(!)
        elseif comp == "all"
            gr1 = graph
            edg_list = edges(graph)
        end
        emap = vertices(gr1)
        ids = map(x -> get_prop(gr1, x, :id), emap)
        n_deg_int = map(x -> length(neighbors(gr1, x)), emap)
        n_deg_int = n_deg_int[sortperm(ids, rev=false, by=custom_cmp)]
        ids = sort(ids, by=custom_cmp)
        if comp_vis
            n_deg = convert.(UInt64, n_deg_int);
            nd_h = Histograms.Histogram(n_deg, scale=:log);
            p = plot(nd_h, yscale=:log10, xscale=:log10, line=false, markershape=:auto)
            x = range(1, 30, length=10); y = 1.0*x.^(-3);
            p = plot!(x, y, label="a = -3")
            gui(p)
        end
        return n_deg_int, ids, edg_list
    end

end

##########################################################################################

module Annotate_SDs

    using DataFrames, Main.SDmoduleInit, Graphs, MetaGraphs
    
    export Read_genes_and_gaps_2_DF, construct_DF!

    function Read_genes_and_gaps_2_DF(filename::String)
        df = DataFrame(chr=String[], coor_s=Int64[], coor_e=Int64[])
        open(filename) do sdfile
            lines = readlines(sdfile);
            coors_str = map(x -> strip(x) |> split |> x -> x[[1,2,3]], lines)
            map(x -> push!(df ,[replace(x[1], r"chr" => "") [parse.(Int64, x[2:3])...]']), coors_str)
        end
        return df
    end


    function construct_DF!(df::DataFrame, segments_obj::Main.SDmoduleInit.Interval_object, column::Symbol; border::Int=-1, gap::Int=0, comp::String="mean", comp2::String="both")
        for line in eachrow(df)
            chr, coor_s, coor_e = line[:chr], line[:coor_s], line[:coor_e]
            if border == -1
                if column in [:genes, :cpgisl_in]
                    if chr in keys(segments_obj.int_tree)
                        over = collect(intersect(segments_obj.int_tree[chr], (coor_s, coor_e)))
                        for ov in over
                            line[column] += 1 
                        end
                    end
                elseif column in [:repli_in, :dnase_in, :repli_deriv, :repli_vari, :recomb_in]
                    if chr in keys(segments_obj.int_tree)
                        over = collect(intersect(segments_obj.int_tree[chr], (coor_s, coor_e)))
                        over = collect(map(x -> x.value, over))
                        if length(over) > 0  
                            if comp == "mean"
                                line[column] = round(mean(over), digits=3)
                            elseif comp == "derivative"
                                line[column] = maximum(over) - minimum(over)
                            elseif comp == "variance"
                                if length(over) > 1
                                    line[column] = round(std(over), digits=3)
                                end
                            end
                        end
                    end
                end
            else
                if column in [:gaps, :cpgisl_bor, :ctcf, :genes_bor]
                    if coor_s - gap <= coor_e + gap
                        if (comp2 == "both") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_s - gap - border, coor_s - gap)))
                            append!(over, collect(intersect(segments_obj.int_tree[chr], (coor_e + gap, coor_e + gap + border))))
                        elseif (comp2 == "left") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_s - gap - border, coor_s - gap)))
                        elseif (comp2 == "right") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_e + gap, coor_e + gap + border)))
                        elseif !(chr in keys(segments_obj.int_tree))
                            over = []
                        end
                    else
                        if (comp2 == "both") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_s - border, coor_s)))
                            append!(over, collect(intersect(segments_obj.int_tree[chr], (coor_e, coor_e + border))))
                        elseif (comp2 == "left") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_s - border, coor_s)))
                        elseif (comp2 == "right") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_e, coor_e + border)))
                        elseif !(chr in keys(segments_obj.int_tree))
                            over = []
                        end
                    end
                    for ov in over
                        line[column] += 1 
                    end
                elseif column in [:repli_bor, :repli_bor_deriv, :dnase_bor, :recomb_bor]
                    if coor_s - gap <= coor_e + gap
                        if (comp2 == "both") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_s - gap - border, coor_s - gap)))
                            append!(over, collect(intersect(segments_obj.int_tree[chr], (coor_e + gap, coor_e + gap + border))))
                        elseif (comp2 == "left") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_s - gap - border, coor_s - gap)))
                        elseif (comp2 == "right") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_e + gap, coor_e + gap + border)))
                        elseif !(chr in keys(segments_obj.int_tree))
                            over = []
                        end
                    else
                        if (comp2 == "both") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_s - border, coor_s)))
                            append!(over, collect(intersect(segments_obj.int_tree[chr], (coor_e, coor_e + border))))
                        elseif (comp2 == "left") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_s - border, coor_s)))
                        elseif (comp2 == "right") & (chr in keys(segments_obj.int_tree))
                            over = collect(intersect(segments_obj.int_tree[chr], (coor_e, coor_e + border)))
                        elseif !(chr in keys(segments_obj.int_tree))
                            over = []
                        end
                    end                                        
                    if (length(over) > 0) & (column == :repli_bor_deriv)
                        line[column] = round(maximum(collect(map(x -> x.value, over))), digits=3)
                    elseif length(over) > 0
                        line[column] = round(mean(collect(map(x -> x.value, over))), digits=3)
                    end
                end
            end
        end
    end

end


############################################################################################

module Script_funcs

    using Graphs, CSV, MetaGraphs, DataFrames, ArgParse, Main.MST_from_SD, GraphIO, Main.Annotate_SDs, Main.SDmoduleInit
    
    export handle_suspicious_scr, Run_MST_scr!, parse_arguments, save_all_dataframes

    function handle_suspicious_scr(b, graph_obj, graph_obj2, filt_sharp, keep_gaps, border; gaps_filename=nothing)
        if filt_sharp == 1
            dif = collect(edges(difference(graph_obj.mg, graph_obj2.mg)))
            if keep_gaps == 1
                if (gaps_filename != nothing) & isfile(gaps_filename)
                    gaps_df = Annotate_SDs.Read_genes_and_gaps_2_DF(gaps_filename)
                    gaps_tree = SDmoduleInit.Itree_from_DF_ape_lineage(gaps_df, nothing)
                    gaps_obj = SDmoduleInit.Interval_object(gaps_df, gaps_tree)
                    b.int_DF[!, :gaps] = zeros(Int, size(b.int_DF)[1])
                    Annotate_SDs.construct_DF!(b.int_DF, gaps_obj, :gaps, border=border, gap=0)
                    ids = b.int_DF[b.int_DF[!, "gaps"] .> 0.0, "ids"]
                    
                    suspicious_edges = []
                    for e in dif
                        if !((get_prop(graph_obj.mg, e.src, :id) in ids) | (get_prop(graph_obj.mg, e.dst, :id) in ids))
                            append!(suspicious_edges, [e])
                        end
                    end
                    println("Number of suspicious edges $(length(suspicious_edges)) (filtered out)")
                else
                    println(stderr, "ERROR: When keep_gaps == 1 assembly annotation file is needed. Can't find one.")
                    exit(-1)
                end
            elseif keep_gaps == 0
                suspicious_edges = deepcopy(dif)
                println("Number of suspicious edges $(length(suspicious_edges))")
            else
                println(stderr, "ERROR: Unknown keep_gaps value: $keep_gaps")
                exit(-1)
            end
        elseif filt_sharp == 0
            suspicious_edges = nothing
        else
            println(stderr, "ERROR: Unknown filt_sharp value: $filt_sharp")
            exit(-1)
        end
        return suspicious_edges
    end


    function Run_MST_scr!(b, graph_obj, suspicious_edges, N)
        b.int_DF[!, :jumps] = zeros(Float64, size(b.int_DF)[1])
        println("Do $N iterations....")
        for i in 1:N
            if suspicious_edges == nothing
                MST_from_SD.assign_weights_edges!(graph_obj.mg, comp="normal_noise")
            else
                MST_from_SD.assign_weights_edges!(graph_obj.mg, comp="sharp_out", suspicious_edges=suspicious_edges)
            end
            n_deg, ids_hui, edg_lis = graph_nosec_from_graph(graph_obj.mg, comp="mst_krus")
            global edg_lis
            for (i, val) in enumerate(ids_hui)
                b.int_DF[b.int_DF[!,:ids] .== val, :jumps] .+= n_deg[i]
            end
        end
        b.int_DF[!, :jumps] = round.(b.int_DF[!, :jumps]/N)
        return edg_lis
    end


    function parse_arguments()
        s = ArgParseSettings("The script finds MST of duplication events from the given list of SDs. It also builds the SD network. \n
        It can be run in 3 possible mods: \n
        1) without filtering any sharp edges (-f 0 -k 0); \n
        2) filtering all 'suspicious' edges with matching duplication breakpoints (-f 1 -k 0);\n
        3) filtering suspicious edges which are not proximal to assembly gaps (-f 1 -k 1).\n
        It is possible to generate an ensemble of MSTs from the input network. It leads to less biased prediction of number of duplications that each node has undergone. Use '-N' parameter to run several iterations. Mean rounded number of duplications then reported as an output + single MST of duplications.")
        @add_arg_table s begin
            "--input_SDs", "-i"
                help = "input SDs file"
                arg_type = String
                required = true
            "--outdir", "-o"
                help = "output directory"
                arg_type = String
                default = "./outputs/MST_out_scr"
            "--N", "-N"
                help = "number of MST iterations"
                arg_type = Int
                default = 1
            "--border", "-b"
                help = "Flanking windows length"
                arg_type = Int
                default = 50
            "--filt_sharp", "-f"
                help = "Filter sharp borders out (1/0)"
                arg_type = Int
                default = 0
            "--gaps_file", "-g"
                help = "input assembly gaps annotation"
                arg_type = String
                default = nothing
            "--keep_gaps", "-k"
                help = "Keep duplicated regions associated with gaps unfiltered (1/0)"
                arg_type = Int
                default = 0
        end
        args = parse_args(s)
        N = args["N"]
        border = args["border"]
        filt_sharp = args["filt_sharp"]
        keep_gaps = args["keep_gaps"]
        input_SDs = args["input_SDs"]
        gaps_file = args["gaps_file"]
        outdir = args["outdir"]
        return [N, border, filt_sharp, keep_gaps, input_SDs, gaps_file, outdir]
    end


    function save_all_dataframes(b, edg_lis, graph_obj, outdir)
        savegraph("$outdir/SD_network_gen.lg", graph_obj.mg)
        Dup_df = DataFrame(chr = String[], coor_s = Int64[], coor_e = Int64[], ids1 = String[], chr2 = String[], coor_s2 = Int64[], coor_e2 = Int64[], ids2 = String[], identity=Float64[])
        for e in edg_lis
            id1, id2 = get_prop(graph_obj.mg, e.src, :id), get_prop(graph_obj.mg, e.dst, :id)
            df = b.int_DF[b.int_DF[!, "ids"] .== id1, ["chr", "coor_s", "coor_e", "ids"]]
            l1 = collect(df[1, :])
            df = b.int_DF[b.int_DF[!, "ids"] .== id2, ["chr", "coor_s", "coor_e", "ids"]]
            l2 = collect(df[1, :])
            identity = graph_obj.edges_ident[(id1, id2)]
            push!(Dup_df, vcat(l1, l2, [identity]))
        end
        CSV.write("$outdir/Duplicated_df_jumps.tab", b.int_DF[!, ["chr", "coor_s", "coor_e", "ids", "jumps"]], delim="\t")
        CSV.write("$outdir/Primary_edges_dupregs.tab", Dup_df, delim="\t")
    end
        
end
                                                                                                