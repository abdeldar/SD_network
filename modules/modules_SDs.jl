
module SDmoduleInit
    
    using IntervalTrees, DataFrames

    export Initialize_SD_object, Initialize_pyr_object, Interval_object, Itree_from_DF
    
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
                        
    function filter_SEDEF_DF!(a_df)
        a_df = a_df[a_df[!, :identity] .>= 0.9, :]
        inds = min.(a_df[!,:coor_e] - a_df[!,:coor_s], a_df[!,:coor_e2] - a_df[!,:coor_s2]) .>= 999
        a_df = a_df[inds, :]
        return a_df
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
    
end

#####################################################################################################

module NetworkAnalysis
    
    using LightGraphs, MetaGraphs, IntervalTrees, DataFrames, Main.SDmoduleInit, GraphPlot, StatsBase
    using Colors, Compose, Plots, Histograms, GLM, Statistics, HypothesisTests, Distributions, ArndtLabTheme
    using LaTeXStrings
    pyplot()
    theme(:arndtlab)
    
    export Graph_pars_object, Initialize_SD_network, Plot_SD_network, Node_degree_Comp_size_distr, Stat_tests_SD_net_features, Spliting_network_diff_ways
                
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
            
            if has_edge(mg, Edge(vert1, vert2)) && (vert1 != vert2)
                edges_double[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))] = get(edges_double, (get_prop(mg, vert1, :id), get_prop(mg, vert2, :id)), 0) + 1
                #edges_ident[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))] = max(edges_ident[(get_prop(mg, vert1, :id), get_prop(mg, vert2, :id))], line[:identity])
            elseif !has_edge(mg, Edge(vert1, vert2)) && (vert1 != vert2)
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
                            
    function Plot_SD_network(graph_obj::Graph_pars_object, comp::String = "normal")
        mg, self_loops, edges_double, edges_type, edges_tandem, edges_ident = Extract_fields_from_Graph_object(graph_obj)
        if comp == "normal"
            gp = gplot(mg)
            draw(PDF("./pics/SD_net.pdf", 120cm, 120cm), gp)
        elseif comp == "self_loops"
            verts = map(x -> mg[x, :id], collect(keys(self_loops)))
            membership = ones(Int, nv(mg))
            membership[verts] .+= 1
            nodecol = [colorant"lightseagreen", colorant"red"][membership]
            gp = gplot(mg, nodefillc=nodecol)
            draw(PDF("./pics/SD_net_self_loops.pdf", 120cm, 120cm), gp)
        elseif comp == "edges_double"
            edgs = map(x -> Edge(mg[x[1], :id], mg[x[2], :id]), collect(keys(edges_double)))
            membership = map(x -> x in edgs ? 2 : 1, collect(edges(mg)))
            edgecol = [colorant"grey", colorant"red"][membership]
            gp = gplot(mg, edgestrokec=edgecol)
            draw(PDF("./pics/SD_net_double_edges.pdf", 120cm, 120cm), gp)                                    
        elseif comp == "edges_tandem"
            edgs = map(x -> Edge(mg[x[1], :id], mg[x[2], :id]), collect(keys(edges_tandem)))
            membership = map(x -> x in edgs ? 2 : 1, collect(edges(mg)))
            edgecol = [colorant"grey", colorant"red"][membership]
            gp = gplot(mg, edgestrokec=edgecol)
            draw(PDF("./pics/SD_net_tandem_edges.pdf", 120cm, 120cm), gp)                                               elseif comp == "edges_intra"
            edgs = map(x -> edges_type[x] == 1 ? Edge(mg[x[1], :id], mg[x[2], :id]) : 0, collect(keys(edges_type)))
            edgs = edgs[edgs .!= 0]
            membership = map(x -> x in edgs ? 2 : 1, collect(edges(mg)))
            edgecol = [colorant"grey", colorant"red"][membership]
            gp = gplot(mg, edgestrokec=edgecol)
            draw(PDF("./pics/SD_net_intra_edges.pdf", 120cm, 120cm), gp)
        elseif comp == "both_selfloops_intra"
            edgs = map(x -> edges_type[x] == 1 ? Edge(mg[x[1], :id], mg[x[2], :id]) : 0, collect(keys(edges_type)))
            edgs = edgs[edgs .!= 0]
            membership = map(x -> x in edgs ? 2 : 1, collect(edges(mg)))
            edgecol = [colorant"grey", colorant"red"][membership]                                               
            verts = map(x -> mg[x, :id], collect(keys(self_loops)))
            membership = ones(Int, nv(mg))
            membership[verts] .+= 1
            nodecol = [colorant"lightseagreen", colorant"red"][membership]
            gp = gplot(mg, edgestrokec=edgecol, nodefillc=nodecol);                                            
            draw(PDF("./pics/SD_both_selfloops_intra.pdf", 120cm, 120cm), gp);                                
        end
                                
    end
    
    function Node_degree_Comp_size_distr(mgs::Array, legends::Array{String,1}, comp=nothing; cc_bin::Int=20, nd_bin::Int=20, thr::Int=5, max_x::Int64=13000)                                         
        colors = ArndtLabTheme.arndtlab_palette
        colors2 = copy(colors)
        colors2[4], colors2[5] = colors[5], colors[4]
        colors2[5], colors2[6] = colors[6], colors[4]
        markers = [:circle, :circle, :circle, :circle, :circle, :circle, :circle, :circle, :circle]  # fivth might be :xcross
        m_sizes = [6, 6, 6, 6, 6, 6, 6, 6, 6]   # fivth might be 5
        h = Histograms.Histogram([UInt64(30000)], n=cc_bin, scale=:log)
        cc_p, nd_p, rat_p, rat_m_p = plot(), plot(), plot(), plot()
        for (i, mg) in enumerate(mgs)
            ccs = collect(connected_components(mg))
            cc_len = map(length, ccs);
            ucc_len = convert(Array{UInt,1}, cc_len)
            cc_h = Histograms.Histogram(ucc_len, breaks=h.breaks, scale=:log);
            plot!(cc_p, cc_h, xscale=:log10, yscale=:log10, markershape=markers[i], color=colors2[i],
                    label=legends[i], line=false, markersize=m_sizes[i], xtickfont=Plots.font(13),
                    ytickfont=Plots.font(13), xlim=[1, max_x], legendfontsize=12);
            x = range(2, 1800, length=10); y = 1.5*x.^(-2.7);
            plot!(cc_p, x, y, label="", grid=false, color=colors2[6])
            if comp != nothing
                savefig(comp)
            end

            n_deg = map(x -> length(neighbors(mg, x)), collect(vertices(mg)));
            n_deg = convert.(Float64, n_deg);
            nd_h = Histograms.Histogram(n_deg, N=nd_bin, scale=:identity);
            plot!(nd_p, nd_h, xscale=:identity, yscale=:log10, markershape=markers[i], color=colors2[i], label=legends[i], grid=false, line=false, markersize=m_sizes[i], xtickfont=Plots.font(13), ytickfont=Plots.font(13));    # or legends=false
            #nd_p = plot(nd_h, yscale=:log10, xscale=:log10, markershape=:auto, legend=false, grid=false, line=false, xlim=[1,170]);
            x = range(1, 150, length=10); y = 0.01*exp.(-0.056x);
            plot!(nd_p, x, y, label="", grid=false, color=colors2[6], legendfontsize=12)
            #savefig("./pics/node_degree_lin_SD.pdf")
               
            edgs::Array{Int64,1} = []
            for cc in ccs
                push!(edgs, sum([length(neighbors(mg, node)) for node in cc])/2)
            end                                                                                          
            plot!(rat_p, cc_len, edgs, xscale=:log10, yscale=:log10, color=colors2[i], markersize=m_sizes[i], markershape=markers[i], line=false, label="", grid=false, xtickfont=Plots.font(13), ytickfont=Plots.font(13))
            b, k, rat_p = Regression_from_2_arrays(cc_len, edgs, ((cc_len .> thr)), rat_p, [1, 2000], comp="log", color=colors2[i])                                 
            m_edgs = [mean(edgs[cc_len .== vals]) for vals in sort(unique(cc_len))]
            m_cc_len = sort(unique(cc_len))
            plot!(rat_m_p, m_cc_len, m_edgs, xscale=:log10, yscale=:log10, color=colors2[i], markershape=markers[i], line=false, grid=false, label="", markersize=m_sizes[i], xtickfont=Plots.font(13), ytickfont=Plots.font(13))
            xs = range(1, 2000, length=10);
            ys = (10^b)*xs.^k;
            plot!(rat_m_p, xs, ys, label=legends[i]*" (α = $k)", color=colors2[i], legendfontsize=15, legend=:bottomright)
            #b, k, rat_m_p = Regression_from_2_arrays(m_cc_len, m_edgs, ((m_cc_len .> 5)), rat_m_p, [1, 2000], comp="log", color=colors2[i])
            #savefig("./pics/WORM_comps_nodges.pdf")
        end
        
        display(cc_p)
        #savefig(cc_p, "./pics/review_cc_flanks_distr.pdf")
        display(nd_p)
        #savefig(nd_p, "./pics/review_nd_flanks_distr.pdf")
        display(rat_p)
        gui(rat_m_p)
        #savefig(rat_m_p, "./pics/animals_m_n_e_distr.pdf")
    end 
    
    function Regression_from_2_arrays(arrX::Array, arrY::Array, cond::BitArray{1}, plotobj::Plots.Plot, xlim::Array; comp::String = "linear", color=:red)
        arrX = arrX[cond]; arrY = arrY[cond]
        if comp == "linear"
            data = DataFrame(X=arrX, Y=arrY)
        elseif comp == "log"
            data = DataFrame(X=log10.(arrX), Y=log10.(arrY))
        end
        slm = lm(@formula(Y ~ X), data)
        b, k = coef(slm)                                                   
        x = range(xlim[1], xlim[2], length=10);
        if comp == "linear"
           y = k*x .+ b;                                             
        elseif comp == "log"                                                
           y = (10^b)*x.^k;                                             
        end
        k_round = round(k, digits=2)
        plot!(plotobj, x, y, label="α = $k_round", color=color)
        display(k)
        return b, k_round, plotobj                            
    end
                                                            
    function Stat_tests_SD_net_features(graph_obj::Graph_pars_object, limits::Array{Int64,1}=[10, 1000])
        mg, self_loops, edges_double, edges_type, edges_tandem, edges_ident = Extract_fields_from_Graph_object(graph_obj)
        big_nodes::Array{String,1} = []
        other_nodes::Array{String,1} = []
        middle_ccs::Array{Array{Int64,1},1} = []
        middle_intra::Array{Int64, 1} = []
        middle_lens::Array{Int64, 1} = []
                                                                
        ccs = collect(connected_components(mg))
        bc_size = maximum(map(length, ccs))
        map(x -> length(x) == bc_size ? append!(big_nodes, Node_array_to_back_ids(mg, x)) : append!(other_nodes, Node_array_to_back_ids(mg, x)), ccs)
        
        f_plus, f_minus = 0, 0
        map(x -> x in big_nodes ? f_plus += 1 : f_minus += 1, collect(keys(self_loops)))
        pv_self = pvalue(FisherExactTest(f_plus, f_minus, length(big_nodes), length(other_nodes)))
        display([f_plus/length(big_nodes) f_minus/length(other_nodes)])
        display("Fisher test p-value for self_loops in big/other components = $pv_self")

        f_plus, f_minus = 0, 0
        arr_intra = collect(filter(x -> edges_type[x] == 1, collect(keys(edges_type))))
        map(x -> (x[1] in big_nodes) && (x[2] in big_nodes) ? f_plus += 1 : f_minus += 1, arr_intra)
        ne_big = sum(map(x -> neighbors(mg, mg[x, :id]) |> length, big_nodes))/2
        ne_other = ne(mg) - ne_big
        pv_intra = pvalue(FisherExactTest(f_plus, f_minus, Int(ne_big), Int(ne_other)))
        display([f_plus/Int(ne_big) f_minus/Int(ne_other)])
        display("Fisher test p-value for intra_edges in big/other components = $pv_intra")
        
        f_plus, f_minus = 0, 0
        map(x -> (x[1] in big_nodes) && (x[2] in big_nodes) ? f_plus += 1 : f_minus += 1, collect(keys(edges_double)))
        ne_big = sum(map(x -> neighbors(mg, mg[x, :id]) |> length, big_nodes))/2
        ne_other = ne(mg) - ne_big
        pv_double = pvalue(FisherExactTest(f_plus, f_minus, Int(ne_big), Int(ne_other)))
        display([f_plus/Int(ne_big) f_minus/Int(ne_other)])
        display("Fisher test p-value for double_edges in big/other components = $pv_double")
                                                                
        middle_ccs = filter(x -> (length(x) >= limits[1]) && (length(x) <= limits[2]) , ccs)
        for comp in middle_ccs
            comp_intra = 0
            gr = induced_subgraph(mg, comp)[1]
            set_indexing_prop!(gr, :id)
            push!(middle_lens, ne(gr))
            for edg in edges(gr)
                edge = (gr[edg.src, :id], gr[edg.dst, :id])
                comp_intra += edges_type[edge]        
            end
            push!(middle_intra, comp_intra)
        end
        chi_val = Calculate_Chi_square(middle_lens, middle_intra)
        simul_chi_vals = Simul_intra_distribution(middle_lens, sum(middle_intra))
        #display(middle_intra./middle_lens)
        num_more = count(x -> x >= chi_val, simul_chi_vals)
        display("There are $num_more simulations with bigger intra- fraction value than observed in middle size components out of $(length(simul_chi_vals))")
        
    end
                                                            
    function Node_array_to_back_ids(mg::MetaGraph, node_arr::Array)
        [mg[i, :id] for i in node_arr]                                                                
    end
                                                            
    function Calculate_Chi_square(ns::Array{Int64,1}, ks::Array{Int64,1})
        p = sum(ks)/sum(ns)
        sum(((p*ns .- ks).^2)./(p*ns))                       
    end
                                                            
    function Simul_intra_distribution(ns::Array{Int64,1}, num::Int64, iter::Int64=1000)
        vals::Array{Float64,1} = []
        prs = ns/sum(ns)
        distr = Multinomial(num, prs)
        for i in 1:iter
            push!(vals, Calculate_Chi_square(ns, rand(distr)))
        end
        return vals
    end                                                           
                                                            
    function Spliting_network_diff_ways(graph_obj::Graph_pars_object, comp::String)
        mg, self_loops, edges_double, edges_type, edges_tandem, edges_ident = Extract_fields_from_Graph_object(graph_obj)
        cc = connected_components(mg)
        lens = map(length, cc)
        big_nodes = Node_array_to_back_ids(mg, cc[lens .== maximum(lens)][1])
        if comp == "intra" || comp == "intra_nobig"
            arr_intra = collect(filter(x -> edges_type[x] == 1, collect(keys(edges_type))))
            edges_intra = map(x -> Edge(mg[x[1], :id], mg[x[2], :id]), arr_intra)
            gr = induced_subgraph(mg, edges_intra)[1]                                                        
            set_indexing_prop!(gr, :id)                                                        
            if comp == "intra_nobig"
                node_arr = filter(ver -> !(gr[ver, :id] in big_nodes), vertices(gr))
                gr = induced_subgraph(gr, node_arr)[1]
                set_indexing_prop!(gr, :id)
            end
            return gr
        elseif comp == "inter" || comp == "inter_nobig"
            arr_inter = collect(filter(x -> edges_type[x] == 0, collect(keys(edges_type))))
            edges_inter = map(x -> Edge(mg[x[1], :id], mg[x[2], :id]), arr_inter)
            gr = induced_subgraph(mg, edges_inter)[1]                                                        
            set_indexing_prop!(gr, :id)
            if comp == "inter_nobig"
                node_arr = filter(ver -> !(gr[ver, :id] in big_nodes), vertices(gr))
                gr = induced_subgraph(gr, node_arr)[1]
                set_indexing_prop!(gr, :id)
            end
            return gr
        end                                                        
    end
end
                                                            
###########################################################################################################
                                                           
module PyrSDTypes
    
    using StatsBase, Statistics, DataFrames, Histograms, Plots, MetaGraphs, LightGraphs, StatsPlots, GLM, Main.SDmoduleInit, Main.NetworkAnalysis, DecisionTree, Random, ArndtLabTheme
    pyplot()
    theme(:arndtlab)
                                                            
    export Real_SD_length_distr, Pyr_length_node_degree!, Modify_DF_for_RandForset!
    
    function Real_SD_length_distr(Itree_SDs::Dict, Itree_pyrs::Dict, sds_df::DataFrame)
        real_SDs::Array{String, 1} = []
        for chr in keys(Itree_SDs)
            sing_SD = filter(x -> length(collect(intersect(Itree_SDs[chr], x))) == 1, collect(Itree_SDs[chr]))
            append!(real_SDs, [i.value for i in sing_SD])
        end
        real_SDs_df = sds_df[map(x -> x in real_SDs, sds_df[!, :ids]),:]
        #real_SDs_df = real_SDs_df[0.98 .< real_SDs_df[!, :identity] .< 0.9999, :]         
        lens::Array{Float64,1} = real_SDs_df[!, :coor_e] - real_SDs_df[!, :coor_s]
        cc_h = Histograms.Histogram(lens, N=30, scale=:log10);
        cc_p = plot(cc_h, xscale=:log10, yscale=:log10, markershape=:auto, line=false);
        x = range(10^3, 10^5.5, length=10); y = 900*x.^(-2);
        cc_p = plot!(x, y, label="a = -2")
        display(mean(lens))
        display(cc_p)
    end
                                                             
    function Pyr_length_node_degree!(pyrs_df::DataFrame, mg::MetaGraph)
        node_len::Dict = Dict()
        pyrs_df[!, :degree] = zeros(size(pyrs_df)[1])
        pyrs_df[!, :length] = convert.(Float64, pyrs_df[!, :coor_e] .- pyrs_df[!, :coor_s])
        for i in vertices(mg)
            neis = convert(Float64, length(neighbors(mg, i)))
            len = pyrs_df[pyrs_df[!, :ids] .== get_prop(mg, i, :id), :length][1]
            pyrs_df[pyrs_df[!, :ids] .== get_prop(mg, i, :id), :degree] = neis
            node_len[neis] = push!(get(node_len, neis, Float64[]), len)
        end
        
        nd = sort(collect(keys(node_len)))
        mlen = map(x -> mean(node_len[x]), nd)
        p = plot(nd[nd .< 30], mlen[nd .< 30], markershape=:auto, line=false, legend=false, markersize=6, xtickfont=Plots.font(13), xticks=[0,5,10,15,20,25,30], xlim=[0,30], ytickfont=Plots.font(13))
        p = plot!(nd[nd .< 30], 3625 .+ 3700*nd[nd .< 30])
        display(p)
        #savefig("./pics/length_nd_SD.pdf")
        p = plot(nd, mlen, markershape=:auto, line=false)
        p = plot!(nd[nd .< 30], 3625 .+ 3700*nd[nd .< 30])
        display(p)
        q = plot(log10.(pyrs_df[!, :degree]), log10.(pyrs_df[!, :length]), markershape=:auto, line=false, markersize = 1.0)
        gui(q)
    end
                                                            
    function Pyr_nd_bins_length(pyrs_df::DataFrame, fr_starts::Array, fr_ends::Array)
        low_freq_len = pyrs_df[fr_starts[1] .<= pyrs_df[:degree] .<= fr_ends[1] , :length]
        mid_freq_len = pyrs_df[fr_starts[2] .<= pyrs_df[:degree] .<= fr_ends[2] , :length]
        high_freq_len = pyrs_df[fr_starts[3] .<= pyrs_df[:degree] .<= fr_ends[3] , :length]
        #p = histogram(low_freq_len, bins=:scott, weights=repeat(1:5, outer=200))
    end                                                        

                                                            
    function Modify_DF_for_RandForset!(a::Interval_object, b::Interval_object, segs_int::Interval_object, graph_obj::Graph_pars_object, mean_covs::Array{Float64,1})
        
        df = DataFrame(length=Float64[], self=Int64[], double=Float64[], degree=Float64[], comp_s=Float64[], mean_cov=Float64[], breaks = Int64[], max_cov = Int64[], intra_frac = Float64[], max_len = Int64[])

        d_e = Dict()
        breaks = Dict()
        ls = Int64[]
        b.int_DF[!, :mean_cov] = mean_covs
                                                                
        b.int_DF[!, :comp_size] = zeros(size(b.int_DF)[1])
        ccs = connected_components(graph_obj.mg)
        for com in ccs
            len = length(com)
            map(x -> get_prop(graph_obj.mg, x, :id) |> x -> b.int_DF[b.int_DF[!, :ids] .== x, :comp_size] = len, com)
        end

        for i in keys(graph_obj.edges_double)
            d_e[i[1]] = get(d_e, i[1], 0) + graph_obj.edges_double[i]
            d_e[i[2]] = get(d_e, i[2], 0) + graph_obj.edges_double[i]
        end
    
        for i in keys(b.int_tree)
            for j in collect(intersect(b.int_tree[i], a.int_tree[i]))
                breaks[j[1].value] = get(breaks, j[1].value, reshape(Int64[], (0, 2)))
                breaks[j[1].value] = vcat(breaks[j[1].value], [j[2].first j[2].last])
            end
        end

        for i in keys(breaks)
            l = 0
            for j in 0:10
                l += length(unique(breaks[i][map(x -> x + j - 5 in breaks[i][:,2], breaks[i][:,1]),2]))
            end
            push!(ls, l)
        end

        b.int_DF[!, :pyr_breaks] = ls

        for i in eachrow(b.int_DF)
            deg = i[:degree]
            len = i[:length]
            mcov = i[:mean_cov]
            brks = i[:pyr_breaks]
            comp_s = i[:comp_size]
            self = get(graph_obj.self_loops, i[:ids], 0)
            double = get(d_e, i[:ids], 0)
            sds = collect(intersect(a.int_tree[i[:chr]], (i[:coor_s], i[:coor_e])))
            max_len = maximum(map(x -> x.last - x.first + 1, sds))
            intra_frac = mean(map(x -> i[:chr] == a.int_DF[a.int_DF[!, :ids] .== x.value, :chr2][1] ? 1.0 : 0.0, sds))
            segs = collect(intersect(segs_int.int_tree[i[:chr]], (i[:coor_s], i[:coor_e])))
            max_cov = maximum(map(x -> x.value, segs))
            push!(df, [len self double deg comp_s mcov brks max_cov intra_frac max_len])
        end
                                                                
        b.int_DF[!, :intra_frac] = df[:, :intra_frac]

        lms = lm(@formula(length ~ self + double + degree + comp_s + mean_cov + breaks + max_cov + intra_frac), df)
        display(r2(lms))
        display(lms)
        return df
    end
                                                            
end                                                       
                                                            
##############################################################################################
                                                            
module SegmentBlocks                                                            

    using IntervalTrees, Main.SDmoduleInit, DataFrames, StatsBase, Statistics, GenomicFeatures, IntervalSets, Histograms, Plots, ArndtLabTheme
    pyplot()
    theme(:arndtlab)

    export pyramid_coverage_calc, plot_segments_distribution, itree_from_df
                                                        
    function find_segments_from_it(int_SDs::Interval_object, int_pyrs::Interval_object)
        dic_starts, dic_ends = Dict(), Dict()
        for chr in keys(int_pyrs.int_tree)
            for it in int_pyrs.int_tree[chr]
                sds_set = intersect(int_SDs.int_tree[chr], it)
                for sd in sds_set                                                        
                    dic_starts[chr] = get(dic_starts, chr, Int64[])
                    dic_ends[chr] = get(dic_ends, chr, Int64[])
                    push!(dic_starts[chr], sd.first)
                    push!(dic_ends[chr], sd.last)
                end                                       
            end                                                        
        end
        out_df = find_segments_second(dic_starts, dic_ends)
        itree = itree_from_df(out_df)
        return Interval_object(out_df, itree)
    end
    
    
    function find_segments_second(dic_starts::Dict, dic_ends::Dict)
        out_df = DataFrame(chr = String[], coor_s = Int64[], coor_e = Int64[], length = Int64[], cover = Int64[])
        map(x -> sort!(dic_starts[x]), collect(keys(dic_starts)))
        map(x -> sort!(dic_ends[x]), collect(keys(dic_ends)))
        for chr in keys(dic_starts)
            st_ind, en_ind, cover = 1, 1, 0
            covers = Int64[]                                        
            push!(dic_starts[chr], 5e10)
            push!(dic_ends[chr], 5e10)
            while en_ind <= length(dic_ends[chr]) - 1
                if dic_starts[chr][st_ind] < dic_ends[chr][en_ind]
                    cover += 1
                    st_ind += 1
                elseif dic_starts[chr][st_ind] >= dic_ends[chr][en_ind]
                    cover -= 1
                    en_ind += 1                                                
                end
                push!(covers, cover)
            end
            
            coors = sort(append!(dic_starts[chr], dic_ends[chr]))
            for i in 1:(length(coors) - 3)                          # because added 2 big values;
                push!(out_df, (chr, coors[i], coors[i + 1], coors[i + 1] - coors[i], covers[i]))            
            end
        end
        out_df = out_df[(out_df[!, :length] .!= 0) .& (out_df[!, :cover] .!= 0), :]
        return out_df                                                        
    end
                                                            
                                                            
    function itree_from_df(cnvs_df::DataFrame, comp::String = "cover")
        Itree::Dict{String, IntervalTrees.IntervalBTree} = Dict()

        for line in eachrow(cnvs_df)
            Itree[line[:chr]] = get(Itree, line[:chr], IntervalTree{Int, IntervalValue{Int, Int}}())
            if comp == "cover"
                push!(Itree[line[:chr]], IntervalValue(line[:coor_s], line[:coor_e], line[:cover]))
            elseif comp == "nothing"
                push!(Itree[line[:chr]], IntervalValue(line[:coor_s], line[:coor_e], 0))
            end
        end
        
        return Itree
    end
                                                         
    
    function pyramid_coverage_calc(int_pyrs::Interval_object, int_sds::Interval_object)
        mean_covs = Float64[]
        all_sds::Array{GenomicFeatures.Interval, 1} = []
        
        for line in eachrow(int_pyrs.int_DF)
            len = line[:length]
            chr, coor1, coor2 = line[:chr], line[:coor_s], line[:coor_e]
                        
            sds = collect(intersect(int_sds.int_tree[chr], (coor1, coor2)))
            sds2 = [IntervalSets.intersect(coor1..coor2, i.first..i.last) for i in sds]
            sds3 = [(i.left - coor1)/len..(i.right - coor1)/len for i in sds2]
            sds4 = [GenomicFeatures.Interval("1", 1 + floor(Int, 10000*i.left), 1 + floor(Int, 10000*i.right)) for i in sds3]
            
            if length(sds4) > 0
                append!(all_sds, sds4)
                cov_loc = coverage(sort(sds4))
                mean_cov = sum(map(x -> (rightposition(x) - leftposition(x) + 1)*metadata(x), cov_loc))/10001
                push!(mean_covs, mean_cov)
            else
                push!(mean_covs, 0.0)
            end
        end
        
        covers = coverage(sort(all_sds))
        return mean_covs, covers
    end
          

    function plot_segments_distribution(a::Interval_object, b::Interval_object; len_thr::Int64=1000)

        segs_int = SegmentBlocks.find_segments_from_it(a, b)

        sf = segs_int.int_DF[segs_int.int_DF[!, :length] .> len_thr, :]

        cc_h = Histograms.Histogram(collect(sf[!, :cover]), N=40, scale=:log10)
        cc_p = plot(cc_h, xscale=:log10, yscale=:log10, markershape=:auto, line=false)
        x = range(1, 100, length=10); y = 2*x.^(-2);
        cc_p = plot!(x, y)
        gui(cc_p)

        display("Only long segments genome fraction:")
        display(sum(segs_int.int_DF[segs_int.int_DF[!, :length] .> len_thr, :length])/2.88e9)
        display("All segments genome fraction:")
        display(sum(segs_int.int_DF[!, :length])/2.88e9)
        
        return segs_int
    end
                                                        
end                                       
                                                        
###############################################################################################

module ApeLineage
    
    using StatsBase, DataFrames, Plots, IntervalTrees, ArndtLabTheme, Statistics
    pyplot()
    theme(:arndtlab)
                                                            
    export Read_ape_SDs_2_DF, Itree_from_DF_ape_lineage, Extract_from_ape_lineage
                                                                    
    function Read_ape_SDs_2_DF(filename::String, col::Int)
        df = DataFrame(chr=String[], coor_s=Int64[], coor_e=Int64[], lineage=String[], mean=Float64[], val=Float64[])
        open(filename) do sdfile
            lines = readlines(sdfile);
            coors_str = map(x -> strip(x) |> split |> x -> x[[1,2,3,4,end,end-col]], lines)
            map(x -> push!(df ,[replace(x[1], r"chr" => "") [parse.(Int64, x[2:3])...]' x[4] [parse.(Float64, x[5:6])...]']), coors_str)
        end
        return df
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

    function Extract_from_ape_lineage(SD_df::DataFrame, ape_lineage::Main.SDmoduleInit.Interval_object, comp::String)
        out_arr = Float64[]
        for (i, line) in enumerate(eachrow(SD_df))
            chr, coor_s, coor_e = line[:chr], line[:coor_s], line[:coor_e]
            over = collect(intersect(ape_lineage.int_tree[chr], (coor_s, coor_e)))
            if length(over) > 0
                append!(out_arr, [i.value for i in over])
            end
        end
        out_arr = sort(out_arr)
        if length(out_arr) == 0
            return []
        end
        if comp == "all"
            return out_arr
        elseif comp == "restrict"
            #perc = Int(round(length(out_arr)*0.05)) + 1
            q25, q75 = quantile!(out_arr, 0.25), quantile!(out_arr, 0.75)
            iq = q75 - q25
            return out_arr[q25 - 3*iq .< out_arr .< q75 + 3*iq]     # out_arr[1:end-perc]
        end
    end
end

###############################################################################################

module CNVspart                                                            

    using DataFrames, Main.SDmoduleInit, IntervalTrees, IntervalSets, StatsBase, GenomicFeatures, Plots, ArndtLabTheme
    pyplot()
    theme(:arndtlab)
    
    export read_CNVs_to_DF, shuffle_DF_to_Intobject, overlap_CNVs_pyrs!, CNV_rare_fixed_overlap, find_gene_overlap!, itree_from_df
                                                        
    function read_CNVs_to_DF(file_SDs::String = "../CNV_pop_cancer/CNVs_only_popul_hg38.txt", ac_thr::Array{Int64,1} = [1, 2504], samplesize::Int64=5008)
        
        cnvs_df = DataFrame(chr = String[], coor_s = Int64[], coor_e = Int64[], ac = Int64[])

        open(file_SDs) do sdfile
            lines = readlines(sdfile)
            acs = []
            for l in lines
                a = split(l)
                cn_vals = [x.match for x in eachmatch(r"<CN[0-9]+>", a[5])]
                m = match(r"AC=([0-9,]+);", a[8])
                ac_vals = parse.(Int64, split(m.captures[1], ","))
                ac = sum(ac_vals[cn_vals .!= "<CN0>"])
                push!(acs, ac <= samplesize/2 ? ac : samplesize - ac)
            end
                                                                        
            #acs = map(x -> strip(x) |> x -> split(x)[8] |> x -> split(x, r",|;")[1] |> x -> split(x, "=")[2] |> x -> parse(Int64, x) |> x -> x <= samplesize/2 ? x : samplesize - x, lines)
            chrs = map(x -> split(x)[1], lines)
            
            coor_s = map(x -> split(x)[2] |> x -> parse(Int64, x), lines)
            coor_e = map(x -> split(x)[8] |> x -> match(r"END=[0-9]+", x).match |> x -> split(x, "=")[2] |> x -> parse(Int64, x), lines)
            chrs = map(x -> x in ["X", "Y", "M"] ? "0" : x, chrs)                   # maybe change to "M"
            [push!(cnvs_df, [chrs[i] coor_s[i] coor_e[i] acs[i]]) for i in 1:length(chrs)];
            cnvs_df = cnvs_df[(cnvs_df[!, :chr] .!== "0") .& (ac_thr[1] .<= cnvs_df[!, :ac] .<= ac_thr[2]), :]
        end
        
        Itree = itree_from_df(cnvs_df)
                                                            
        return Interval_object(cnvs_df, Itree)                                      
    end
                                                        
    function overlap_CNVs_pyrs!(cnv_int::Interval_object, pyrs_int::Interval_object, comp::String = "real"; pnum::Int64 = 1000)
        
        if comp == "real" 
            pyrs_int.int_DF[!, :cnv_over] = zeros(Int64, size(pyrs_int.int_DF)[1])
        elseif comp == "shuffle"
            pyrs_int.int_DF[!, :shuff_over] = zeros(Int64, size(pyrs_int.int_DF)[1])
        elseif comp == "perm_over"
            pyrs_int.int_DF[!, :cnv_over] = zeros(Int64, size(pyrs_int.int_DF)[1])
            pyrs_int.int_DF[!, :over_len] = zeros(Float64, size(pyrs_int.int_DF)[1])
            pyrs_int.int_DF[!, :nor_left] = zeros(Float64, size(pyrs_int.int_DF)[1])
            pyrs_int.int_DF[!, :nor_right] = zeros(Float64, size(pyrs_int.int_DF)[1])
            perm_overs = reshape(Float64[], (0, pnum))
        elseif comp == "shuffle2"
            pyrs_int.int_DF[!, :shuff_over] = zeros(Int64, size(pyrs_int.int_DF)[1])
            pyrs_int.int_DF[!, :nor_left] = zeros(Float64, size(pyrs_int.int_DF)[1])
            pyrs_int.int_DF[!, :nor_right] = zeros(Float64, size(pyrs_int.int_DF)[1])
        end
        
        for k in keys(pyrs_int.int_tree)
            if k in keys(cnv_int.int_tree)
                for inter in intersect(pyrs_int.int_tree[k], cnv_int.int_tree[k])
                    if comp == "real"
                        pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :cnv_over] = pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :cnv_over][1] + 1
                    elseif comp == "shuffle"
                        pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :shuff_over] = pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :shuff_over][1] + 1
                    elseif comp == "perm_over"
                        l_over, lens_exp, nor_left, nor_right = calc_overlap_fraction_predict(inter[1].first..inter[1].last, inter[2].first..inter[2].last, pnum)
                        if pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :cnv_over][1] == 0
                            pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :nor_left] = nor_left
                            pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :nor_right] = nor_right
                        else 
                            pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :nor_left] = 0.0
                            pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :nor_right] = 1.0
                        end
                        pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :cnv_over] = pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :cnv_over][1] + 1
                        pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :over_len] = pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :over_len][1] + l_over
                        perm_overs = vcat(perm_overs, lens_exp')
                    elseif comp == "shuffle2"
                        l_over, lens_exp, nor_left, nor_right = calc_overlap_fraction_predict(inter[1].first..inter[1].last, inter[2].first..inter[2].last, 1)
                        if pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :shuff_over][1] == 0
                            pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :nor_left] = nor_left
                            pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :nor_right] = nor_right
                        else 
                            pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :nor_left] = 0.0
                            pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :nor_right] = 1.0
                        end
                        pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :shuff_over] = pyrs_int.int_DF[pyrs_int.int_DF[!, :ids] .== inter[1].value, :shuff_over][1] + 1
                    end
                end
            end
        end

        if comp == "perm_over"
            sums = [sum(perm_overs[!, i]) for i in 1:size(perm_overs)[2]]
            pval = countmap(sums .> sum(pyrs_int.int_DF[!, :over_len]))[:true]/pnum
            display("p-value = $(pval) it means nothing")
        end
        
        if comp == "shuffle2" || comp == "perm_over"
            only_over = pyrs_int.int_DF[pyrs_int.int_DF[!, :nor_right] - pyrs_int.int_DF[!, :nor_left] .> 0.0, :]
            ints = [GenomicFeatures.Interval("1", 1 + floor(Int, 10000*i[:, :nor_left]), 1 + floor(Int, 10000*i[:, :nor_right])) for i in eachrow(only_over)]
            covers = coverage(sort(ints))
            return covers, sort(ints)
        end
    end                                                                                                      
    
    function shuffle_DF_to_Intobject(in_df::DataFrame, chr_file::String = "chrom.size.hg38")
        ch_size = Dict()
        df_sim = DataFrame(chr = String[], coor_s = Int64[], coor_e = Int64[])

        open(chr_file) do chr_f
            lines = readlines(chr_f);
            map(x -> strip(x) |> split |> x -> ch_size[strip(x[1], ['c', 'h', 'r'])] = parse(Int64, x[2]), lines)
        end

        for line in eachrow(in_df)
            lensd = 1 + Int(line[:coor_e] - line[:coor_s])
            stpos = rand(1:(ch_size[line[:chr]] - lensd), 1)[1]
            enpos = stpos + lensd
            chr = line[:chr]
            push!(df_sim, [chr stpos enpos])                                                    
        end
        
        Itree = itree_from_df(df_sim, "no_id")
        
        return Interval_object(df_sim, Itree)                                                   
    end
    
    
    function itree_from_df(cnvs_df::DataFrame, comp::String = "ac_id")
        Itree::Dict{String, IntervalTrees.IntervalBTree} = Dict()

        for line in eachrow(cnvs_df)
            if comp == "ac_id"
                Itree[line[:chr]] = get(Itree, line[:chr], IntervalTree{Int, IntervalValue{Int, Int}}())
                push!(Itree[line[:chr]], IntervalValue(line[:coor_s], line[:coor_e], line[:ac]))
            elseif comp == "no_id"
                Itree[line[:chr]] = get(Itree, line[:chr], IntervalTree{Int, IntervalValue{Int, Int}}())
                push!(Itree[line[:chr]], IntervalValue(line[:coor_s], line[:coor_e], 0))
            else comp == "id_id"
                Itree[line[:chr]] = get(Itree, line[:chr], IntervalTree{Int, IntervalValue{Int, String}}())
                push!(Itree[line[:chr]], IntervalValue(line[:coor_s], line[:coor_e], line[:ids]))
            end
        end
        return Itree
    end
    
    
    function calc_overlap_fraction_predict(pyr::IntervalSets.Interval, cnv::IntervalSets.Interval, pnum::Int64)
        l1, l2 = duration(pyr), duration(cnv)
        over = IntervalSets.intersect(pyr, cnv)
        nor_left = (over.left - pyr.left)/l1
        nor_right = 1.0 - (pyr.right - over.right)/l1
                                                                                
        l_over = duration(over)
        lens_exp::Array{Float64,1} = []
        
        for i in 1:pnum
            p1 = rand(1:(l1+l2))
            p2 = p1 + l1
            l_perm = float(duration(IntervalSets.intersect(p1..p2, l1..(l1+l2))))
            push!(lens_exp, l_perm)
        end
        
        return float(l_over), lens_exp, nor_left, nor_right
    end
                                                                                                
    function plot_coverage(covers, num::Int64 = 1)
        xs, ys = [], []
        for i in covers
            append!(xs, [leftposition(i), rightposition(i)])
            append!(ys, [metadata(i), metadata(i)])
        end
        p = plot(xs, ys/num)
        return p
    end
                                                                                                
    function CNV_rare_fixed_overlap(b::Interval_object, CNV_int::Interval_object, comp::String = "norm", comp2::String = "both", n::Int64 = 100)
        
        if comp2 == "both"
            piece = deepcopy(b)
        elseif comp2 == "gene"
            piece_df = deepcopy(b.int_DF[b.int_DF[!, :gene_over] .== 1.0, :])
            piece_tree = itree_from_df(piece_df, "id_id")
            piece = Interval_object(piece_df, piece_tree)
        elseif comp2 == "nogene"
            piece_df = deepcopy(b.int_DF[b.int_DF[!, :gene_over] .== 0.0, :])
            piece_tree = itree_from_df(piece_df, "id_id")
            piece = Interval_object(piece_df, piece_tree)
        end
        
        mer = zeros(size(piece.int_DF)[1])
        bins = Int64[1 1; 2 5; 6 200]
        overs_df = DataFrame(low_f = Float64[], middle_f = Float64[], high_f = Float64[])
        
        if comp == "shuff"
            for i in 1:n
                cnv_shuff = CNVspart.shuffle_DF_to_Intobject(CNV_int.int_DF)
                CNVspart.overlap_CNVs_pyrs!(cnv_shuff, piece, "shuffle")
                mer .+= piece.int_DF[!, :shuff_over]
                push!(overs_df, CNVspart.binning_lin_frequency(piece, bins, :shuff_over))
            end
            mer = mer/n
           # b.int_DF[!, :shuff_over] = mer
        elseif comp == "norm"
            push!(overs_df, CNVspart.binning_lin_frequency(piece, bins, :cnv_over))
        end
        
        return overs_df
    end
                                                                                                
    function binning_lin_frequency(b::Interval_object, bins::Array{Int64,2}, col)
        c_mer = Float64[]
        for i in eachrow(bins)
            spl = b.int_DF[i[1] .<= b.int_DF[!, :degree] .<= i[2], :]
            push!(c_mer, sum(spl[!, col])/size(spl)[1])
        end
        return c_mer
    end
                                                                                                    
    function find_gene_overlap!(any_int::Interval_object, filename::String)
        gene_int = read_gene_file(filename)
        any_int.int_DF[!, :gene_over] = zeros(size(any_int.int_DF)[1])
        for row in eachrow(any_int.int_DF)
            chr, coor_s, coor_e = row[:chr], row[:coor_s], row[:coor_e]
            over = collect(intersect(gene_int.int_tree[chr], (coor_s, coor_e)))
            if length(over) > 0
                row[:gene_over] = 1.0
            end
        end
    end
                                                                                                        
    function read_gene_file(filename::String)
        genes_df = DataFrame(chr = String[], coor_s = Int64[], coor_e = Int64[])
        open(filename) do genefile
            lines = readlines(genefile)
            chrs = map(x -> split(x)[1], lines)
            coors = map(x -> split(x)[2:3] |> x -> parse.(Int64, x), lines)
            chrs = map(x -> replace(x, r"chr" => "") |> x -> x in ["X", "Y"] ? "0" : x, chrs)    # X or Y are excluded here;
            [push!(genes_df, [chrs[i] coors[i][1] coors[i][2]]) for i in 1:length(chrs)];
            genes_df = genes_df[genes_df[!, :chr] .!== "0",:]
        end
        Itree = itree_from_df(genes_df, "no_id")
        
        return Interval_object(genes_df, Itree)
    end
    
end
                                                                    
###############################################################################################          
                                                                                                
module MST_from_SD

    using DataFrames, Main.SDmoduleInit, LightGraphs, MetaGraphs, SparseArrays, GraphPlot, Distributions, Histograms, Statistics, StatsBase, Plots
    
    export graph_nosec_from_graph, assign_weights_edges!, wrapper_MST, filt_unneeded, test_sharpness_distribution


    function modify_adj_matrix(graph::MetaGraph)
        Ad_pr = float(deepcopy(adjacency_matrix(graph)))
        for i in vertices(graph)
            neis = neighbors(graph, i)
            for n in neis
                neis2 = neighbors(graph, n)
                mat_inds = findall(x -> x in neis, neis2)
                wei = (length(mat_inds) + 1)/length(neis)
                Ad_pr[i, n] = wei
            end
        end
        Ad_pr = sparse(Matrix(Ad_pr) ./ sum(Ad_pr, dims=2))
        return Ad_pr
    end


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
                if LightGraphs.weights(graph)[i, n] == 1.0
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
                        wei = minimum([1.001 - (length(mat_inds) + 1)/length(neis), LightGraphs.weights(graph)[i, n]])
                    elseif comp == "normal_noise"
                        wei = minimum([1.001 - (length(mat_inds) + 1)/length(neis), LightGraphs.weights(graph)[i, n]]) + rand(Uniform(-0.00005, 0.00005))
                    elseif comp == "sharp_out"
                        if (Edge(i, n) in suspicious_edges) | (Edge(n, i) in suspicious_edges)
                            wei = 100.0
                        else
                            wei = minimum([1.001 - (length(mat_inds) + 1)/length(neis), LightGraphs.weights(graph)[i, n]]) + rand(Uniform(-0.00005, 0.00005))
                        end
                    elseif comp == "betweenness_weight"
                        wei = minimum([1.001 - bws[n]/sum_bws, LightGraphs.weights(graph)[i, n]])
                    elseif comp == "shuffle"
                        wei = minimum([1.001 + rand(), LightGraphs.weights(graph)[i, n]])
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


    function samarkand(arr_in)
        arr_out = Int64[]
        for i in arr_in
            append!(arr_out, i.src)
            append!(arr_out, i.dst)
        end
        sort!(arr_out)
        return counts(arr_out)
    end


    function wrapper_MST(copy_gr, edg_copy_real; verbose=false, meanover=1)
        for i in 1:nv(copy_gr)
            id = "id$i"
            set_prop!(copy_gr, i, :id, id)
        end
        copy_gr_copy = deepcopy(copy_gr)
        assign_weights_edges!(copy_gr, comp="normal_noise")       # "shuffle", "normal_noise", betweenness_weight"
        edg_copy_mst = kruskal_mst(copy_gr)
        edg_copy_mst = map(x -> x.src < x.dst ? x : Edge(x.dst, x.src), edg_copy_mst)
    
        assign_weights_edges!(copy_gr_copy, comp="shuffle")
        edg_copy_mst_shuff = kruskal_mst(copy_gr_copy)
        edg_copy_mst_shuff = map(x -> x.src < x.dst ? x : Edge(x.dst, x.src), edg_copy_mst_shuff)

        brr = collect(map(x -> x in edg_copy_mst ? 1 : 0, edg_copy_real))
        frac_Kruskal = mean(brr)
        brr = collect(map(x -> x in edg_copy_mst_shuff ? 1 : 0, edg_copy_real))
        frac_random_MST = mean(brr)
        if verbose
            println("Fraction of edges match (MST Kruskal):", frac_Kruskal)
            println("Fraction of edges match (Random MST):", frac_random_MST)
        end
        j, n_deg_mst = 1, []
        while j < meanover
            assign_weights_edges!(copy_gr, comp="normal_noise")
            append!(edg_copy_mst, kruskal_mst(copy_gr))
            j += 1
        end
        n_deg_mst = samarkand(edg_copy_mst)
        n_deg_mst_shuff = samarkand(edg_copy_mst_shuff)
        n_deg_simple = samarkand(edges(copy_gr))
        n_deg_real = samarkand(edg_copy_real)
        n_deg_pagerank = pagerank(copy_gr)   # or other centralities from https://docs.juliahub.com/LightGraphs/Xm08G/1.3.5/centrality/
    
        if verbose
            println("Variance explained for MST Kruskal:", cor(n_deg_real, n_deg_mst)^2);
            println("Variance explained for Random MST:", cor(n_deg_real, n_deg_mst_shuff)^2);
            println("Variance explained for Simple:", cor(n_deg_real, n_deg_simple)^2);
            println("Variance explained for PageRank:", cor(n_deg_real, n_deg_pagerank)^2);
        end
        n_deg_mst = n_deg_mst .+ rand(length(n_deg_mst)) .- 0.5
        n_deg_real = n_deg_real .+ rand(length(n_deg_real)) .- 0.5
        return transpose([frac_Kruskal, frac_random_MST, cor(n_deg_real, n_deg_mst)^2, cor(n_deg_real, n_deg_mst_shuff)^2, cor(n_deg_real, n_deg_simple)^2, cor(n_deg_real, n_deg_pagerank)^2])
        #p = plot(n_deg_real, n_deg_mst, line=false, markershape=:auto, markersize=1.0, aspect_ratio=:equal);
        #gui(p)
    end


    function filt_unneeded(edgs, graph; nodes_out=nothing)
        ccs = connected_components(graph)
        lens = map(length, ccs)
        bcs = ccs[lens .== maximum(lens)][1]
        bcs = filter(x -> length(neighbors(graph, x)) > 1, bcs)
        edgs = filter(x -> (x.src in bcs) & (x.dst in bcs), edgs)
        if !(nodes_out == nothing)
            edgs = filter(x -> !((get_prop(graph, x.src, :id) in nodes_out) | (get_prop(graph, x.dst, :id) in nodes_out)), edgs)
        end
        return edgs
    end


    function test_sharpness_distribution(gr1, gr2, edg_lis, N, mode; nodes_out=nothing, edges_ident=nothing)
        dif = collect(edges(difference(gr1, gr2)))
        dif_bc = MST_from_SD.filt_unneeded(dif, gr1; nodes_out)
        edg_lis_bc = MST_from_SD.filt_unneeded(edg_lis, gr1; nodes_out)
        edg_all_bc = MST_from_SD.filt_unneeded(collect(edges(gr1)), gr1; nodes_out)
        if mode == "sharp edges"
            match = map(x -> x in dif_bc ? 1 : 0, edg_lis_bc)
            val = sum(match)/length(edg_lis_bc)
        elseif mode == "identity"
            if !(nodes_out == nothing)
                edg_lis_f = filter(x -> !((get_prop(gr1, x.src, :id) in nodes_out) | (get_prop(gr1, x.dst, :id) in nodes_out)), edg_lis)
                vals = [edges_ident[(get_prop(gr1, i.src, :id), get_prop(gr1, i.dst, :id))] for i in edg_lis_f]
                all_edges_f = filter(x -> !((get_prop(gr1, x.src, :id) in nodes_out) | (get_prop(gr1, x.dst, :id) in nodes_out)), collect(edges(gr1)))
            else
                vals = [edges_ident[(get_prop(gr1, i.src, :id), get_prop(gr1, i.dst, :id))] for i in edg_lis]
            end
            val = mean(vals .> 0.99)
        end
        display(val)
        out_lis = [] 
        for i in 1:N
            if mode == "sharp edges"
                shuff_dif = sample(edg_all_bc, length(edg_lis_bc), replace=false)
                match = map(x -> x in dif_bc ? 1 : 0, shuff_dif)
                append!(out_lis, sum(match)/length(shuff_dif))
            elseif mode == "identity"
                if !(nodes_out == nothing)
                    shuff_dif = sample(all_edges_f, length(edg_lis_f), replace=false)
                else
                    shuff_dif = sample(collect(edges(gr1)), length(edg_lis), replace=false)
                end
                vals = [edges_ident[(get_prop(gr1, i.src, :id), get_prop(gr1, i.dst, :id))] for i in shuff_dif]
                append!(out_lis, mean(vals .> 0.99))
            end
        end
        if mode == "sharp edges"
            display("percantage of random permutations with less sharp borders:")
            display(length(out_lis[out_lis .<= val])/length(out_lis))
        else
            display("percantage of random permutations with higher fraction of hihgh identity edges:")
            display(1.0 - length(out_lis[out_lis .<= val])/length(out_lis))
        end
        return val, out_lis
    end

end

##########################################################################################

module Annotate_SDs

    using DataFrames, Main.SDmoduleInit, LightGraphs, MetaGraphs, DataFramesGenomics, Printf, Statistics
    
    export Read_genes_and_gaps_2_DF, Read_repli_recomb_dnase_2_DF, Read_repeats_2_DF, construct_DF!, go_thru_repeat_overlaps!, go_thru_thru!, add_all_columns!, modify_dfs_gaps!, modify_dfs_means!, main_modification, finzalize_df_cols!, GC_flanks_measure!

    function Read_genes_and_gaps_2_DF(filename::String)
        df = DataFrame(chr=String[], coor_s=Int64[], coor_e=Int64[])
        open(filename) do sdfile
            lines = readlines(sdfile);
            coors_str = map(x -> strip(x) |> split |> x -> x[[1,2,3]], lines)
            map(x -> push!(df ,[replace(x[1], r"chr" => "") [parse.(Int64, x[2:3])...]']), coors_str)
        end
        return df
    end


    function Read_repli_recomb_dnase_2_DF(filename::String)
        df = DataFrame(chr=String[], coor_s=Int64[], coor_e=Int64[], value=Float64[])
        open(filename) do sdfile
            lines = readlines(sdfile);
            coors_str = map(x -> strip(x) |> split |> x -> x[[1,2,3,4]], lines)
            map(x -> push!(df ,[replace(x[1], r"chr" => "") [parse.(Int64, x[2:3])...]' parse(Float64, x[4])]), coors_str)
        end
        return df
    end


    function Read_repeats_2_DF(filename::String)
        df = DataFrame(chr=String[], coor_s=Int64[], coor_e=Int64[], rep_class=String[], rep_family=String[])
        open(filename) do sdfile
            lines = readlines(sdfile);
            coors_str = map(x -> strip(x) |> split |> x -> x[[1,2,3,4,5]], lines)
            map(x -> push!(df ,[replace(x[1], r"chr" => "") [parse.(Int64, x[2:3])...]' x[4] x[5]]), coors_str)
        end
        return df
    end


    function construct_DF!(df::DataFrame, segments_df::DataFrame, column::Symbol)
        # for telo/centro distance
        for line in eachrow(df)
            coors = Int[]
            chr, coor_s, coor_e = line[:chr], line[:coor_s], line[:coor_e]
            coors = vcat(segments_df[segments_df[!, :chr] .== chr, :coor_s], segments_df[segments_df[!, :chr] .== chr, :coor_e])
            vals = vcat(coor_s .- coors, coors .- coor_e)
            line[column] = minimum(vals[vals .>= 0])
        end
    end


    function construct_DF!(pyrs_obj::Main.SDmoduleInit.Interval_object, sds_obj::Main.SDmoduleInit.Interval_object, column::Symbol)
        # for intrafraction
        for i in eachrow(pyrs_obj.int_DF)
            sds = collect(intersect(sds_obj.int_tree[i[:chr]], (i[:coor_s], i[:coor_e])))
            intra_frac = mean(map(x -> i[:chr] == sds_obj.int_DF[sds_obj.int_DF[!, :ids] .== x.value, :chr2][1] ? 1.0 : 0.0, sds))
            i[column] = round(intra_frac, digits=3)
        end
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


    function construct_DF!(df::DataFrame, border::Int, gap::Int)
        df[!, :used_coor_l_s], df[!, :used_coor_l_e] = zeros(Int, size(df)[1]), zeros(Int, size(df)[1])
        df[!, :used_coor_r_s], df[!, :used_coor_r_e] = zeros(Int, size(df)[1]), zeros(Int, size(df)[1])
        for line in eachrow(df)
            chr, coor_s, coor_e = line[:chr], line[:coor_s], line[:coor_e]
            if coor_s - gap <= coor_e + gap
                line[:used_coor_l_s], line[:used_coor_l_e] = coor_s - gap - border, coor_s - gap
                line[:used_coor_r_s], line[:used_coor_r_e] = coor_e + gap, coor_e + gap + border
            else
                line[:used_coor_l_s], line[:used_coor_l_e] = coor_s - border, coor_s
                line[:used_coor_r_s], line[:used_coor_r_e] = coor_e, coor_e + border
            end
        end
    end


    function my_chr2int(vec)
        max_val = maximum(skipmissing(chr2int.(vec)))
        chr2int_dic, int2chr_dic = Base.ImmutableDict("" => 0), Base.ImmutableDict(0 => "")
        for i in vec
            if ismissing(chr2int(i))
                max_val += 1
                chr2int_dic = Base.ImmutableDict(chr2int_dic, i => max_val)
                int2chr_dic = Base.ImmutableDict(int2chr_dic, max_val => i)
            else
                chr2int_dic = Base.ImmutableDict(chr2int_dic, i => chr2int(i))
                int2chr_dic = Base.ImmutableDict(int2chr_dic, chr2int(i) => i)
            end
        end
        return chr2int_dic, int2chr_dic
    end


    function my_int2chr(vec, int2chr_dic)
        return [int2chr_dic[x] for x in vec]
    end


    function go_thru_repeat_overlaps!(df, repeats_df, col1::Symbol, col2::Symbol, comp::String; first2=:coor_s, last2=:coor_e)
        vec_chroms = unique(vcat(df.chr, repeats_df.chr))
        chr2int_dic, int2chr_dic = my_chr2int(vec_chroms)
        df[!, :chr] = [chr2int(i, chr2int_dic) for i in df.chr]
        repeats_df[!, :chr] = [chr2int(i, chr2int_dic) for i in repeats_df.chr]
        sort!(df, [:chr, col1, col2])
        over = eachoverlap(df, repeats_df, first1=col1, last1=col2, first2=first2, last2=last2)
        nams = ["DNA", "LINE", "LTR", "SINE", "Low_complexity", "Retroposon", "Satellite", "Simple_repeat",
                    "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", "RC"]
        nams2 = ["L1", "L2", "MIR", "Alu", "Satellite"]
        # One may need to change repeat families and subfamilies namesa when working with other species;
        
        for (x, y) in over
            col_c, col_f = repeats_df[y, :rep_class], repeats_df[y, :rep_family]
            if col_c in nams
                if comp == "left"
                    df[x, Symbol(col_c*"_l")] += 1
                elseif comp == "right"
                    df[x, Symbol(col_c*"_r")] += 1
                end
                if col_f in nams2
                    if comp == "left"
                        df[x, Symbol(col_f*"_s_l")] += 1
                    elseif comp == "right"
                        df[x, Symbol(col_f*"_s_r")] += 1
                    end
                end
            end
        end
        df[!, :chr] = my_int2chr(df.chr, int2chr_dic)
        repeats_df[!, :chr] = my_int2chr(repeats_df.chr, int2chr_dic)
    end


    function go_thru_thru!(df_both, repeats_df, border, gap_value)
        t1 = time()
        df_both[!, :DNA_l], df_both[!, :LINE_l], df_both[!, :LTR_l], df_both[!, :SINE_l] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        df_both[!, :Low_complexity_l], df_both[!, :Retroposon_l], df_both[!, :Satellite_l], df_both[!, :Simple_repeat_l] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        df_both[!, :rRNA_l], df_both[!, :snRNA_l], df_both[!, :scRNA_l], df_both[!, :srpRNA_l] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        df_both[!, :tRNA_l], df_both[!, :RC_l] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]);

        df_both[!, :DNA_r], df_both[!, :LINE_r], df_both[!, :LTR_r], df_both[!, :SINE_r] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        df_both[!, :Low_complexity_r], df_both[!, :Retroposon_r], df_both[!, :Satellite_r], df_both[!, :Simple_repeat_r] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        df_both[!, :rRNA_r], df_both[!, :snRNA_r], df_both[!, :scRNA_r], df_both[!, :srpRNA_r] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        df_both[!, :tRNA_r], df_both[!, :RC_r] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]);

        df_both[!, :L1_s_l], df_both[!, :L2_s_l] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        df_both[!, :MIR_s_l], df_both[!, :Alu_s_l], df_both[!, :Satellite_s_l] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
    
        df_both[!, :L1_s_r], df_both[!, :L2_s_r] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        df_both[!, :MIR_s_r], df_both[!, :Alu_s_r], df_both[!, :Satellite_s_r] = zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1]), zeros(Int, size(df_both)[1])
        construct_DF!(df_both, border, gap_value)
        go_thru_repeat_overlaps!(df_both, repeats_df, :used_coor_r_s, :used_coor_r_e, "right");
        go_thru_repeat_overlaps!(df_both, repeats_df, :used_coor_l_s, :used_coor_l_e, "left");
        display(time() - t1)
    end


    function add_all_columns!(df_in)
        df_in[!, :length] = df_in[!, :coor_e] .- df_in[!, :coor_s]
        df_in[!, :centro], df_in[!, :telo], df_in[!, :gaps] = zeros(Int, size(df_in)[1]), zeros(Int, size(df_in)[1]), zeros(Int, size(df_in)[1])
        df_in[!, :genes], df_in[!, :intra_frac] = zeros(Int, size(df_in)[1]), zeros(Float64, size(df_in)[1])
        df_in[!, :cpgisl_in], df_in[!, :cpgisl_bor], df_in[!, :ctcf] = zeros(Int, size(df_in)[1]), zeros(Int, size(df_in)[1]), zeros(Int, size(df_in)[1])

        df_in[!, :repli_in], df_in[!, :repli_bor] = zeros(Float64, size(df_in)[1]), zeros(Float64, size(df_in)[1]);
        df_in[!, :repli_bor_deriv] = zeros(Float64, size(df_in)[1])
        df_in[!, :repli_deriv], df_in[!, :repli_vari] = zeros(Float64, size(df_in)[1]), zeros(Float64, size(df_in)[1]);
        df_in[!, :recomb_in], df_in[!, :recomb_bor] = zeros(Float64, size(df_in)[1]), zeros(Float64, size(df_in)[1]);
        df_in[!, :dnase_in], df_in[!, :dnase_bor] = zeros(Float64, size(df_in)[1]), zeros(Float64, size(df_in)[1]);
    end


    function modify_dfs_gaps!(df_in, border, gap_val, cpgisl_obj, dnase_obj, recomb_obj, repli_obj, gaps_obj, ctcf_obj, repli_deriv_obj)
        construct_DF!(df_in, gaps_obj, :gaps, border=border, gap=gap_val);
        construct_DF!(df_in, cpgisl_obj, :cpgisl_bor, border=border, gap=gap_val);
        construct_DF!(df_in, dnase_obj, :dnase_bor, border=border, gap=gap_val);
        construct_DF!(df_in, ctcf_obj, :ctcf, border=border, gap=gap_val);
        construct_DF!(df_in, recomb_obj, :recomb_bor, border=border, gap=gap_val);
        construct_DF!(df_in, repli_obj, :repli_bor, border=border, gap=gap_val);
        construct_DF!(df_in, repli_deriv_obj, :repli_bor_deriv, border=border, gap=gap_val);
        
        m_repb = round(mean(df_in[df_in[!, :repli_bor] .!= 0, :repli_bor]), digits=3)
        m_recb = round(mean(df_in[df_in[!, :recomb_bor] .!= 0, :recomb_bor]), digits=3)
        m_dnab = round(mean(df_in[df_in[!, :dnase_bor] .!= 0, :dnase_bor]), digits=3)
         
        df_in[df_in[!, :repli_bor] .== 0, :repli_bor] .= m_repb;
        df_in[df_in[!, :recomb_bor] .== 0, :recomb_bor] .= m_recb;
        df_in[df_in[!, :dnase_bor] .== 0, :dnase_bor] .= m_dnab;
    end


    function modify_dfs_means!(df_in, column, var)
        df_in[df_in[!, column] .== 0, column] .= var;
    end


    function main_modification(df_in, border, telo_df, centro_df, genes_obj, cpgisl_obj, dnase_obj, recomb_obj, repli_obj, gaps_obj, ctcf_obj, repli_deriv_obj)
        construct_DF!(df_in, telo_df, :telo);
        construct_DF!(df_in, centro_df, :centro);
        construct_DF!(df_in, genes_obj, :genes, border=-1);
        construct_DF!(df_in, cpgisl_obj, :cpgisl_in, border=-1);
        construct_DF!(df_in, dnase_obj, :dnase_in, border=-1);
        construct_DF!(df_in, recomb_obj, :recomb_in, border=-1);
        construct_DF!(df_in, repli_obj, :repli_in, border=-1);
        construct_DF!(df_in, repli_obj, :repli_deriv, border=-1, comp="derivative");
        construct_DF!(df_in, repli_obj, :repli_vari, border=-1, comp="variance");
        
        m_repi = round(mean(df_in[df_in[!, :repli_in] .!= 0, :repli_in]), digits=3)
        m_repder = round(mean(df_in[df_in[!, :repli_deriv] .!= 0, :repli_deriv]), digits=3)
        m_repvar = round(mean(df_in[df_in[!, :repli_vari] .!= 0, :repli_vari]), digits=3)
        m_reci = round(mean(df_in[df_in[!, :recomb_in] .!= 0, :recomb_in]), digits=3)
        m_dnai = round(mean(df_in[df_in[!, :dnase_in] .!= 0, :dnase_in]), digits=3)
        
        modify_dfs_means!(df_in, :repli_in, m_repi);
        modify_dfs_means!(df_in, :repli_deriv, m_repder);
        modify_dfs_means!(df_in, :repli_vari, m_repvar);
        modify_dfs_means!(df_in, :recomb_in, m_reci);
        modify_dfs_means!(df_in, :dnase_in, m_dnai);
    end


    function finzalize_df_cols!(df_array, b, border, gaps_arr, telo_df, centro_df, genes_obj, cpgisl_obj, dnase_obj, recomb_obj, repli_obj, gaps_obj, ctcf_obj, repli_deriv_obj)
        t1 = time()
        for (i, s) in enumerate(df_array)
            modify_dfs_gaps!(s, border, gaps_arr[i], cpgisl_obj, dnase_obj, recomb_obj, repli_obj, gaps_obj, ctcf_obj, repli_deriv_obj);
            s[!, :intra_frac] = deepcopy(b.int_DF[!, :intra_frac])
            main_modification(s, border, telo_df, centro_df, genes_obj, cpgisl_obj, dnase_obj, recomb_obj, repli_obj, gaps_obj, ctcf_obj, repli_deriv_obj)
        end
        display(time() - t1)
    end


    function GC_flanks_measure!(df, GC_dict, AT_dict)
        l_s, l_e = df[!, "used_coor_l_s"], df[!, "used_coor_l_e"]
        r_s, r_e = df[!, "used_coor_r_s"], df[!, "used_coor_r_e"]
        chr = "chr".*df[!, "chr"]
        GC_l = [(GC_dict[chr[i]][l_e[i]] - GC_dict[chr[i]][l_s[i]]) for i in 1:size(df)[1]]
        GC_r = [(GC_dict[chr[i]][r_e[i]] - GC_dict[chr[i]][r_s[i]]) for i in 1:size(df)[1]]
        AT_l = [(AT_dict[chr[i]][l_e[i]] - AT_dict[chr[i]][l_s[i]]) for i in 1:size(df)[1]]
        AT_r = [(AT_dict[chr[i]][r_e[i]] - AT_dict[chr[i]][r_s[i]]) for i in 1:size(df)[1]]
        GC_in = [(GC_dict[chr[i]][r_s[i]] - GC_dict[chr[i]][l_e[i]]) for i in 1:size(df)[1]]
        AT_in = [(AT_dict[chr[i]][r_s[i]] - AT_dict[chr[i]][l_e[i]]) for i in 1:size(df)[1]]
        Not_N = GC_l .+ AT_l
        Not_N[Not_N .== 0.0] .= 1.0
        GC_l = GC_l./Not_N
        Not_N = GC_r .+ AT_r
        Not_N[Not_N .== 0.0] .= 1.0
        GC_r = GC_r./Not_N
        Not_N = GC_in .+ AT_in
        Not_N[Not_N .== 0.0] .= 1.0
        GC_in = GC_in./Not_N
        df[!, "CG_frac_l"], df[!, "CG_frac_r"], df[!, "CG_frac_in"] = GC_l, GC_r, GC_in
    end

end


############################################################################################
                                                                                                