module MainNetworkGrowth
    
    using LightGraphs, Distributions

    export Network_grow

    function Network_grow(delta::Float64, frac::Float64, fin_elem::Int64, in_comp::String, out_comp::String; lim_comp::String="nodes", start_cnt::Int64=2, nd_lim1::Int64=0, nd_lim2::Int64=1000000, rem_prob::Float64=0.0)
        gr = Graph(start_cnt, Int(start_cnt*(start_cnt - 1)/2))
        neis = ones(Int64, start_cnt)*(start_cnt - 1)
        last = start_cnt
        T = 0
        if lim_comp == "nodes"
            lim_val = start_cnt
        elseif lim_comp == "comps"
            lim_val = 1
        else lim_comp == "time"
            lim_val = 0
        end
        
        while (lim_val < fin_elem)
            event, dt = Choose_event(gr, neis, delta, in_comp)
            if event == 0
                last += 2
                add_vertices!(gr, 2)
                add_edge!(gr, last-1, last)
                append!(neis, [1, 1])
            else
                last += 1
                add_vertices!(gr, 1)
                gr, neis = Edges_inherit!(gr, neis, frac, event, last)
            end
            
            if rem_prob > 0
                rem_edges = filter(x -> rand() < rem_prob, collect(LightGraphs.edges(gr)))
                for e in rem_edges
                    rem_edge!(gr, e)
                    neis[e.src] -= 1
                    neis[e.dst] -= 1
                end
            end
            
            T += dt
            ccs = connected_components(gr)
            lens = map(length, ccs)
            if lim_comp == "nodes"
                lim_val = last          #  - length(ccs[lens .== 1])        # because of deletion
            elseif lim_comp == "comps"
                lim_val = length(ccs)   #  - length(ccs[lens .== 1])        # because of deletion
            else lim_comp == "time"
                lim_val = T
            end
        end
        
        if out_comp == "cc_distr"
            cc = connected_components(gr)
            lens = sort(map(length, cc))
            return lens[lens .> 1]                    # because of deletion
        elseif out_comp == "big_comp"
            cc = connected_components(gr)
            lens = map(length, cc)
            return maximum(lens)
        elseif out_comp == "node_deg"
            cc = connected_components(gr)
            lens = sort(map(length, cc))
            max_size = lens[end]
            if nd_lim1 < max_size < nd_lim2
                nd = map(x -> length(neighbors(gr, x)), vertices(gr))
            else
                nd = [1, 1]
            end
            return sort(nd[nd .> 0])                        # because of deletion
        elseif out_comp == "graph"
            vr = collect(vertices(gr))
            nd = map(x -> length(neighbors(gr, x)), vr)
            map(x -> rem_vertex!(gr, x), vr[nd .== 0])      # because of deletion
            return gr, T
        elseif out_comp == "com_nodes_edges"
            cc = connected_components(gr)
            lens = sort(map(length, cc))
            cc = cc[lens .> 1]                               # because of deletion
            nodes = map(length, cc)
            edges = map(x -> sum([length(neighbors(gr, i)) for i in x])/2, cc)
            nodes_m = sort(unique(nodes))
            edges_m = map(x -> mean(edges[nodes .== x]), sort(unique(nodes)))
            return nodes_m, edges_m
        end
        
    end

    function Choose_event(gr::LightGraphs.SimpleGraph, neis::Array{Int64, 1}, delta::Float64, in_comp::String)
        if in_comp == "mod1"
            rates = delta*ones(length(neis))
        elseif in_comp == "mod2"
            rates = delta*neis
        end
        push!(rates, 1.0)
        events = collect(vertices(gr))
        push!(events, 0)
        event = wsample(events, rates)
        dt = log(1/rand())/sum(rates)
        return event, dt
    end
    
    function Edges_inherit!(gr::LightGraphs.SimpleGraph, neis::Array{Int64, 1}, frac::Float64, event::Int64, last::Int64)
        edges = neighbors(gr, event)
        inh_edges = filter(x -> rand() < frac, edges)
        push!(inh_edges, event)
        map(x -> add_edge!(gr, x, last), inh_edges)
        map(x -> neis[x] += 1, inh_edges)
        push!(neis, length(inh_edges))
        return gr, neis
    end
end