module SimulGrowth
    
    using StatsBase, Plots, Histograms, LightGraphs, MetaGraphs
                                                            
    export SystemStates, Run_ss_object_growth, main_growth_PL

    pyplot()
                
    mutable struct SystemStates
        t::Float64
        nodes::Array{Float64,1}
    end
    
    function SS_initialize()
        ss_obj = SystemStates(0.0, Float64[2.0]) 
    end                                
    
    function SS_time_step(ss_obj::SystemStates, step::Float64, delta::Float64, frac::Float64, comp::String="int")
        ss_obj.t += step
        nodes_prew = ss_obj.nodes[end]
        nodes_exp = step*delta*(nodes_prew^(1.0+frac))
        if comp == "int"
            nodes_add1, nodes_prob = Split_to_Int_Float(nodes_exp)
            nodes_add2 = sample([0.0, 1.0], Weights([1.0 - nodes_prob, nodes_prob]))
            nodes_new = nodes_prew + nodes_add1 + nodes_add2
            push!(ss_obj.nodes, nodes_new)
        elseif comp == "float"
            nodes_new = nodes_prew + nodes_exp
            push!(ss_obj.nodes, nodes_new)
        end
        return ss_obj
    end
    
    function Main_evolve_ss_object(time::Float64, delta::Float64, frac::Float64, iter_num::Int)
        t_steps, t_last = SimulGrowth.Split_to_Int_Float(time)
        out_arr::Array{Float64,1} = []
        for iter in 1:iter_num
            ss_objects = []
            for i in 1:t_steps
                for (j, ss) in enumerate(ss_objects)
                    ss_objects[j] = SimulGrowth.SS_time_step(ss, 1.0, delta, frac, "int")
                end
                ss_new = SimulGrowth.SS_initialize()
                push!(ss_objects, ss_new)   
            end

            for (j, ss) in enumerate(ss_objects)
                ss_objects[j] = SimulGrowth.SS_time_step(ss, t_last, delta, frac, "int")
            end
            append!(out_arr, SimulGrowth.Extract_nodes_from_array(ss_objects))
        end
        return out_arr
    end
                                                            
    function Main_evolve_ss_object(node_stop::Float64, delta::Float64, frac::Float64, t_step::Float64=0.1)
        ss = SS_initialize()
        while ss.nodes[end] < node_stop
            ss = SS_time_step(ss, t_step, delta, frac, "float")
        end
        return ss
    end
    
    function Run_ss_object_growth(delta::Float64, frac::Float64, node_stop::Float64=1000.0, t_step::Float64=0.1)
        con = (-1/frac)*2^(-frac)
        s_obj = Main_evolve_ss_object(node_stop, delta, frac, t_step)
        ts = collect(0.0:t_step:s_obj.t)
        
        ser_plt = plot(log10.(ts), log10.(s_obj.nodes), markershape=:auto, line=false)
        pred_nodes = []
        ts2 = []
        for t in ts
            if -frac*con .<= delta*frac*t
                break
            end
            push!(pred_nodes, (-frac*con - delta*frac*t)^(-1/frac))
            push!(ts2, t)
        end
        ser_plt = plot!(log10.(ts2), log10.(pred_nodes))
        display(ser_plt)

        cc_h = Histograms.Histogram(s_obj.nodes, N=40, scale=:log10);
        cc_p = plot(cc_h, xscale=:log10, yscale=:log10, markershape=:auto, line=false);
        x = range(1.0, 10.0^3.0, length=10); y = 0.18*x.^(-2);
        cc_p = plot!(x, y, label="a = -2");
        x = range(1.0, 10.0^3.0, length=10); y = 0.18*x.^(-1.0 - frac);
        cc_p = plot!(x, y, label="a = $(-1.0 - frac)");
        x = range(1.0, 10.0^3.0, length=10); y = 0.18*x.^(-1);
        cc_p = plot!(x, y, label="a = -1")
        display(cc_p)

        cc_count = countmap(floor.(s_obj.nodes))
        cc_size = sort(collect(keys(cc_count)))
        cc_vals = map(x -> cc_count[x], cc_size)
        p = plot(log10.(cc_size), log10.(cc_vals), markershape=:auto, line=false)
        x = range(0, 2.5, length=10); y = -2*x .+ 5.6;
        p = plot!(x, y, label="a = -2")
        x = range(0, 2.5, length=10); y = (-1.0 - frac)*x .+ 3.6;
        p =plot!(x, y, label="a = $(-1.0 - frac)")
        gui(p)
    end
    
    function Split_to_Int_Float(val::Float64)
        out_main = floor(val)
        out_rest = val - floor(val)
        [out_main, out_rest]
    end
    
    function Extract_nodes_from_array(in_arr::Array)                                                        
        out_arr = [ss.nodes[end] for ss in in_arr]
        out_arr = filter(x -> (x < Inf) && (!isnan(x)), out_arr)
        return out_arr
    end
    
    function simul_PL(T::Int64, pls::Array{Int64,1})
        objs = Int64[1]

        for t in 1:T
            inds = Int64[]
            pl = rand(pls)
            for i in 1:pl
                push!(inds, sample(1:length(objs), Weights(objs))[1])
            end
            map(x -> objs[x] += 1, inds)
            push!(objs, 1)
        end
        return objs
    end
    
    function main_growth_PL(T::Int64, pls::Array{Int64,1})
        objs = SimulGrowth.simul_PL(T, pls)
        cc_h = Histograms.Histogram(objs, N=20, scale=:log10)
        cc_p = plot(cc_h, xscale=:log10, yscale=:log10, markershape=:auto, line=false)
        x = range(1, 1000, length=10); y = 2*x.^(-2.0);
        cc_p = plot!(x, y, label="a = -2")
        gui(cc_p)
    end
    
end