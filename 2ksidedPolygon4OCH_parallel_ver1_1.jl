using LinearAlgebra
using Random
using Base.Threads
include("utils.jl")
@inline function starting_vertices_parallel_ver1_1(X,k)
    Is = [[] for _ in 1:(4k+1)]
    len = size(X,1)

    @threads for j in 1:2k
        if j == 1
            d = pi*(j-1)/2k
            e_x = cos(d)
            e_y = sin(d)
            max_j = e_x*X[1,1]+e_y*X[1,2]
            p_max = [X[1,1],X[1,2]]
            p_min = [X[1,1],X[1,2]]

            min_j = e_x*X[1,1]+e_y*X[1,2]
            p_max2 = [X[1,1],X[1,2]]
            p_min2 = [X[1,1],X[1,2]]

            for i in 2:len
                x = X[i,1]
                y = X[i,2]
                t = e_x*x+ e_y*y
                #I^j
                max_j_eq_t = max_j == t
                max_j_lt_t = max_j < t

                y_gt_y_max = y > p_max[2]
                y_lt_y_min = y < p_min[2]

                max_j = max_j_lt_t ? t : max_j
                p_max[1] = max_j_lt_t ? x : p_max[1]
                p_min[1] = max_j_lt_t ? x : p_min[1]
                p_max[2] = max_j_lt_t ? y : p_max[2]
                p_min[2] = max_j_lt_t ? y : p_min[2]

                p_min[1] = (max_j_eq_t & y_lt_y_min) ? x : p_min[1]
                p_min[2] = (max_j_eq_t & y_lt_y_min) ? y : p_min[2]
                p_max[1] = (max_j_eq_t & y_gt_y_max) ? x : p_max[1]
                p_max[2] = (max_j_eq_t & y_gt_y_max) ? y : p_max[2]

                #I^(j+2k)
                min_j_eq_t = min_j == t 
                min_j_gt_t =  min_j > t
                
                y_gt_y_max2 = y > p_max2[2]
                y_lt_y_min2 = y < p_min2[2]

                min_j = min_j_gt_t ? t : min_j
                p_min2[1] = min_j_gt_t ? x : p_min2[1]
                p_max2[1] = min_j_gt_t ? x : p_max2[1]
                p_min2[2] = min_j_gt_t ? y : p_min2[2]
                p_max2[2] = min_j_gt_t ? y : p_max2[2]

                p_min2[1] = (min_j_eq_t & y_lt_y_min2) ? x : p_min2[1]
                p_min2[2] = (min_j_eq_t & y_lt_y_min2) ? y : p_min2[2]
                p_max2[1] = (min_j_eq_t & y_gt_y_max2) ? x : p_max2[1]
                p_max2[2] = (min_j_eq_t & y_gt_y_max2) ? y : p_max2[2]
        
            end
            Is[j] = [p_min,p_max]
            Is[j+2k] = [p_max2,p_min2]
        elseif j == k+1
            d = pi*(j-1)/2k
            e_x = cos(d)
            e_y = sin(d)
            max_j = e_x*X[1,1]+e_y*X[1,2]
            p_max = [X[1,1],X[1,2]]
            p_min = [X[1,1],X[1,2]]

            min_j = e_x*X[1,1]+e_y*X[1,2]
            p_max2 = [X[1,1],X[1,2]]
            p_min2 = [X[1,1],X[1,2]]
            for i in 2:len
                x = X[i,1]
                y = X[i,2]
                t = e_x*x+ e_y*y
                #I^j
                max_j_eq_t = max_j == t
                max_j_lt_t = max_j < t

                x_gt_x_max = x > p_max[1]
                x_lt_x_min = x < p_min[1]

                max_j = max_j_lt_t ? t : max_j
                p_max[1] = max_j_lt_t ? x : p_max[1]
                p_min[1] = max_j_lt_t ? x : p_min[1]
                p_max[2] = max_j_lt_t ? y : p_max[2]
                p_min[2] = max_j_lt_t ? y : p_max[2]

                p_min[1] = (max_j_eq_t & x_lt_x_min) ? x : p_min[1]
                p_min[2] = (max_j_eq_t & x_lt_x_min) ? y : p_min[2]
                p_max[1] = (max_j_eq_t & x_gt_x_max) ? x : p_max[1]
                p_max[2] = (max_j_eq_t & x_gt_x_max) ? y : p_max[2]

                #I^(j+2k)
                min_j_eq_t = min_j == t 
                min_j_gt_t =  min_j > t
                
                x_gt_x_max2 = x > p_max2[1]
                x_lt_x_min2 = x < p_min2[1]

                min_j = min_j_gt_t ? t : min_j
                p_min2[1] = min_j_gt_t ? x : p_min2[1]
                p_max2[1] = min_j_gt_t ? x : p_max2[1]
                p_min2[2] = min_j_gt_t ? y : p_min2[2]
                p_max2[2] = min_j_gt_t ? y : p_max2[2]

                p_min2[1] = (min_j_eq_t & x_lt_x_min2) ? x : p_min2[1]
                p_min2[2] = (min_j_eq_t & x_lt_x_min2) ? y : p_min2[2]
                p_max2[1] = (min_j_eq_t & x_gt_x_max2) ? x : p_max2[1]
                p_max2[2] = (min_j_eq_t & x_gt_x_max2) ? y : p_max2[2]
            end
            Is[j] = [p_max,p_min]
            Is[j+2k] = [p_min2,p_max2]

        else
            d = pi*(j-1)/2k
            e_x = cos(d)
            e_y = sin(d)
            max_j = e_x*X[1,1]+e_y*X[1,2]
            p = [X[1,1],X[1,2]]

            min_j = e_x*X[1,1]+e_y*X[1,2]
            p2 = [X[1,1],X[1,2]]
            for i in 2:len
                x = X[i,1]
                y = X[i,2]
                t = e_x*x+ e_y*y
                max_j_lt_t = max_j < t
                max_j = max_j_lt_t ? t : max_j
                p[1] = max_j_lt_t ? x : p[1]
                p[2] = max_j_lt_t ? y : p[2]

                min_j_gt_t =  min_j > t

                min_j = min_j_gt_t ? t : min_j
                p2[1] = min_j_gt_t ? x : p2[1]
                p2[2] = min_j_gt_t ? y : p2[2]
            end
            Is[j] = [p]
            Is[j+2k] = [p2]
        end          
    end
    Is[4k+1] = [Is[1][1]]
    return Is
end

@inline function find_P_parallel_ver1_1(X, Is,k)
    Xx = X[:,1]
    Xy = X[:,2]
    Pp = [zeros(0,2) for _ in 1:4k]
    @threads for i in 1:4k 
        if i == 1
            Pp[i] = vcat(X[(Xx .> Is[i+1][1][1]) .& (Xy .> Is[i][2][2]),:],[Is[i][2][1] Is[i][2][2];Is[i+1][1][1] Is[i+1][1][2]])
        elseif i > 1 && i <= k
            Pp[i] = vcat(X[(Xx .> Is[i+1][1][1]) .& (Xy .> Is[i][1][2]),:],[Is[i][1][1] Is[i][1][2];Is[i+1][1][1] Is[i+1][1][2]])
        elseif i == k+1
            Pp[i] = vcat(X[(Xx .< Is[i][2][1]) .& (Xy .> Is[i+1][1][2]),:],[Is[i][2][1] Is[i][2][2];Is[i+1][1][1] Is[i+1][1][2]])
        elseif  i > k+1 && i<=2k
            Pp[i] = vcat(X[(Xx .< Is[i][1][1]) .& (Xy .> Is[i+1][1][2]),:],[Is[i][1][1] Is[i][1][2];Is[i+1][1][1] Is[i+1][1][2]])
        elseif i == 2k+1
            Pp[i] = vcat(X[(Xx .< Is[i+1][1][1]) .& (Xy .< Is[i][2][2]),:],[Is[i][2][1] Is[i][2][2];Is[i+1][1][1] Is[i+1][1][2]])
        elseif i > 2k+1 && i <=3k
            Pp[i] = vcat(X[(Xx .< Is[i+1][1][1]) .& (Xy .< Is[i][1][2]),:],[Is[i][1][1] Is[i][1][2];Is[i+1][1][1] Is[i+1][1][2]])
        elseif i == 3k+1
            Pp[i] = vcat(X[(Xx .> Is[i][2][1]) .& (Xy .< Is[i+1][1][2]),:],[Is[i][2][1] Is[i][2][2];Is[i+1][1][1] Is[i+1][1][2]])
        else
            Pp[i] = vcat(X[(Xx .> Is[i][1][1]) .& (Xy .< Is[i+1][1][2]),:],[Is[i][1][1] Is[i][1][2];Is[i+1][1][1] Is[i+1][1][2]])
        end
    end
    return Pp
end

@inline function o_1_parallel_ver1_1(P)
    P = sortslices(P, dims = 1, by = x -> -x[1])
    ohull = [P[1,:]]
    x_min = P[1,1]
    y_max = P[1,2]
    len = size(P,1)
    for i in 2:len
        x = P[i,1]
        y = P[i,2]
        if y > y_max
            if x != x_min
                push!(ohull,P[i,:])
                x_min = x
                y_max = y
            else
                ohull[end] = P[i,:]
                y_max = y
            end
        end
    end
    return ohull
end

@inline function o_2_parallel_ver1_1(P)
    P = sortslices(P, dims = 1, by = x -> x[1])
    ohull = [P[1,:]]
    x_max = P[1,1]
    y_max = P[1,2]
    len = size(P,1)
    for i in 2:len
        x = P[i,1]
        y = P[i,2]
        if y > y_max
            if x != x_max
                pushfirst!(ohull,P[i,:])
                x_max = x
                y_max = y
            else
                ohull[1] = P[i,:]
                y_max = y
            end
        end
    end
    return ohull
end

@inline function o_3_parallel_ver1_1(P)
    P = sortslices(P, dims = 1, by = x -> x[1])
    ohull = [P[1,:]]
    x_max = P[1,1]
    y_min = P[1,2]
    len = size(P,1)
    for i in 2:len
        x = P[i,1]
        y = P[i,2]
        if y < y_min 
            if x != x_max
                push!(ohull,P[i,:])
                x_max = x
                y_min = y
            else
                ohull[end] = P[i,:]
                y_min = y
            end
        end
    end
    return ohull
end

@inline function o_4_parallel_ver1_1(P)
    P = sortslices(P, dims = 1, by = x -> -x[1])
    ohull = [P[1,:]]
    x_min = P[1,1]
    y_min = P[1,2]
    len = size(P,1)
    for i in 2:len
        x = P[i,1]
        y = P[i,2]
        if y < y_min 
            if x != x_min
                pushfirst!(ohull,P[i,:])
                x_min = x
                y_min = y
            else
                ohull[1] = P[i,:]
                y_min = y
            end 
        end
    end
    return ohull
end

@inline function alg_2ksidePolygon4OCH_parallel_ver1_1(X,k)
    Is = starting_vertices_parallel_ver1_1(X,k)
    P = find_P_parallel_ver1_1(X,Is,k)
    V = [[] for _ in 1:4k]
    @threads for i in 1:4k
        if i > 0 && i <= k
            V[i] = o_1_parallel_ver1_1(P[i])
        elseif i > k && i <= 2k
            V[i] = o_2_parallel_ver1_1(P[i])
        elseif i > 2k && i <= 3k
            V[i] = o_3_parallel_ver1_1(P[i])
        else
            V[i] = o_4_parallel_ver1_1(P[i])
        end

    end
    ohull = V[1]
    for i in 2:4k
        if V[i][1] == V[i-1][end]
            append!(ohull,V[i][2:end])
        else
            append!(ohull,V[i])
        end
    end 
    if ohull[1] == ohull[end]
        pop!(ohull)
    end
    return ohull
end

@inline function alg_2ksidePolygon4OCH_parallel_ex_ver1_1(X,k, exportFile)
    ohull = alg_2ksidePolygon4OCH_parallel_ver1_1(X,k)
    exportResult(ohull, exportFile)
end

X = create_discs(1000)
ohull = @time alg_2ksidePolygon4OCH_parallel_ver1_1(X,8)
