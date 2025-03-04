using LinearAlgebra
# using CSV
# using DataFrames
function find_special_points_graham(points)
    maxY = maximum(points[:,2])
    minY = minimum(points[:,2])
    maxX = maximum(points[:,1])
    minX = minimum(points[:,1])

    rightPoints = points[findall(row -> row[1]==maxX,eachrow(points)),:]
    leftPoints = points[findall(row -> row[1]==minX,eachrow(points)),:]
    topPoints = points[findall(row -> row[2] == maxY,eachrow(points)),:]
    bottomPoints = points[findall(row -> row[2] == minY,eachrow(points)),:]

    if size(topPoints, 1) == 1
        top = [topPoints[1, :]]
    else
        idx = sortperm(topPoints[:, 1])
        topPoints = topPoints[idx,:]
        top = [topPoints[1, :], topPoints[end, :]]
    end
    
    if size(bottomPoints, 1) == 1
        bottom = [bottomPoints[1, :]]
    else
        idx = sortperm(bottomPoints[:,1],rev = true)
        bottomPoints = bottomPoints[idx,:]
        bottom = [bottomPoints[1, :], bottomPoints[end, :]]
    end
    
    if size(rightPoints, 1) == 1
        right = [rightPoints[1, :]]
    else
        idx = sortperm(rightPoints[:,2],rev = true)
        rightPoints = rightPoints[idx,:]
        right = [rightPoints[1, :], rightPoints[end, :]]
    end
    
    if size(leftPoints, 1) == 1
        left = [leftPoints[1, :]]
    else
        idx = sortperm(leftPoints[:,2])
        leftPoints = leftPoints[idx,:]
        left = [leftPoints[1, :], leftPoints[end, :]]
    end
    
    if length(top) == 1
        p2_1 = p2_2 = top[1]
    else
        p2_2 = top[1]
        p2_1 = top[2]
    end
    p1_2 = right[1]
    if length(right) == 1
        p1_1 = right[1]
    else
        p1_1 = right[2]
    end
    p4_2 = bottom[1]
    if length(bottom) == 1
        p4_1 = bottom[1]
    else
        p4_1 = bottom[2]
    end
    p3_2 = left[1]
    if length(left) == 1
        p3_1 = left[1]
    else
        p3_1 = left[2]
    end
    return p1_1,p1_2, p2_1,p2_2, p3_1,p3_2, p4_1,p4_2
end

function find_sets_graham(points)
    p1_1,p1_2, p2_1,p2_2, p3_1,p3_2, p4_1,p4_2 = find_special_points_graham(points)
    set1 = points[(points[:,1] .>= p2_1[1]) .& (points[:,2] .>=p1_2[2]),:]
    set2 = points[(points[:,1] .<= p2_2[1]) .& (points[:,2] .>=p3_1[2]),:]
    set3 = points[(points[:,1] .<= p4_1[1]) .& (points[:,2] .<=p3_2[2]),:]
    set4 = points[(points[:,1] .>= p4_2[1]) .& (points[:,2] .<=p1_1[2]),:]
    return set1,set2,set3,set4
end

function o_graham_1(set1)
    set1 = sortslices(set1, dims = 1 , by = x -> (-x[1],-x[2]))
    b = [true]
    for i in 1:(size(set1,1)-1)
        if set1[i,1] == set1[i+1,1]
            append!(b, false)
        else
            append!(b, true)
        end
    end
    set1 = set1[b,:]
    return set1

end
function o_graham_2(set2)
    set2 = sortslices(set2, dims = 1, by = x -> (-x[2], x[1]))
    b = [true]
    for i in 1:(size(set2,1)-1)
        if set2[i,2] == set2[i+1,2]
            append!(b, false)
        else
            append!(b, true)
        end
    end
    set2 = set2[b,:]
    return set2
end
function o_graham_3(set3)
    set3 = sortslices(set3, dims = 1, by = x -> (x[1], x[2]))
    b = [true]
    for i in 1:(size(set3,1)-1)
        if set3[i,1] == set3[i+1,1]
            append!(b, false)
        else
            append!(b, true)
        end
    end
    set3 = set3[b,:]
    return set3
end
function o_graham_4(set4)
    set4 = sortslices(set4, dims = 1, by = x -> (x[2], -x[1]))
    b = [true]
    @inbounds @simd for i in 1:(size(set4,1)-1)
        if set4[i,1] == set4[i+1,1]
            append!(b, false)
        else
            append!(b, true)
        end
    end
    set4 = set4[b,:]
    return set4
end

function find_hull_o_graham(points)
    set1,set2,set3,set4 = find_sets_graham(points)
    set1 = o_graham_1(set1)
    set2 = o_graham_2(set2)
    set3 = o_graham_3(set3)
    set4 = o_graham_4(set4)
    hull = [set1[1,:]]
    len1 = size(set1,1)
    for i in 1:len1 
        cand = hull[end]
        if set1[i,2] > cand[2]
            push!(hull,set1[i,:])
        end
    end
    len2 = size(set2,1)
    for i in 1:len2 
        cand = hull[end]
        if set2[i,1] < cand[1]
            push!(hull, set2[i,:])
        end
    end

    len3 = size(set3,1)
    for i in 1:len3 
        cand = hull[end]
        if set3[i,2] < cand[2]
            push!(hull,set3[i,:])
        end
    end
    len4 = size(set4,1)
    for i in 1:len4 
        cand = hull[end]
        if set4[i,1] > cand[1]
            push!(hull, set4[i,:])
        end
    end
    push!(hull, hull[1])
    return hull
end

function find_hull_o_graham_ex(points, exportFile)
    hull = find_hull_o_graham(points)
    exportResult(hull, exportFile)
end

# using Plots
include("utils.jl")
points = create_discs(1000)
hull = @time find_hull_o_graham(points)
# unique(hull)
# scatter(points[:,1],points[:,2],color=:red, markersize=2,legend=false)
# x = []
# y = []
# for i in 1:length(hull) 
#     push!(x,hull[i][1])
#     push!(y,hull[i][2])
# end
# s =Vector{Tuple{Int, Vector{Float64}}}(undef,0)
# for i in 1:(length(hull)-1)
#     if hull[i+1][1] > hull[i][1] && hull[i+1][2] > hull[i][2]
#         push!(s,(i,[hull[i][1],hull[i+1][2]]))
#     elseif hull[i+1][1] > hull[i][1] && hull[i+1][2] < hull[i][2]
#         push!(s,(i,[hull[i+1][1],hull[i][2]]))
#     elseif hull[i+1][1] < hull[i][1] && hull[i+1][2] < hull[i][2]
#         push!(s,(i,[hull[i][1],hull[i+1][2]]))
#     elseif hull[i+1][1] < hull[i][1] && hull[i+1][2] > hull[i][2]
#         push!(s,(i,[hull[i+1][1],hull[i][2]]))        
#     end
# end

# for i in 1:length(s) 
#     insert!(x,s[i][1]+i,s[i][2][1])
#     insert!(y,s[i][1]+i,s[i][2][2])
# end

# plot!(x,y,color=:green, markersize=8)

