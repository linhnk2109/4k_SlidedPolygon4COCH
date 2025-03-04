using Random
using BenchmarkTools
include("2ksidedPolygon4OCH.jl") 
include("2ksidedPolygon4OCH_parallel_ver1_1.jl") 
include("2ksidedPolygon4OCH_parallel_ver1_2.jl") 
include("2ksidedPolygon4OCH_parallel_ver2_1.jl") 
include("2ksidedPolygon4OCH_parallel_ver2_2.jl") 

include("O_Graham.jl")
include("O_Quickhull.jl")


function main()
    sizes = [1_000_000, 10_000_000, 20_000_000, 30_000_000, 40_000_000, 50_000_000, 60_000_000, 70_000_000, 80_000_000, 90_000_000, 100_000_000] 
    setK = [2, 3, 5, 8, 10, 12]
   
    
    setNumbers = 1
    resultDirectory = "result/"
    # set benchmarking = true if want to benchmark
    benchmarking = true
    # benchmarking = false
    # only export results if exportResult = true and benchmarking = false
    exportResult = true
    # dataType: 
        #1 for discs type, 
        #2 for hollow discs type
        #3 for square type
        #4 for hollow square type
        #5 for sun type
        #6 for hollow sun type
        #7 for circles type
    dataType = 7
    dataTypeName = ["Discs", "HollowDiscs","Square","HollowSquare","Sun","HollowSun","Circles"]
    
    instanceNames = Vector{String}(undef, length(sizes)*setNumbers)

    runningTime2kOCH = [Matrix{Float64}(undef, length(instanceNames),7) for _ in 1:length(setK)]
    runningTime2kOCHParallel_ver1_1 = [Matrix{Float64}(undef, length(instanceNames),7) for _ in 1:length(setK)]
    runningTime2kOCHParallel_ver1_2 = [Matrix{Float64}(undef, length(instanceNames),7) for _ in 1:length(setK)]
    runningTime2kOCHParallel_ver2_1 = [Matrix{Float64}(undef, length(instanceNames),7) for _ in 1:length(setK)]
    runningTime2kOCHParallel_ver2_2 = [Matrix{Float64}(undef, length(instanceNames),7) for _ in 1:length(setK)]
    runningTimeOGraham = Matrix{Float64}(undef, length(instanceNames),7)
    runningTimeOGrahamOp = Matrix{Float64}(undef, length(instanceNames),7)
    runningTimeOQhull = Matrix{Float64}(undef, length(instanceNames),7)
    runningTimeOQhullOp = Matrix{Float64}(undef, length(instanceNames),7)

    Random.seed!(42)
    for i in 1:length(sizes)
        for j in 1:setNumbers
            k = (i-1)*setNumbers+j
            instanceNames[k] = string(dataTypeName[dataType], "_", sizes[i], "_", j)

            println()
            println("Consider instance ", instanceNames[k])

            # create random data
            points = Matrix{Float64}(undef, sizes[i], 2)
            if dataType == 1
                points = create_discs(sizes[i])
            elseif dataType == 2
                points = create_hollowDiscs(sizes[i])
            elseif dataType == 3
                points = create_square(sizes[i])
            elseif dataType == 4
                points = create_hollowSquare(sizes[i])
            elseif dataType == 5
                points = create_sun(sizes[i])
            elseif dataType == 6
                points = create_hollowSun(sizes[i])
            else
                points = create_circles(sizes[i])
            end

            if benchmarking
                lenk = length(setK)
                for l in 1:lenk
                    K_ = setK[l]
                    println("2k_Sided_Polygon_4_OCH; k = ",K_)
                    bm2kOCH = run(@benchmarkable alg_2ksidePolygon4OCH($points,$K_) samples=7 seconds=10000)
                    runningTime2kOCH[l][k,:] = report(bm2kOCH, 2)

                    println("4k_Sided_Polygon_4_OCH_Parallel_ver1_1; k = ",K_)
                    bm2kOCHParallel_ver1_1 = run(@benchmarkable alg_2ksidePolygon4OCH_parallel_ver1_1($points,$K_) samples=7 seconds=10000)
                    runningTime2kOCHParallel_ver1_1[l][k,:] = report(bm2kOCHParallel_ver1_1, 2)

                    println("2k_Sided_Polygon_4_OCH_Parallel_ver1_2; k = ",K_)
                    bm2kOCHParallel_ver1_2 = run(@benchmarkable alg_2ksidePolygon4OCH_parallel_ver1_2($points,$K_) samples=7 seconds=10000)
                    runningTime2kOCHParallel_ver1_2[l][k,:] = report(bm2kOCHParallel_ver1_2, 2)

                    println("2k_Sided_Polygon_4_OCH_Parallel_ver2_1; k = ",K_)
                    bm2kOCHParallel_ver2_1 = run(@benchmarkable alg_2ksidePolygon4OCH_parallel_ver2_1($points,$K_) samples=7 seconds=10000)
                    runningTime2kOCHParallel_ver2_1[l][k,:] = report(bm2kOCHParallel_ver2_1, 2)

                    println("2k_Sided_Polygon_4_OCH_Parallel_ver2_2; k = ",K_)
                    bm2kOCHParallel_ver2_2 = run(@benchmarkable alg_2ksidePolygon4OCH_parallel_ver2_2($points,$K_) samples=7 seconds=10000)
                    runningTime2kOCHParallel_ver2_2[l][k,:] = report(bm2kOCHParallel_ver2_2, 2)

                end
                
                println("O_Quickhull")
                bmOQhull = run(@benchmarkable find_o_quickhull($points) samples=7 seconds=10000)
                runningTimeOQhull[k,:] = report(bmOQhull, 2)
                
                println("O_Graham")
                bmOGraham = run(@benchmarkable find_hull_o_graham($points) samples=7 seconds=10000)
                runningTimeOGraham[k,:] = report(bmOGraham, 2)

             else
                
                 exportFileOQhullParallel = string(resultDirectory, instanceNames[k], "_OQHullParallel")
                 find_o_quickhull_parallel_nr_ex(points, exportFileOQhullParallel)

                 exportFileOQhull = string(resultDirectory, instanceNames[k], "_OQhull")
                 find_o_quickhull_ex(points, exportFileOQhull)
            end
        end
    end

    if benchmarking
        baseName = string(resultDirectory, dataTypeName[dataType], "_")
        for l in 1:length(setK)
            exportReport(instanceNames, runningTime2kOCH[l], string(baseName, "4kOCH_", setK[l], "_running_time.csv"))
            exportReport(instanceNames, runningTime2kOCHParallel_ver1_1[l], string(baseName, "4kOCH_", setK[l], "_running_time.csv"))
            exportReport(instanceNames, runningTime2kOCHParallel_ver1_2[l], string(baseName, "2ksidedPolygon4OCHParallel_ver1_2_", setK[l], "_running_time.csv"))
            exportReport(instanceNames, runningTime2kOCHParallel_ver2_1[l], string(baseName, "2ksidedPolygon4OCHParallel_ver2_1_", setK[l], "_running_time.csv"))
            exportReport(instanceNames, runningTime2kOCHParallel_ver2_2[l], string(baseName, "2ksidedPolygon4OCHParallel_ver2_2_", setK[l], "_running_time.csv"))
        end
        exportReport(instanceNames, runningTimeOQhull, string(baseName,"OQhull.csv"))
        exportReport(instanceNames, runningTimeOGraham, string(baseName,"OGraham.csv"))
    end
end

main()