using BenchmarkTools, JSON, .Threads

function parse_tuple(str)
    str = replace(str, '(' => "", ')' => "")
    parts = split(str, ", ")
    tuple = Tuple(parse(Int, part) for part in parts)
    
    tuple
end

function main()
    n_dps = 40_000

    file_content = read("Data/p$(n_dps)/nu_M1-1-200_K1-1-600_s1-1-1_p$(n_dps).json",String)
    nu_str = JSON.parse(file_content)
    nu = Dict(parse_tuple(index) => value for (index,value) in nu_str)
    println("Opened file.")
    empty!(nu_str)

    levelc = Dict{Integer,String}()
    for ((M,K,s),value) in nu 
        if M == Int(ceil(K/2))
            levelc[K] = value
        end
    end
    println("Extracted level curve.")

    json_levelc = JSON.json(levelc)
    open("Data/p$(n_dps)/nu_levelc_MceilKdiv2_p$(n_dps).json","w") do f
        write(f,json_levelc)
    end
    println("Exported JSON file.")
end

main()
