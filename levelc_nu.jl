using BenchmarkTools, JSON, .Threads

function parse_tuple(str)
    str = replace(str, '(' => "", ')' => "")
    parts = split(str, ", ")
    tuple = Tuple(parse(Int, part) for part in parts)
    
    tuple
end

function main()
    n_dps = 10_000
    file_content = Dict(i => "" for i in 1:16)
    nu_str = Dict(i => Dict() for i in 1:16)
    nu_part = Dict(i => Dict() for i in 1:16)
    @threads for i in 1:16
       file_content[i] = read("Data/p$(n_dps)/nu_M1-1-1500_K1-1-3000_s1-1-1_p$(n_dps)_part_$(i).json",String)
       nu_str[i] = JSON.parse(file_content[i])
       nu_part[i] = Dict(parse_tuple(index) => value for (index,value) in nu_str[i])
    end
    nu = merge(values(nu_part)...)
    println("Opened files.") 

    levelc = Dict{Integer,String}()
    for ((M,K,s),value) in nu 
        if M == K
            levelc[K] = value
        end
    end
    println("Extracted level curve.")

    json_levelc = JSON.json(levelc)
    open("Data/p$(n_dps)/nu_levelc_MeqK_p$(n_dps).json","w") do f
        write(f,json_levelc)
    end
    println("Exported JSON file.")
end

main()
