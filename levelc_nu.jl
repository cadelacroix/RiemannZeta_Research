using BenchmarkTools, JSON, .Threads

n_dps = 10_000
max_M = 1_300
increments = [0,1,2,10,20,50,100,200,500,800,1000]
nu_str = Dict{String,String}()
nu = Dict{Tuple{Integer,Integer,Integer},String}()


function parse_tuple(str)
    str = replace(str, '(' => "", ')' => "")
    parts = split(str, ", ")
    tuple = Tuple(parse(Int, part) for part in parts)
    
    tuple
end

function main()
    for i in 1:6
        file_content = read("/home/c.delacruz/JuliaYuri/M$(max_M)_p$(n_dps)/nu_N1-2-2599_K2-2-3000_s1_part$(i).json",String)
        data = JSON.parse(file_content)
        merge!(nu_str, data)
    end

    nu = Dict(parse_tuple(index) => value for (index,value) in nu_str)
    empty!(nu_str)

    @threads for incr in increments
        levelc = Dict{Integer,String}()
        for ((N,K,s),value) in nu 
            if N == K + 1 + incr
                levelc[K] = value
            end
        end
        json_levelc = JSON.json(levelc)
        open("nu_levelc_K+1+$(incr).json","w") do f
            write(f,json_levelc)
        end
    end
end

main()