"""
This script extracts level curves from large nu files.
"""

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

    incr = -500:10:500
    levelc = Dict((i,m) => "" for i in incr for m in (i >= 0 ? 1 : 1-i):(i>=0 ? 1500-i : 1500))
    @threads for i in incr
        for ((M,K,s),value) in nu 
            if K == 2 * (M-i)
                levelc[i,M-i] = value
            end
        end
    end
    println("Extracted level curves.")

    json_levelc = JSON.json(levelc)
    open("Data/p$(n_dps)/nu_levelc_M+i_2M_im500-10-500_p$(n_dps).json","w") do f
        write(f,json_levelc)
    end
    println("Exported JSON file.")
end

main()
