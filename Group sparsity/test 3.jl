using BenchmarkTools

function f()
    i=1
    while true
        if i==50
            break
        end

        i+=1
    end
end

function f₂()
    for i=1:100
        if i==50
            break
        end
    end
end

function repeval(f:: Function)
    for i in 1:10000
        f()
    end
end

@btime repeval(f)
@btime repeval(f₂)