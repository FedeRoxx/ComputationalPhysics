using Distributions
# import Pkg; Pkg.add("Plots")
using Plots
using CurveFit

function random_walk()
    step = rand(Uniform(-1, 1))
    x=step
    t = 1
    while x * sign(step) >= 0
        x += rand(Uniform(-1, 1))
        t += 1
    end
    return t
end

function pdf(maxdim)
    mylist = [0 for _ in 1:maxdim ]
    for _ in 1:1000000
        mylist[min(random_walk(), maxdim)] += 1
    end
    return mylist ./ 1000000
end



maxdim = 100
mylist = pdf(maxdim)


plotdim = 99

x_val = 2:plotdim

plot(x_val, mylist[2:plotdim])

y_fit = [el for el in mylist[2:maxdim-1] if el > 0]
x_fit = [float(i+1) for (i, el) in enumerate(mylist[2:maxdim-1])  if el > 0]

@show a = curve_fit(PowerFit, x_fit, y_fit)

y0 = a.(x_val)
plot(x_val, y0)
plot!(x_val, mylist[x_val])
# plot!(x_val,  mylist[x_val])
# println(a)