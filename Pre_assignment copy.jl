using Distributions
# import Pkg; Pkg.add("Plots")
using Plots
plotly()
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
    for _ in 1:100000
        mylist[min(random_walk(), maxdim)] += 1
    end
    return mylist ./ 100000
end

# function probdf(maxdim)
#     x = []
#     y = []
#     for _ in 1:maxdim

# end

maxdim = 1000  # Walks higher than maxdim will be set at 10000
plotdim = 100   # To only plot first plotdim points


println("Starting...")
mylist = pdf(maxdim)
println("Running finished")

x_val = 2:plotdim
plot(x_val, mylist[2:plotdim])

# For fitting, exclude first values and all zero values
println("Fitting p(t) = A t^α")
y_fit = [el for el in mylist[2:maxdim-1] if el > 0]
x_fit = [float(i+1) for (i, el) in enumerate(mylist[2:maxdim-1])  if el > 0]
@show fitting = curve_fit(PowerFit, x_fit, y_fit)

y0 = fitting.(x_fit)
plot(x_fit, y0)
plot!(x_fit, y_fit)