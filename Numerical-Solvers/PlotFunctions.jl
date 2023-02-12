function PlotNucleusPlot(solutions)
    p = scatter(legend=false)
    for X = solutions
        scatter!(p, X[1,[1]], X[2,[1]], X[3,[1]], color=:blue)
        scatter!(p, X[1,[N]], X[2,[N]], X[3,[N]], color=:red)
        plot!(p, X[1,[1,N]], X[2,[1,N]], X[3,[1,N]], color=:black)
    end
    display(p)
end

function PlotTrajectories(solutions)
    p = plot(legend=false)
    for X = solutions
        plot!(p, X[1,:], X[2,:], X[3,:]) #, color=:blue) 
    end
    display(p)
end

function PlotTrajectoriesInterpolated(solutions)
    p = plot(legend=false)
    for X = solutions
        plot!(p, X, vars=(1,2,3)) #, color=:blue) 
    end
    display(p)
end