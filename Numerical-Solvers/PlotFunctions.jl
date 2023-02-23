
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

function PlotTrajectoriesOneVariable(solutions, index, legend=false)
    p = plot(legend=legend)
    for X = solutions
        plot!(p, X[index,:])
    end
    display(p)
end

function PlotPipPlot(M,Δt;title="Correlation plot",clim=(-1,1),cmap=:grays, yaxis=false)
    r,c = size(M)
    y = [1:1:r]
    x = [0:Δt:Δt*(c-1)]
    heatmap(x,y,M,xlabel="time t", title=title,clim=clim,cmap=cmap,yaxis=yaxis)
end