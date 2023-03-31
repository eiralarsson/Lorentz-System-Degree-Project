
function PlotNucleusPlot(p,solutions)
    for X = solutions
        scatter!(p, X[1,[1]], X[2,[1]], X[3,[1]], color=:blue)
        scatter!(p, X[1,[N]], X[2,[N]], X[3,[N]], color=:red)
        plot!(p, X[1,[1,N]], X[2,[1,N]], X[3,[1,N]], color=:black)
    end
    return p
end

function PlotTrajectories(p,solutions; label = "")
    for X = solutions
        plot!(p, X[1,:], X[2,:], X[3,:], label=label) 
    end
    return p
end

function PlotTrajectoriesInterpolated(p,solutions)
    for X = solutions
        plot!(p, X, vars=(1,2,3))
    end
    return p
end

function PlotTrajectoriesOneVariable(p,solutions, index)
    for X = solutions
        plot!(p, X[index,:])
    end
    return p
end

function PlotPipPlot(M,Δt;title="Correlation plot",clim=(-1,1),cmap=:grayC, yaxis=false)
    r,c = size(M)
    y = [1:1:r]
    x = [0:Δt:Δt*(c-1)]
    heatmap(x,y,M,xlabel="time t", title=title,clim=clim,cmap=cmap,yaxis=yaxis, dpi=500)
end