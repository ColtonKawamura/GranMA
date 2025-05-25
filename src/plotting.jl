using GLMakie
using LaTeXStrings






function     plotRealImagEigenValues(plotData, pressureArray, dampingArray; loglog=false, semilogx=false, semilogy=false, scaling=false)
    

    for pressure in pressureArray
        for damping in dampingArray
            
        pressure



    f = Figure()
    ax = Axis(f[1, 1],
        xlabel = L"Im(\lambda)/\sqrt{\hat{P}}",
        ylabel = L"Re(\lambda)/(\sqrt{\hat{P}\hat{\gamma}})",
        title = L"\text{My LaTeX Title}"
    )

    # Example with LaTeX tick labels
    ax.xticks = ([0, 1, 2, 3], [L"0", L"1", L"2", L"3"])
    ax.yticks = ([0, 0.5, 1.0], [L"0", L"\frac{1}{2}", L"1"])

    scatter!(ax, rand(10), rand(10))
    f
end

