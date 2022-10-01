using Plots,Plots.Measures,Printf
default(size=(1200,400),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)

@views function diffusion_1D()
    # physics
    lx   = 20.0
    c   = 1.0
    # numerics
    nx   = 200
    nvis = 2
    # derived numerics
    dx   = lx/nx
    dt   = dx^2/c/2
    nt   = nx^2 ÷ 100
    xc   = LinRange(dx/2,lx-dx/2,nx)
    # array initialisation
    Pr    = @. 0.5cos(9π*xc/lx)+0.5; C_i = copy(Pr)
    Vx   = zeros(Float64, nx-1)
    # time loop
    @gif for it = 1:nt
        Vx          .-= c.*diff(Pr )./dx .*dt #update Vx
        Pr[2:end-1] .-=   dt.*diff(Vx)./dx
        Pr[1] = Pr[2]
        Pr[end] = Pr[end-1]
        plot(xc,[C_i,Pr];xlims=(0,lx), ylims=(-0.1,1.1),
                        xlabel="lx", ylabel="Concentration",
                        title="time = $(round(it*dt,digits=1))")
    end every 10
end

diffusion_1D()
