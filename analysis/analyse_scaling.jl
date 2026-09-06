using LinearAlgebra
using DataFrames
using CSV
using Glob
using Plots
using LsqFit


begin
    options=[:dpi=>200,:size=>(400,300),:fontfamilys => "Computer Modern",:markersize=>4.,:markerstrokewidth=>0,:tickfontsize=>10,:framestyle=>:box]
    files=glob("*.csv","project_eigen/data")
    df=DataFrame()
    for file in files
        data=CSV.read(file,DataFrame)
        append!(df,data)
    end
    rename!(df,[:n,:p,:k,:hT,:qT,:T])
    sort!(df,[:n,:p,:k])
    gdf=groupby(df,[:p])
end

# Making the plots

begin
    vp = unique(df.p)
    plot_folder=mkpath("project_eigen/plots_scaling")*"/"
    p1=scatter(;xaxis="\$n\$",yaxis="\$T\$ (s)",legend=:bottomright,options...)
    p4=scatter(;xaxis="Work per proc \$\\approx\$ \$n^3/p\$",yaxis="\$T\$ (s)",legend=:outerright,options...)
    p5=scatter(;xaxis="Work per proc \$\\approx\$ \$n^3/p\$",yaxis="Speedup",legend=:outerright,options...)
    for (ip,p) in enumerate(vp)
        vn=unique(gdf[(p=p,)].n)
        T_av=[ (x->sum(x)/length(x))(gdf[(p=p,)][gdf[(p=p,)].n .==n .&& gdf[(p=p,)].k .!=0. ,:T]) for n in vn]
        err(x)= (x.- sum(x)/length(x)).^2 |> y->sqrt(sum(y)/length(y))
        # T_err=[err(gdf[(p=p,)][gdf[(p=p,)].n .==n .&& gdf[(p=p,)].k .!=0. ,:T]) for n in vn]
        plot!(p1,vn,T_av;m=:circ,label="\$P=$p\$",xaxis=:log10,yaxis=:log10,color=ip,options...)
        @show T_av
        @show p
        @show vn
        p2=scatter(vn,T_av,axis="\$n\$",yaxis="\$T\$ (s)")
        savefig(p2,plot_folder*"T_vs_n_p=$p.pdf")
        work=vn.^3
        plot!(p4,work /p,T_av;m=:circ,label="\$p=$p\$",xaxis=:log10,yaxis=:log10,options...)
        # speedup=T_av .* p ./T_av
    end
    plot!(p1,x->x^3/1e10,color=:black,label="\$a\\times n^3\$")
    display(p1)
    display(p4)
    savefig(p1,plot_folder*"T_vs_n.pdf")
    savefig(p4,plot_folder*"work.pdf")
end

begin
    gdf2=groupby(df,[:n])
    p1=scatter(;xaxis="\$p\$",yaxis="\$T\$ (s)",legend=:outerright,yticks=10. .^ (-6:4),options...)
    p2=scatter(;xaxis="\$p\$",yaxis="Speedup=\$S(p)=\\frac{p_{0}T_{p_0}}{T_p}\$",legend=:outerright,options...)
    p3=scatter(;xaxis="\$p\$",yaxis="Efficiency",legend=:outerright,options...)
    
    vn = unique(df.n)
    for (in,n) in enumerate(vn[1:end-4])
        vp=unique(gdf2[(n=n,)].p)
        @show n
        @show vp
        T_av=[ (x->sum(x)/length(x))(gdf2[(n=n,)][gdf2[(n=n,)].p .==p .&& gdf2[(n=n,)].k .!=0. ,:T]) for p in vp]
        @show T_av
        plot!(p1,vp,T_av;m=:circ,label="\$n=$n\$",xaxis=:log10,yaxis=:log10,color=in,options...)
        plot!(x->vp[1]*T_av[1]/x;label="",color=in,linestyle=:dash,options...)
        speedup=T_av[1].*vp[1] ./T_av
        plot!(p2,vp,speedup;color=in,m=:circ,label="\$n=$n\$",options...)
        
        eff = speedup ./ vp
        plot!(p3,vp,eff;m=:circ,label="\$n=$n\$",options...)
        @. amdahl(x,param) = x/(1+(param[1]*x/n)*(1-1/x))
        amdahl2(x,param) = x/(1+(param[1]*x/n)*(1-1/x))
        fit=curve_fit(amdahl,1. .* vp,1. .* speedup,[0.1])
        @show coef(fit)[1]
        plot!(p2,x->amdahl2(x,coef(fit)[1]),color=in,label="\$\\mu=$(round(coef(fit)[1],digits=2))\$",style=:dash)
    end
    plot!(p1,x->1/x*1e-2,color=:black,label="\$1/p\$")
    plot!([36,36],[1e-3,1e3],color=:black,linestyle=:dash,label="\$p=36\$")
    # display(p1)
    display(p2)
    # display(p3)
    # display(p4)
    savefig(p1,plot_folder*"T_vs_p.pdf")
    savefig(p2,plot_folder*"speedup.pdf")
    savefig(p3,plot_folder*"efficiency.pdf")
    savefig(p4,plot_folder*"T_vs_work.pdf")
end

begin
    df1=gdf[(p=1,)]
    vn=unique(df1.n)
    T_av=[ (x->sum(x)/length(x))(df1[df1.n .==n .&& df1.k .!=0. ,:T]) for n in vn]
    scatter(vn,T_av;label="",xaxis="\$n\$",yaxis="\$T\$",options...)
    plot!(x->T_av[end]/vn[end]^3*x^3)
end