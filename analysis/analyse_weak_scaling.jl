using LinearAlgebra
using DataFrames
using CSV
using Glob
using Plots

begin
    options=[:dpi=>200,:size=>(400,300),:fontfamilys => "Computer Modern",:markersize=>4.,:markerstrokewidth=>0,:tickfontsize=>10,:framestyle=>:box]
    files=glob("*.csv","project_eigen/data_weak_scaling")
    df=DataFrame()
    for file in files
        data=CSV.read(file,DataFrame)
        append!(df,data)
    end
    rename!(df,[:n,:p,:k,:hT,:qT,:T])
    sort!(df,[:n,:p,:k])
end

# Making the plots

begin
    vp = unique(df.p)
    plot_folder=mkpath("project_eigen/plots_weak_scaling")*"/"
    p1=scatter(;xaxis="\$p\$",yaxis="\$T\$ (s)",legend=:bottomright,options...)
    
    
    vn=unique(df.n)
    vp=unique(df.p)
    T_av=[ (x->sum(x)/length(x))(df[df.n .==n .&& df.k .!=0. ,:T]) for n in vn]
    plot!(p1,vp,T_av;m=:circ,label="\$n=(W\\times p)^{1/3}\$",options...)
    display(p1)
    savefig(p1,plot_folder*"T_vs_p.pdf")
    speedup=T_av[1] .* vp ./T_av
    efficiency=T_av[1] ./T_av
    p1=plot(vp,speedup;m=:circ,xaxis="\$p\$",yaxis="Speedup = \$S(p)=\\frac{p T_1}{T_p}\$",label="\$n=(W\\times p)^{1/3}\$",options...)
    savefig(p1,plot_folder*"Speedup_vs_p.pdf")
    p1=plot(vp,efficiency;m=:circ,xaxis="\$p\$",yaxis="Effciency = \$E(p)=\\frac{T_1}{T_p}\$",label="\$n=(W\\times p)^{1/3}\$",options...)
    savefig(p1,plot_folder*"Efficiency_vs_p.pdf")
    # plot!(p1,x->x^3/1e10,color=:black,label="\$a\\times n^3\$")
    # display(p1)
    # display(p4)
    
end

# begin
#     gdf2=groupby(df,[:n])
#     p1=scatter(;xaxis="\$p\$",yaxis="\$T\$ (s)",legend=:outerright,yticks=10. .^ (-6:4),options...)
#     p2=scatter(;xaxis="\$p\$",yaxis="Speedup",legend=:outerright,options...)
#     p3=scatter(;xaxis="\$p\$",yaxis="Efficiency",legend=:outerright,options...)
    
#     vn = unique(df.n)
#     for (in,n) in enumerate(vn)
#         vp=unique(gdf2[(n=n,)].p)
#         @show n
#         @show vp
#         T_av=[ (x->sum(x)/length(x))(gdf2[(n=n,)][gdf2[(n=n,)].p .==p .&& gdf2[(n=n,)].k .!=0. ,:T]) for p in vp]
#         @show T_av
#         plot!(p1,vp,T_av;m=:circ,label="\$n=$n\$",xaxis=:log10,yaxis=:log10,color=in,options...)
#         plot!(x->vp[1]*T_av[1]/x;label="",color=in,linestyle=:dash,options...)
#         speedup=T_av[1].*vp[1] ./T_av
#         plot!(p2,vp,speedup;m=:circ,label="\$n=$n\$",options...)
        
#         eff = speedup ./ vp
#         plot!(p3,vp,eff;m=:circ,label="\$n=$n\$",options...)
        

#     end
#     plot!(p1,x->1/x*1e-2,color=:black,label="\$1/p\$")
#     plot!([36,36],[1e-3,1e3],color=:black,linestyle=:dash,label="\$p=36\$")
#     display(p1)
#     display(p2)
#     display(p3)
#     display(p4)
#     savefig(p1,plot_folder*"T_vs_p.pdf")
#     savefig(p2,plot_folder*"speedup.pdf")
#     savefig(p3,plot_folder*"efficiency.pdf")
#     savefig(p4,plot_folder*"T_vs_work.pdf")
# end

# begin
#     df1=gdf[(p=1,)]
#     vn=unique(df1.n)
#     T_av=[ (x->sum(x)/length(x))(df1[df1.n .==n .&& df1.k .!=0. ,:T]) for n in vn]
#     scatter(vn,T_av;label="",xaxis="\$n\$",yaxis="\$T\$",options...)
#     plot!(x->T_av[end]/vn[end]^3*x^3)
# end