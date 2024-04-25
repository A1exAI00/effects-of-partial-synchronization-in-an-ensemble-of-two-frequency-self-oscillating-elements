# import Pkg; Pkg.add("LSODA")
# import Pkg; Pkg.add("DiffEqBase")
# import Pkg; Pkg.add("Plots")
# import Pkg; Pkg.add("PyPlot")
# import Pkg; Pkg.add("Sundials")
# import Pkg; Pkg.add("SciPyDiffEq")
# import Pkg; Pkg.add("OrdinaryDiffEq")
# import Pkg; Pkg.add("ODEInterfaceDiffEq")
# import Pkg; Pkg.add("DifferentialEquations")
# import Pkg; Pkg.add("Gtk")
# import Pkg; Pkg.add("JLD2")
# import Pkg; Pkg.add("PlotlyJS")
# import Pkg; Pkg.add("LaTeXStrings")




using LinearAlgebra
using DifferentialEquations
using ODEInterfaceDiffEq
using LSODA, DiffEqBase
using SciPyDiffEq
using OrdinaryDiffEq
using Sundials
# using Plots
using Gtk
using JLD2
using Dates
using LaTeXStrings
using Printf
# using PlotlyJS
using PyPlot

function fhn!(xdot, x, p, t)
    for k in 1:M:M*K
        if k==1
            xdot[k]   = -x[k]*(x[k]^2-u1^2)*(x[k]^2-u2^2)*(x[k]^2-u3^2)-x[k+1]+d*(x[M*K-M+1]-2*x[k]+x[k+2])
            xdot[k+1] = eps1*(x[k]-E0)
        end

        if k==(M*K-M+1)
            xdot[k]   = -x[k]*(x[k]^2-u1^2)*(x[k]^2-u2^2)*(x[k]^2-u3^2)-x[k+1]+d*(x[k-2]-2*x[k]+x[1])
            xdot[k+1] = eps1*(x[k]-E0)
        end

        if 1<k<(M*K-M+1)
            xdot[k]   = -x[k]*(x[k]^2-u1^2)*(x[k]^2-u2^2)*(x[k]^2-u3^2)-x[k+1]+d*(x[k-2]-2*x[k]+x[k+2])
            xdot[k+1] = eps1*(x[k]-E0)
        end            
    end
nothing  
end


K,M,eps1,d,E0,u1,u2,u3,x0=load_object("d:\\Shapin_new\\Papers\\2023_china\\Александр Акинин\\7-power\\ini\\01 04 24 09 04 19_par_ini.jld2")


d=.006
u1=0.32; u2=0.79; u3=1.166
E0=0
eps1=0.001

K=200 # количество элементов
M=2 # размерность элемента
t_end = 1000 # период непрерывного интегрирования системы
cycle_num=1 # количество периодов непрерывного интегрирования системы
pont_number=t_end # время через которое обновляется мгновенный снимок фаз
t_show_phi= t_end*(cycle_num-1) # время после которого записывается пространственно-временная диаграмма для фаз
phi_time=.1 # интервал времени через который записывается пространственно-временная диаграмма для фаз

nado=1 # если не равно нулю, то рисуются все переменные от времени и парциальные проекции элементов

Ft = t -> -t*(t^2-u1^2)*(t^2-u2^2)*(t^2-u3^2)



# ss=1
# kk=1
# mu=0.01
# sr=range(0.21,0.38,Int(floor(K/8+2)))
# sn=range(-0.21,-0.38,Int(floor(K/8+2)))
# br=range(1.06,1.208,Int(floor(K/8+2)))
# bn=range(-1.06,-1.208,Int(floor(K/8+2)))
# x0 = zeros(K*M) # начальные условия
# for s in 1:M:M*K
#     global ss,sr,sn,br,bn,mu,kk
  
#     if s<M*K/4
#         x0[s]=-1.25
#         x0[s+1]=0.13-0.52*(M*K-2*s)/(M*K)
        
#     end

#     if (s>=M*K/4)*(s<Int(M*K/2-2))
#         x0[s]=1.25#+0.5*rand()
#         x0[s+1]=0.39-0.52*(M*K-2*s)/(M*K)
#     end

#     if s>=Int(M*K/2-2)
#         if mod(kk,4)==0
#         ss+=1
#         end
#         kk+=1
        

#         if mod(s,8)==1
#             x0[s]=sr[ss]+mu
#             x0[s+1]=Ft(sr[ss])
#         end
        
#         if mod(s,8)==3
#             x0[s]=bn[ss]-mu
#             x0[s+1]=Ft(bn[ss])
#         end
        
#         if mod(s,8)==5
#             x0[s]=sn[ss]-mu
#             x0[s+1]=Ft(sn[ss])
#         end
      
#         if mod(s,8)==7
#             x0[s]=br[ss]+mu
#             x0[s+1]=Ft(br[ss])
#         end
#     end
# end   


# figure(8)
# clf()
# for i in K+1:M:M*K
#     plt.axis([-1.3,1.3, -0.16, 0.16])
#     scatter(x0[i], x0[i+1])
# end

tspan = (0., t_end)
T0=zeros(M*K) # время предыдущего пересечения порога
T1=zeros(M*K) # время последнего пересечения порога
w_aver=zeros(K)  # средние частоты 
W_aver=zeros(M*K) #  сумма средних частот
T_aver=zeros(M*K) # сумма периодов [k] и их количество [k+1]
Amp_aver=zeros(K) # средние амплитуды
A_aver=zeros(M*K) # cумма амплитутд [k] и их количество [k+1]
Tphi=0 # моменты времени запси мгновенных фаз
Phiplot=zeros(K) # мгновенные фазы всех элементов, записанные в моменты времени Tphi
Prod=0 # произведение моментов времён первого пересечения порогв всеми элементми (если Prod==0 значит не все пересекли порог)
Phi=zeros(K) # мгновенные фазы

for number in 1:cycle_num
global  T0,T1,Prod,Phi,sts,x0,tspan,Tt,y,prob,w_aver,T_aver,Amp_aver,A_aver,phi_time,Tphi,Phiplot,t_show_phi,W_aver
prob = ODEProblem(fhn!,x0,tspan) # ,reltol= 1.e-9, abstol = 1.e-10,adaptive=false,dt=0.001)
sts=@timed y = solve(prob,RK4())
# sts=@timed y = solve(prob,lsoda())

# sts=@timed y = solve(prob,FBDF())
# sts=@timed y = solve(prob,rodas())
# sts=@timed y = solve(prob,RadauIIA5())  
# sts=@timed y = solve(prob,CVODE_BDF())  
# sts=@timed y = solve(prob,Rodas5P()) 
print(sts.time,'\n')
print(number,'\n')

Tt=y.t[end]
tspan = (number*t_end,(number+1)*t_end)
x0=y[:,end]




    if number == cycle_num
        figure(1)
        clf()
        for i in 1:M:M*K
            subplot(3,1,1)
            plt.axis([-1.3,1.3, -0.16, 0.16])
            plot(y[i,:], y[i+1,:])
            subplot(3,1,2)
            plot(y.t[:], y[i,:])
            subplot(3,1,3)
            plot(y.t[:], y[i+1,:])
        end
        subplot(3,1,1)
        zz=-1.22:.01:1.22
        rr=-.15:.01:0.15
        plot(zz, Ft.(zz))
        plot(ones(size(rr))*E0, rr)
    end

    for i in 1:1:size(y.t,1)-1    
            for k in 1:M:M*K   
            
                if y[k,i]*(y[k,i+1])*(y[k+1,i]>0) <0
                    if T1[k]==0
                        T1[k]=y.t[i]
                        # T1[k+1]=i
                        else
                        T0[k]=T1[k]
                        T1[k]=y.t[i]

                        T_aver[k]+=T1[k]-T0[k]
                        T_aver[k+1]+=1

                        W_aver[k]+=1/(T1[k]-T0[k])
                        W_aver[k+1]+=1

                        A_aver[k]+=y[k+1,i]
                        A_aver[k+1]+=1
                    end

                    Prod=T0[1]           
                    for k in M+1:M:M*K
                    rod=Prod*T0[k]
                    end
                end
            end


            if Prod>0
                if mod(y.t[i],pont_number)>mod(y.t[i+1],pont_number)
                    Phi[1]=0
                    w_aver[1]=W_aver[1]/W_aver[2]
                    Amp_aver[1]=A_aver[1]/A_aver[2]
                    for k in M+1:M:M*K
                        Phi[Int((k+1)/M)]=mod((2*pi)*(((y.t[i]-T1[k])/(T1[k]-T0[k]))-((y.t[i]-T1[1])/(T1[1]-T0[1])))-pi,2*pi)-pi
                        w_aver[Int((k+1)/M)]=W_aver[k]/W_aver[k+1]
                        Amp_aver[Int((k+1)/M)]=A_aver[k]/A_aver[k+1]
                    end
                    # if number == cycle_num
                        figure(2)
                        clf()
                        subplot(3,1,1)
                        plt.axis([1,K, -pi, pi])
                        scatter(1:1:K,Phi) #-ones(size(Phi,1))*pi
                        subplot(3,1,2)
                        plt.xlim(1,K)
                        scatter(1:1:K,w_aver)
                        subplot(3,1,3)
                        plt.xlim(1,K)
                        scatter(1:1:K,Amp_aver)
                    # end
                end
        
        
                    if (y.t[i]>t_show_phi)*mod(y.t[i],phi_time)>mod(y.t[i+1],phi_time)
                        Phi[1]=0
                        for k in 1:M:M*K
                            Phi[Int((k+1)/M)]=mod((2*pi)*(((y.t[i]-T1[k])/(T1[k]-T0[k]))-((y.t[i]-T1[1])/(T1[1]-T0[1])))-pi,2*pi)-pi
                        end
                        if Tphi==0
                            Tphi=y.t[i]
                            Phiplot=Phi
                            else
                            Tphi=[Tphi y.t[i]]
                            Phiplot=[Phiplot Phi]
                        end
                    end
        end
    end


end


name=Dates.format(Dates.now(), "dd mm yy HH MM SS")
name=string(" K=",K," ","M=", M," " ,"d=",d," ","eps=", eps1," ","E0=", E0,"_")*name
mkpath("d:\\Shapin_new\\Papers\\2023_china\\Александр Акинин\\7-power\\ini\\"*name*"")
savefig("d:\\Shapin_new\\Papers\\2023_china\\Александр Акинин\\7-power\\ini\\"*name*"\\Phases(j)_.png",dpi=600)
save_object("d:\\Shapin_new\\Papers\\2023_china\\Александр Акинин\\7-power\\ini\\"*name*"\\_par_ini.jld2", (K,M,eps1,d,E0,u1,u2,u3,x0))
save_object("d:\\Shapin_new\\Papers\\2023_china\\Александр Акинин\\7-power\\ini\\"*name*"\\Phi_w_Amp_.jld2", (Phi,w_aver,Amp_aver))
# print("file save time ", sts.time,'\n')

# figure(3)
# clf()
# zu=y[1,:]
# for i in M+1:M:M*K
#     global zu
#     zu=[zu y[i,:]]
# end
# contourf(y.t,1:1:K, zu',  levels=10,cmap=get_cmap("winter"))
# plt.colorbar()
# savefig("d:\\Shapin_new\\Papers\\2023_china\\Александр Акинин\\7-power\\ini\\"*name*"\\ ST_u_j_.png",dpi=600)


# figure(4)
# clf()
# contourf(Tphi,1:1:K, Phiplot,  levels=10,cmap=get_cmap("hot"))
# plt.colorbar()
# savefig("d:\\Shapin_new\\Papers\\2023_china\\Александр Акинин\\7-power\\ini\\"*name*"\\ ST_Phases_.png",dpi=600)

# figure(5)
# clf()
# for i in 1:1:Int(K/2)
#     plot(Tphi',Phiplot[i,:])
# end
# savefig("d:\\Shapin_new\\Papers\\2023_china\\Александр Акинин\\7-power\\ini\\"*name*"\\ Phases(t)_sinhr_.png",dpi=600)

# figure(6)
# clf()
# for i in Int(K/2):1:K
#     plot(Tphi',Phiplot[i,:])
# end