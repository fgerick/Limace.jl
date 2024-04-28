module Setup

using .Limace
using ..Bases

export MHDProblem

Base.@kwdef struct MHDProblem{Tu,Tb, NT}
    N::Int = 10
    u::Basis{Tu} = Inviscid(N)
    b::Basis{Tb} = Insulating(N)
    B0 = BasisElement(Basis{Insulating}, Poloidal, (1,0,1), 1.0)
    U0 = BasisElement(Basis{Inviscid}, Toroidal, (1,0,1), 0.0)
    params::NamedTuple{NT} = (Ω=2.0, ν=0.0, η = 0.0)
    threading::Bool = false
end

Base.@kwdef struct HDProblem{Tu,Tb, NT}
    N::Int = 10
    u::Basis{Tu} = Inviscid(N)
    U0 = BasisElement(Basis{Inviscid}, Toroidal, (1,0,1), 0.0)
    params::NamedTuple{NT} = (Ω=2.0, ν=0.0)
    threading::Bool = false
end




end #module