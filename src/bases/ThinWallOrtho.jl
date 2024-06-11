module ThinWallBC

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..UnconstrainedBasis
using ..InviscidBasis
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t

export ThinWall

struct ThinWall; end

function ThinWall(N; σw=1.0, σf = 1.0, h = 1.0, kwargs...)
    params=Dict(:σw => σw, :σf => σf, :h => h)
    return Basis{ThinWall}(;N, V=Sphere(), BC=NoBC(), params,  kwargs...)
end

function s(b::Basis{ThinWall}, l,m,n,r) 
    h, σf, σw = b.params[:h], b.params[:σf], b.params[:σw]
    λ1 = 1.0
    λ2 = 1.0

    c[1] = 8*n*(l+2*n)*(n+3+l)*(2835-7560*n^2*λ1^2-32640*n^6*λ2^2+9072*n^4*λ1^2+44016*n^4*λ2^2-24080*n^2*λ2^2+270*λ1*λ2*l^3+4644*λ1*λ2*l^2+702*λ1*λ2*l+1155*λ2^2*l^4-3672*λ1*λ2+8960*λ2^2*n^8+6048*λ2*n^4-7560*λ2*n^2+7560*λ1*n^2-1890*λ1+1512*λ2+22680*n^2*λ1*λ2+17280*n^6*λ1*λ2-1890*λ1*l+2646*λ2*l^2+4158*λ2*l-1701*λ1^2*l^2-378*λ1^2*l+750*λ2^2*l^3-1704*λ2^2*l-5853*λ2^2*l^2-26880*n^5*λ2^2*l^2-97920*n^5*λ2^2*l+53760*n^6*λ2^2*l^2-8960*n^6*λ2^2*l+35840*n^7*λ2^2*l+44480*n^3*λ2^2*l^2-96800*n^4*λ2^2*l^2+22240*n^4*λ2^2*l+560*λ2^2*l^4*n+35840*λ2^2*l^3*n^5+88032*n^3*λ2^2*l+8960*λ2^2*l^4*n^4-26880*λ2^2*l^3*n^4-30400*λ2^2*l^3*n^3+22800*λ2^2*l^3*n^2-8960*λ2^2*l^4*n^3+1120*λ2^2*l^4*n^2+3560*λ2^2*l^3*n+15984*n*λ1*λ2*l^2-72576*n^3*λ1*λ2*l-25920*n^3*λ1*λ2*l^2+51840*n^4*λ1*λ2*l^2+17280*n^3*λ1*λ2*l^3-12960*n^4*λ1*λ2*l+51840*n^5*λ1*λ2*l+1323*λ1^2+15984*n^2*λ1*λ2*l-35208*n^2*λ1*λ2*l^2+22680*n*λ1*λ2*l-12960*λ1*λ2*l^3*n^2+3744*λ2^2-36288*n^4*λ1*λ2-13448*n^2*λ2^2*l-24080*λ2^2*l*n-13448*λ2^2*l^2*n+47576*λ2^2*l^2*n^2+6048*λ2*n^2*l^2-3024*λ2*n^2*l+12096*λ2*n^3*l-3024*λ2*n*l^2-7560*λ2*n*l+7560*λ1*n*l+1080*λ1*λ2*l^3*n-4536*n*λ1^2*l^2-7560*n*λ1^2*l+9072*n^2*λ1^2*l^2-4536*n^2*λ1^2*l+18144*n^3*λ1^2*l)*(l-1+2*n)*(n+2+l)
    c[2] = -12*(l-1+2*n)*(n+2+l)*(2+l+2*n)*(2835-20034*n^2*λ1^2+18144*n^6*λ1^2-32928*n^6*λ2^2+3213*n*λ1^2+30240*n^4*λ1^2+39312*n^5*λ1^2-22680*n^3*λ1^2-2800*n^4*λ2^2+3744*n*λ2^2-27440*n^3*λ2^2+2448*n^2*λ2^2+74256*n^5*λ2^2-113280*n^7*λ2^2+4725*n+945*l+630*λ1*λ2*l^4+3150*λ1*λ2*l^3+4410*λ1*λ2*l^2+1890*λ1*λ2*l+13545*λ2^2*l^4+17920*λ2^2*n^10+62720*λ2^2*n^9+15360*λ2^2*n^8+1512*λ2*n+15120*λ2*n^4+17640*λ2*n^3+40824*λ2*n^2+26208*λ2*n^5+22680*λ1*n^3+12096*λ2*n^6+15120*λ1*n^4+1890*λ1*n+18900*λ1*n^2+5670*n^2+4410*λ2*l^3+5670*λ1+216*n^2*λ1*λ2+48384*n^6*λ1*λ2+3780*n*λ1*l^2+18144*n^3*λ1^2*l^3-3780*λ1*l+13860*λ2*l^2+1890*λ2*l+3465*λ2^2*l^5-10395*λ1^2*l^2-2835*λ1^2*l^3-4725*λ1^2*l-1890*λ1*l^2+9135*λ2^2*l^3-945*λ2^2*l^2+8960*λ2^2*l^5*n^4+179200*λ2^2*l^3*n^7+304640*λ2^2*l^3*n^6-561440*n^5*λ2^2*l^2-29504*n^5*λ2^2*l-120640*n^6*λ2^2*l^2-416320*n^6*λ2^2*l-1280*n^7*λ2^2*l+223944*n^3*λ2^2*l^2+80496*n^4*λ2^2*l^2+207296*n^4*λ2^2*l+412160*n^7*λ2^2*l^2+28905*λ2^2*l^4*n-200640*λ2^2*l^3*n^5+34560*λ2*n^8*λ1+89600*λ2^2*n^9*l+17920*λ2^2*l^5*n^5+89600*λ2^2*l^4*n^6+2792*n^3*λ2^2*l+98560*λ2^2*l^4*n^5+54432*n^5*λ1^2*l+54432*n^4*λ1^2*l^2-24640*λ2^2*l^5*n^3-121280*λ2^2*l^4*n^4-322560*λ2^2*l^3*n^4+121712*λ2^2*l^3*n^3+111724*λ2^2*l^3*n^2+4480*λ2^2*l^5*n^2-59680*λ2^2*l^4*n^3+3990*λ2^2*l^5*n+48630*λ2^2*l^4*n^2+49104*λ2^2*l^3*n-76608*n^5*λ1*λ2+34560*n^4*λ1*λ2*l^4+26208*n*λ1*λ2*l^2-4176*n^3*λ1*λ2*l-228312*n^3*λ1*λ2*l^2-29808*n^4*λ1*λ2*l^2-81216*n^3*λ1*λ2*l^3+141120*n^4*λ1*λ2*l^3-232128*n^4*λ1*λ2*l+5670*n*l+308160*n^6*λ1*λ2*l+138240*n^7*λ1*λ2*l+336960*n^5*λ1*λ2*l^2+71712*n^5*λ1*λ2*l+26460*λ1*n^2*l+2835*λ1^2+46764*n^2*λ1*λ2*l+32832*n^2*λ1*λ2*l^2+207360*n^6*λ1*λ2*l^2+138240*n^5*λ1*λ2*l^3+2178*n*λ1*λ2*l+14400*λ1*λ2*l^4*n^3-69732*λ1*λ2*l^3*n^2+3060*λ1*λ2*l^4*n-15120*n^4*λ1*λ2+232*n^2*λ2^2*l+744*λ2^2*l*n+27720*n^3*λ1*λ2+84672*n^4*λ1^2*l+21189*λ2^2*l^2*n+47646*λ2^2*l^2*n^2+2268*λ2*n^2*l^2+34272*λ2*n^3*l^2+30240*λ1*n^3*l+36288*λ2*n^5*l+36288*λ2*n^4*l^2+23436*λ2*n^2*l+17136*λ2*n^3*l+17262*λ2*n*l^2+46242*λ2*n*l+13230*λ1*n*l+4032*n^2*λ2*l^3+12096*n^3*λ2*l^3+252*n*λ2*l^3+56448*λ2*n^4*l-28080*λ1*λ2*l^4*n^2+23418*λ1*λ2*l^3*n+15120*n^2*λ1*l^2+179200*λ2^2*n^8*l^2+259840*λ2^2*n^8*l+97920*n^7*λ1*λ2-10962*n*λ1^2*l^3-38367*n*λ1^2*l^2-24192*n*λ1^2*l+6048*n^2*λ1^2*l^3-378*n^2*λ1^2*l^2-53676*n^2*λ1^2*l+51408*n^3*λ1^2*l^2+40824*n^3*λ1^2*l-3672*n*λ1*λ2);  
    c[3] = 6*(l+2*n+3)*(2*l+2*n-1)*(3780+108486*n^2*λ1^2+18144*n^6*λ1^2-315168*n^6*λ2^2+58023*n*λ1^2+105840*n^4*λ1^2+69552*n^5*λ1^2+113400*n^3*λ1^2+32480*n^4*λ2^2+33984*n*λ2^2+322000*n^3*λ2^2+188928*n^2*λ2^2-419664*n^5*λ2^2+128640*n^7*λ2^2+6615*n+1890*l-7560*λ1*λ2*l^4-18900*λ1*λ2*l^3+7560*λ1*λ2*l^2+41580*λ1*λ2*l-18900*λ2^2*l^4+22680*λ1*λ2+17920*λ2^2*n^10+116480*λ2^2*n^9+257280*λ2^2*n^8+29232*λ2*n+65520*λ2*n^4+22680*λ2*n^3-2016*λ2*n^2+46368*λ2*n^5+37800*λ1*n^3+12096*λ2*n^6+15120*λ1*n^4+28350*λ1*n+41580*λ1*n^2+5670*n^2-3780*λ2*l^3+15120*λ1+22680*λ2+229536*n^2*λ1*λ2+330624*n^6*λ1*λ2+11340*n*λ1*l^2+18144*n^3*λ1^2*l^3+7560*λ1*l-7560*λ2*l^2+11340*λ2*l-5670*λ2^2*l^5-1890*λ1^2*l^3+13230*λ1^2*l-5670*λ2^2*l^3+22680*λ2^2*l+30240*λ2^2*l^2+35840*λ2^2*l^5*n^4+179200*λ2^2*l^3*n^7+680960*λ2^2*l^3*n^6-30560*n^5*λ2^2*l^2-1071104*n^5*λ2^2*l+1196480*n^6*λ2^2*l^2+242240*n^6*λ2^2*l+912640*n^7*λ2^2*l-1071336*n^3*λ2^2*l^2-1431504*n^4*λ2^2*l^2-1044304*n^4*λ2^2*l+842240*n^7*λ2^2*l^2-99405*λ2^2*l^4*n+686400*λ2^2*l^3*n^5+34560*λ2*n^8*λ1+89600*λ2^2*n^9*l+17920*λ2^2*l^5*n^5+89600*λ2^2*l^4*n^6-5608*n^3*λ2^2*l+259840*λ2^2*l^4*n^5+54432*n^5*λ1^2*l+54432*n^4*λ1^2*l^2+2240*λ2^2*l^5*n^3+147520*λ2^2*l^4*n^4-305760*λ2^2*l^3*n^4-903088*λ2^2*l^3*n^3-520796*λ2^2*l^3*n^2-32480*λ2^2*l^5*n^2-194080*λ2^2*l^4*n^3-17850*λ2^2*l^5*n-245370*λ2^2*l^4*n^2-129396*λ2^2*l^3*n+245952*n^5*λ1*λ2+34560*n^4*λ1*λ2*l^4-36792*n*λ1*λ2*l^2+26064*n^3*λ1*λ2*l+119448*n^3*λ1*λ2*l^2+726192*n^4*λ1*λ2*l^2+201024*n^3*λ1*λ2*l^3+342720*n^4*λ1*λ2*l^3+423072*n^4*λ1*λ2*l+5670*n*l+590400*n^6*λ1*λ2*l+138240*n^7*λ1*λ2*l+699840*n^5*λ1*λ2*l^2+857952*n^5*λ1*λ2*l+49140*λ1*n^2*l+11340*λ1^2+129924*n^2*λ1*λ2*l-156168*n^2*λ1*λ2*l^2+207360*n^6*λ1*λ2*l^2+138240*n^5*λ1*λ2*l^3+157158*n*λ1*λ2*l+54720*λ1*λ2*l^4*n^3-82332*λ1*λ2*l^3*n^2-24660*λ1*λ2*l^4*n+85680*n^4*λ1*λ2+328672*n^2*λ2^2*l+129264*λ2^2*l*n+153720*n^3*λ1*λ2+160272*n^4*λ1^2*l+47439*λ2^2*l^2*n-185874*λ2^2*l^2*n^2+47628*λ2*n^2*l^2+74592*λ2*n^3*l^2+30240*λ1*n^3*l+36288*λ2*n^5*l+36288*λ2*n^4*l^2+15876*λ2*n^2*l+107856*λ2*n^3*l-6678*λ2*n*l^2-10458*λ2*n*l+32130*λ1*n*l+14112*n^2*λ2*l^3+12096*n^3*λ2*l^3+5292*n*λ2*l^3+106848*λ2*n^4*l+2160*λ1*λ2*l^4*n^2-93762*λ1*λ2*l^3*n+15120*n^2*λ1*l^2+179200*λ2^2*n^8*l^2+501760*λ2^2*n^8*l+178560*n^7*λ1*λ2-3402*n*λ1^2*l^3+16443*n*λ1^2*l^2+77868*n*λ1^2*l+21168*n^2*λ1^2*l^3+67662*n^2*λ1^2*l^2+127764*n^2*λ1^2*l+111888*n^3*λ1^2*l^2+176904*n^3*λ1^2*l+124848*n*λ1*λ2)*(l+2*n);   
    c[4] = -(1+2*n)*(2*l+2*n-1)*(2*l-3+2*n)*(l+2*n+3)*(2+l+2*n)*(2835+46872*n^2*λ1^2+218240*n^6*λ2^2+21168*n*λ1^2+9072*n^4*λ1^2+36288*n^3*λ1^2+181616*n^4*λ2^2+3744*n*λ2^2+25024*n^3*λ2^2+1296*n^2*λ2^2+305920*n^5*λ2^2+71680*n^7*λ2^2+5670*λ1*λ2*l^3+11340*λ1*λ2*l^2+5670*λ1*λ2*l+2835*λ2^2*l^4+8960*λ2^2*n^8+9072*λ2*n+6048*λ2*n^4+24192*λ2*n^3+28728*λ2*n^2+15120*λ1*n+7560*λ1*n^2+5670*λ1+64152*n^2*λ1*λ2+17280*n^6*λ1*λ2+5670*λ1*l+5670*λ2*l^2+5670*λ2*l+2835*λ1^2*l^2+5670*λ1^2*l+5670*λ2^2*l^3+2835*λ2^2*l^2+295680*n^5*λ2^2*l^2+600960*n^5*λ2^2*l+53760*n^6*λ2^2*l^2+241920*n^6*λ2^2*l+35840*n^7*λ2^2*l+463680*n^3*λ2^2*l^2+575200*n^4*λ2^2*l^2+652640*n^4*λ2^2*l+11760*λ2^2*l^4*n+35840*λ2^2*l^3*n^5+272992*n^3*λ2^2*l+8960*λ2^2*l^4*n^4+152320*λ2^2*l^3*n^4+220480*λ2^2*l^3*n^3+128720*λ2^2*l^3*n^2+26880*λ2^2*l^4*n^3+28000*λ2^2*l^4*n^2+29640*λ2^2*l^3*n+103680*n^5*λ1*λ2+75168*n*λ1*λ2*l^2+393984*n^3*λ1*λ2*l+181440*n^3*λ1*λ2*l^2+51840*n^4*λ1*λ2*l^2+17280*n^3*λ1*λ2*l^3+246240*n^4*λ1*λ2*l+51840*n^5*λ1*λ2*l+2835*λ1^2+238896*n^2*λ1*λ2*l+198072*n^2*λ1*λ2*l^2+44280*n*λ1*λ2*l+38880*λ1*λ2*l^3*n^2+222912*n^4*λ1*λ2+23128*n^2*λ2^2*l+9600*λ2^2*l*n+200448*n^3*λ1*λ2+16104*λ2^2*l^2*n+137816*λ2^2*l^2*n^2+6048*λ2*n^2*l^2+33264*λ2*n^2*l+12096*λ2*n^3*l+9072*λ2*n*l^2+22680*λ2*n*l+7560*λ1*n*l+27000*λ1*λ2*l^3*n+13608*n*λ1^2*l^2+37800*n*λ1^2*l+9072*n^2*λ1^2*l^2+49896*n^2*λ1^2*l+18144*n^3*λ1^2*l+3888*n*λ1*λ2);  

    return r^l*sum(c[i]*jacobi(n+1-i,3/2,l+1/2,2r^2-1) for i in 1:4)

end    

function t(::Type{Basis{ThinWall}}, V::Volume, l,m,n,r) 
    t(Basis{Unconstrained}, V, l,m,n,r) 
end

@inline _nrange_p(b::Basis{ThinWall},l) = 0:((b.N-l+1)÷2+1)
@inline _nrange_t(b::Basis{ThinWall},l) = 0:((b.N-l)÷2+1)

@inline function bcs_p(b::Basis{ThinWall}) 
    @inline _s = (l,n,r) -> r*s(Basis{ThinWall}, b.V, l, 0, n, r)
    (; r1) = b.V 
    h, σf, σw = b.params[:h], b.params[:σf], b.params[:σw]
    fs = (
          @inline((l,n) -> σw*h/σf*(∂(r->∂(r->_s(l,n,r),r), r1) - l*(l+1)/r1^2*_s(l,n,r1)) + _s(l,n,r1)*l/r1 + ∂(r->_s(l,n,r),r1)*(1 + l*h/r1)), 
          )
    return fs
end

@inline function bcs_t(b::Basis{ThinWall}) 
    @inline _t = (l,n,r) -> r*t(Basis{ThinWall}, b.V, l, 0, n, r)
    (; r1) = b.V 
    h, σf, σw = b.params[:h], b.params[:σf], b.params[:σw]
    fs = (@inline((l,n) -> σw/σf*h*∂(r->_t(l,n,r),r1) + _t(l,n,r1)), )
    return fs
end


lpmax(b::Basis{ThinWall}) = b.N
ltmax(b::Basis{ThinWall}) = b.N


end