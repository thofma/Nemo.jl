###############################################################################
#
#   acb_calc.jl : Calculus with complex valued functions
#
#   Copyright (C) 2018 Marek Kaluba, Tommy Hofmann
#
###############################################################################

const ARB_CALC_SUCCESS = UInt(0)
const ARB_CALC_NO_CONVERGENCE = UInt(2)

struct FunctionWrapper end

(::FunctionWrapper)(f, args...) = f(args...)

#struct CallWrapper end
#
#(::CallWrapper{Ret})(f, args...)::acb = f(args...)
#
#mutable struct Wrapper
#  ptr::Ptr{Void}
#  objptr::Ptr{Void}
#  obj
#  objT
#
#  function (::Type{Wrapper}){objT}(obj::objT)
#    objref = Base.cconvert(Ref{objT}, obj)
#    new(cfunction(CallWrapper{Ret}(), map_rettype(Ret), get_cfunc_argtype(objT, Args)), Base.unsafe_convert(Ref{objT}, objref), objref, objT)
#  end
#(::Type{FunctionWrapper{Ret,Args}}){Ret,Args}(obj::FunctionWrapper{Ret,Args}) = obj
#                                            end

mutable struct Wrapper end

function integrate(C::AcbField, F, a, b;
                   rel_tol = -1.0,
                   abs_tol = -1.0,
                   deg_limit::Int = 0,
                   eval_limit::Int = 0,
                   depth_limit::Int = 0,
                   use_heap::Int = 0,
                   verbose::Int = 0)

   opts = acb_calc_integrate_opts(deg_limit, eval_limit, depth_limit,
                                  Cint(use_heap), Cint(verbose))

   lower = C(a)
   upper = C(b)

   cgoal = 0

   if rel_tol === -1.0
      cgoal = prec(C)
   else
      t = BigFloat(rel_tol, RoundDown)
      cgoal_clong = Ref{Clong}()
      ccall((:mpfr_get_d_2exp, :libmpfr), Float64, (Ref{Clong}, Ref{BigFloat}, Cint), cgoal_clong, t, Base.MPFR.to_mpfr(RoundDown))
      cgoal = -Int(cgoal_clong[]) + 1
   end

   ctol = mag_struct()
   ccall((:mag_init, :libarb), Void, (Ref{mag_struct},), ctol)

   if abs_tol === -1.0
      ccall((:mag_set_ui_2exp_si, libarb), Void, (Ref{mag_struct}, UInt, Int), ctol, 1, -prec(C))
   else
      t = BigFloat(abs_tol, RoundDown)
      expo = Ref{Clong}()
      d = ccall((:mpfr_get_d_2exp, :libmpfr), Float64, (Ref{Clong}, Ref{BigFloat}, Cint), expo, t, Base.MPFR.to_mpfr(RoundDown))
      ccall((:mag_set_d, libarb), Void, (Ref{mag_struct}, Float64), ctol, d)
      ccall((:mag_mul_2exp_si, libarb), Void, (Ref{mag_struct}, Ref{mag_struct}, Int), ctol, ctol, Int(expo[]))
   end

   res = C()

   #if Base.method_exists(F, (acb, )) && !Base.method_exists(F, (acb, Bool))
   #  G(z, b) = F(z)
   #else
   #  G = F
   #end

   g = function(z, x, para, order, prec)
      xx = unsafe_load(x)
      xx.parent = AcbField(prec)
      w = F(xx, order == 1)
      ccall((:acb_set, :libarb), Ptr{Void}, (Ptr{acb}, Ref{acb}), z, w)
      return zero(Cint)
   end


   F = FunctionWrapper()

   gptr = cfunction(F, Cint, (Ptr{Void}, Ptr{acb}, Ptr{acb}, Ptr{Void}, Int, Int))

   status = ccall((:acb_calc_integrate, :libarb), UInt,
                  (Ref{acb},                       #res
                   Ptr{Void},                      #func
                   Ptr{Void},                      #params
                   Ref{acb},                       #a
                   Ref{acb},                       #b
                   Int,                            #rel_goal
                   Ref{mag_struct},                #abs_tol
                   Ref{acb_calc_integrate_opts},   #opts
                   Int),
   res, gptr, C_NULL, lower, upper, cgoal, ctol, opts, prec(C))

   if status == ARB_CALC_SUCCESS
      nothing
   elseif status == ARB_CALC_NO_CONVERGENCE
      warn("Integration did converge")
   end
   return res
end
