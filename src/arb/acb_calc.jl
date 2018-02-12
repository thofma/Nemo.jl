function acb_calc_func_wrap(res::acb_struct, x::acb, param::Ptr{Void}, order::Int, prec::Int)

    acbF = unsafe_pointer_to_objref(param)::Function
    FF = AcbField(prec)
    z = acbF(FF(x))

    _acb_set(res, acb_struct(z))
    return Cint(0)::Cint
end

acb_calc_func_wrap_c() = cfunction(acb_calc_func_wrap, Cint,
        (Ref{acb_struct}, Ref{acb}, Ptr{Void}, Int, Int))

const ARB_CALC_SUCCESS = UInt(0)
const ARB_CALC_NO_CONVERGENCE = UInt(2)

function _integrate(F, a, b, abs_tol::Int = 53)
    opts = acb_calc_integrate_opts()
    _abs_tol = mag_set_ui_2exp_si(mag_struct(), 1, -abs_tol)
    rel_goal = abs_tol

    res = AcbField(abs_tol)()

    g = function(z, x, para, order, prec)
      xx = unsafe_load(x)
      xx.parent = AcbField(prec)
      w = F(xx)
      ccall((:acb_set, :libarb), Ptr{Void}, (Ptr{acb}, Ref{acb}), z, w)
      return zero(Cint)
    end

    gptr = cfunction(g, Cint, (Ptr{acb}, Ptr{acb}, Ptr{Void}, Int, Int))

    status = ccall((:acb_calc_integrate, :libarb), UInt,
        (Ref{acb},   #res
        Ptr{Void},   #func
        Ptr{Void},   #params
        Ref{acb},    #a
        Ref{acb},    #b
        Int,         #rel_goal
        Ref{mag_struct},    #abs_tol
        Ref{acb_calc_integrate_opts}, #opts
        Int),
        res, gptr, C_NULL, a, b, rel_goal, _abs_tol, opts, abs_tol)

    return res
end

function integrate(F::Function, a::acb, b::acb;
    rel_goal::Int=prec(parent(a)),
    abs_tol::mag_struct=mag_set_ui_2exp_si(mag_struct(), 1, -prec(parent(a))),
    opts::acb_calc_integrate_opts=acb_calc_integrate_opts(),
    precision::Int=prec(parent(a)))

    res = AcbField(precision)()
    ptrF = pointer_from_objref(F)

    status = ccall((:acb_calc_integrate, :libarb), UInt,
        (Ref{acb},   #res
        Ptr{Void},   #func
        Ptr{Void},   #params
        Ref{acb},    #a
        Ref{acb},    #b
        Int,         #rel_goal
        Ref{mag_struct},    #abs_tol
        Ref{acb_calc_integrate_opts}, #opts
        Int),
        res, acb_calc_func_wrap_c(), ptrF, a, b, rel_goal, abs_tol, opts, precision)

    if status == ARB_CALC_SUCCESS
        nothing
    elseif status == ARB_CALC_NO_CONVERGENCE
        warn("Integration did not obtained convergence")
    end
    return res
end
