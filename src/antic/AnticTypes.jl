###############################################################################
#
#   AnticTypes.jl : Antic types
#
###############################################################################

###############################################################################
#
#   AnticNumberField / nf_elem
#
###############################################################################

const AnticNumberFieldID = Dict{Tuple{FmpqPolyRing, fmpq_poly, Symbol}, Field}()

mutable struct AnticNumberField <: Field
   pol_coeffs::Ptr{Void}
   pol_den::Int
   pol_alloc::Int
   pol_length::Int
   pinv_dinv::Ptr{Void}
   pinv_n::Int
   pinv_norm::Int
   powers::Ptr{Void}
   powers_len::Int
   traces_coeffs::Ptr{Void}
   traces_den::Int
   traces_alloc::Int
   traces_length::Int
   flag::UInt
   pol::fmpq_poly
   S::Symbol
   auxilliary_data::Array{Any, 1}

   function AnticNumberField(pol::fmpq_poly, s::Symbol, cached::Bool = false)
      if !cached
         nf = new()
         nf.pol = pol
         ccall((:nf_init, :libantic), Void, 
            (Ref{AnticNumberField}, Ref{fmpq_poly}), nf, pol)
         finalizer(nf, _AnticNumberField_clear_fn)
         nf.S = s
         nf.auxilliary_data = Array{Any}(5)
         return nf
      else
         if haskey(AnticNumberFieldID, (parent(pol), pol, s))
            return AnticNumberFieldID[parent(pol), pol, s]::AnticNumberField
         else
            nf = new()
            nf.pol = pol
            ccall((:nf_init, :libantic), Void, 
               (Ref{AnticNumberField}, Ref{fmpq_poly}), nf, pol)
            finalizer(nf, _AnticNumberField_clear_fn)
            nf.S = s
            nf.auxilliary_data = Array{Any}(5)
            if cached
               AnticNumberFieldID[parent(pol), pol, s] = nf
            end
            return nf
         end
      end
   end
end

function _AnticNumberField_clear_fn(a::AnticNumberField)
   ccall((:nf_clear, :libantic), Void, (Ref{AnticNumberField},), a)
end

mutable struct nf_elem <: FieldElem
   elem_coeffs::Ptr{Void}
   elem_den::Int
   elem_alloc::Int
   elem_length::Int
   parent::AnticNumberField

   function nf_elem(p::AnticNumberField)
      r = new()
      ccall((:nf_elem_init, :libantic), Void, 
            (Ref{nf_elem}, Ref{AnticNumberField}), r, p)
      r.parent = p
      finalizer(r, _nf_elem_clear_fn)
      return r
   end

   function nf_elem(p::AnticNumberField, a::nf_elem)
      r = new()
      ccall((:nf_elem_init, :libantic), Void, 
            (Ref{nf_elem}, Ref{AnticNumberField}), r, p)
      ccall((:nf_elem_set, :libantic), Void,
            (Ref{nf_elem}, Ref{nf_elem}, Ref{AnticNumberField}), r, a, p)
      r.parent = p
      finalizer(r, _nf_elem_clear_fn)
      return r
   end
end

function _nf_elem_clear_fn(a::nf_elem)
   ccall((:nf_elem_clear, :libantic), Void, 
         (Ref{nf_elem}, Ref{AnticNumberField}), a, a.parent)
end
