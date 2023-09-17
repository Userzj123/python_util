module lst_scalar
using LinearAlgebra

export generate_LST

struct chebyshev_basis
    D0 :: Array
    D1 :: Array
    D2 :: Array
    D3 :: Array
    D4 :: Array
end

struct baseflow
    U :: Array
    Up :: Array
    Upp :: Array
    T :: Array
    Tp :: Array
end
struct lst_result
    basis :: chebyshev_basis
    base :: baseflow
    axis :: Vector
    F :: GeneralizedEigen
    A :: Array
    B :: Array
end


    function generate_LST(N::Int, bf::Int)
        """Initialize the LST Analysis

        Args:
            N (int): number of Chebyshev polynomials
            bf (int): set = 1 for Couette, 2 for Poiseuille, 3 for quiescent
            
        Returns:

        """
        # Input parameters 
        R  = 180         # Reynolds number
        kx = 1           # streamwise wavenumber
        kz = 0           # spanwise wavenumber
        Ri = 0.0        # Richardson number
        Pr = 0.71        # Prantl number
        Ra = 8*R*R/Pr*Ri # Rayleigh number
        
        # Set up grid and differentiation matrices
        y_phys = cos.(LinRange(0, pi, N))   # Generate Chebyshev grid for base flow solver
        basis = dmat(N)   # Chebyshev polynomials and derivatives at the Gauss points
        
        
        # Find the base flow
        base = bounded_base(y_phys, N, bf)
        
        # Find eigenvalues of stability operators
        A, B = Operator(kx,kz,R,Pr,Ri, base, basis)
        
        # find eigenvalues
        F = eigen(A, B)
        
        return lst_result(basis, base, y_phys, F, A, B)
    end
    
    function dmat(N) :: chebyshev_basis
        num = N-1
        D0 = cos.((0:N-1)' .*(0:N-1) .* π./ num)
        
        # create higher derivative matrices
        D1 = hcat(zeros((N)), D0[:, 1], 4*D0[:, 2])
        D2 = hcat(zeros((N, 2)),        4*D1[:, 2])
        D3 = zeros((N, 3))
        D4 = zeros((N, 3))

        for j in range(3, num)
            D1 = hcat(D1, 2*j*D0[:, j]+j*D1[:, j-1]/(j-2))
            D2 = hcat(D2, 2*j*D1[:, j]+j*D2[:, j-1]/(j-2))
            D3 = hcat(D3, 2*j*D2[:, j]+j*D3[:, j-1]/(j-2))
            D4 = hcat(D4, 2*j*D3[:, j]+j*D4[:, j-1]/(j-2))
        end

        return chebyshev_basis(D0, D1, D2, D3, D4)
    end

    function bounded_base(y_phys,N,bf) :: baseflow
        if bf == 1 # Couette flow
    
            U   = y_phys
            Up  = ones((N, 1))
            Upp = zeros((N, 1))
            T   = y_phys
            Tp  = zeros((N, 1))
            
        elseif bf == 2  # Poiseuille flow
            
            U   = (1 .- y_phys.^2)
            Up  = -2 .* y_phys
            Upp = -2 .* ones((N, 1))
            T   = zeros((N, 1))
            Tp  = zeros((N, 1))
            
        elseif bf == 3 # quiescent flow
            
            U   = zeros((N, 1))
            Up  = zeros((N, 1))
            Upp = zeros((N, 1))
            T   = y_phys
            Tp  = ones((N, 1))
        else
            print("Need to select bf = 1 or 2")
        end
            
        return baseflow(U,Up,Upp,T,Tp)
    end

    function Operator(kx,kz,R,Pr,Ri, base::baseflow, basis::chebyshev_basis)
        N = size(base.U)[1]
        
        k2 = kx^2 + kz^2 
        M  =  ones((1, N))
        er = -200*1im

        A = -1im*kx .* (base.U*M) .* basis.D0 + (1/R/Pr) .*(basis.D2 .- k2.*basis.D0)
        B = basis.D0 .* 1


        A[1, :] .= 0
        A[N, :] .= 0

        A[1, :] .= er*basis.D0[1, :]
        A[N, :] .= er*basis.D0[N, :]
        
        return A, B
    end

end