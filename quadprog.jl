#
# Julia implementation of primal/dual QP algorithm
#
# Tested on Julia 0.3 (0.3.11)
#


function quadprog{T}(Q::Matrix{T},c::Vector{T},A::Matrix{T},b::Vector{T};tol=1E-8)
	n = length(c)
	m = length(b)

	x = zeros(T,n)
	s = ones(T,m)
	z = ones(T,m)

	dx = Array(T,n)
	ds = Array(T,m)
	dz = Array(T,m)

	rd = c - A'*z
	rp = s - A*x + b
	rs = (s.*z)

	mu = n/m

	iter = 0
	for iter = 1:250
		Gbar = Q + A'*diagm(z./s)*A
		Gdecomp = cholfact(Gbar)

		for i=1:2
			# Compute step
			tmp = (rs - z.*rp)./s
			dx = -(Gdecomp\(rd + A'*tmp))
			ds = A*dx - rp
			dz = -(rs - z.*ds)./s

			# Compute alpha
			az = -z./dz
			az = minimum(az[az .> 0.0])
			sz = -s./ds
			sz = minimum(sz[sz .> 0.0])
			alpha = min(sz,az)
			alpha = min(alpha,1.0)

			if (i==2)
				break
			end
			mu_aff  = dot(s+alpha*ds,z+alpha*dz)/m
			sigma   = (mu_aff/mu)^3

			rs += ds.*dz - sigma*mu
		end

		# Step
		alpha *= 0.95
		x += alpha*dx
		s += alpha*ds
		z += alpha*dz

		# Update residuals
		rd = Q*x + c - A'*z
		rp = s - A*x + b
		rs = s.*z
		mu = dot(s,z)/m

		# Convergence
		if ( (mu < tol) &&
			 (norm(rd) < tol) &&
			 (norm(rs) < tol) )
			break
		end
	end
	return x
end