using Quaternions 

type S3Grid
  depth::Int32
  vertices
  tetra
  tetra_levels
end

function InitS3Grid(s3::S3Grid)
  println("Initializing S3 Grid")
  # matrices are column major - so in order to be able to resize
  # efficiently we have one data-point per column
  s3.vertices = zeros(Float64, 4, 120)
  i = 1
  for a in [-1,1], b in [-1,1], c in [-1,1], d in [-1, 1]
    s3.vertices[:,i] = [a, b, c, d]
    i += 1
  end
  for j in 1:4, a in [-2, 2]
    s3.vertices[j,i] = a
    i+=1
  end
  # Golden Ratio
  phi = (1 + sqrt(5)) * 0.5
  # iterate over all *even* permutations
  # http://mathworld.wolfram.com/EvenPermutation.html
  for perm in [ (1,2,3,4), (1,3,4,2), (1,4,2,3), (2,1,4,3),
    (2,3,1,4), (2,4,3,1), (3,1,2,4), (3,2,4,1), (3,4,1,2),
    (4,1,3,2), (4,2,1,3), (4,3,2,1)]
    for a in [-1.,1.], b in [-1.,1.], c in [-1.,1.]
      s3.vertices[perm[1],i] = a*phi
      s3.vertices[perm[2],i] = b
      s3.vertices[perm[3],i] = c/phi
      s3.vertices[perm[4],i] = 0.
      i+=1
    end
  end
  s3.vertices *= 0.5
  # Tada: all of the vertices are unit length and hence \in S^3
  G = ones(Float64, 120, 120) * 99
  for i in 1:120, j in 1:120
    if j != i
      q1 = Quaternion(s3.vertices[1,i], s3.vertices[2,i],
        s3.vertices[3,i], s3.vertices[4,i])
      q2 = Quaternion(s3.vertices[1,j], s3.vertices[2,j],
        s3.vertices[3,j], s3.vertices[4,j])
      dq = normalize(q1 / q2)
      G[j,i] = angle(dq)
    end
  end
  minAngle = sort(G[:,1])[1]
  G[G .< minAngle - 1e-6] = -1
  G[G .> minAngle + 1e-6] = -1
  s3.tetra = zeros(Int32, 4, 600)
  i = 1
  for comb in combinations(1:120, 4) 
    if G[comb[1], comb[2]] > 0 && G[comb[1], comb[3]] > 0 && 
     G[comb[1], comb[4]] > 0 && G[comb[2], comb[3]] > 0 && 
     G[comb[2], comb[4]] > 0 && G[comb[3], comb[4]] > 0
      s3.tetra[:,i] = [comb[1], comb[2], comb[3], comb[4]]
      i += 1
    end
  end
  i == 601 && println("Wohoo found all the 600 Tetrahedra.")
  push!(s3.tetra_levels, 0)
  push!(s3.tetra_levels, 600)
end

function normalized(a)
  a /= sqrt(sum(a.^2))
end

function SubdivideOnce(s3::S3Grid)
  n_vertices = size(s3.vertices)[2]
  n_tetra = size(s3.tetra)[2]
  s3.vertices = reshape(resize!(s3.vertices[:], (n_vertices+n_tetra*6)*4), 
    4, n_vertices+n_tetra*6)
  s3.tetra = reshape(resize!(s3.tetra[:], (n_tetra*(8+1))*4),
    4, n_tetra*(8+1))
  push!(s3.tetra_levels, n_tetra * (8 + 1))
#  @sync @parallel for i in 1:n_tetra
  for i in 1:n_tetra
    i0 = s3.tetra[1,i]
    i1 = s3.tetra[2,i]
    i2 = s3.tetra[3,i]
    i3 = s3.tetra[4,i]
    i01 = n_vertices + (i-1)*6 + 1 
    i12 = n_vertices + (i-1)*6 + 2
    i20 = n_vertices + (i-1)*6 + 3
    i03 = n_vertices + (i-1)*6 + 4
    i13 = n_vertices + (i-1)*6 + 5
    i23 = n_vertices + (i-1)*6 + 6
    s3.vertices[:,i01] = normalized(s3.vertices[:,i0] + s3.vertices[:,i1])
    s3.vertices[:,i12] = normalized(s3.vertices[:,i1] + s3.vertices[:,i2]) 
    s3.vertices[:,i20] = normalized(s3.vertices[:,i2] + s3.vertices[:,i0]) 
    s3.vertices[:,i03] = normalized(s3.vertices[:,i0] + s3.vertices[:,i3]) 
    s3.vertices[:,i13] = normalized(s3.vertices[:,i1] + s3.vertices[:,i3]) 
    s3.vertices[:,i23] = normalized(s3.vertices[:,i2] + s3.vertices[:,i3]) 
    s3.tetra[:,n_tetra + (i-1)*8 + 1] = [i0, i01, i20, i03]
    s3.tetra[:,n_tetra + (i-1)*8 + 2] = [i01, i1, i12, i13]
    s3.tetra[:,n_tetra + (i-1)*8 + 3] = [i12, i2, i20, i23]
    s3.tetra[:,n_tetra + (i-1)*8 + 4] = [i03, i13, i23, i3]
    s3.tetra[:,n_tetra + (i-1)*8 + 5] = [i01, i13, i03, i20]
    s3.tetra[:,n_tetra + (i-1)*8 + 6] = [i01, i12, i13, i20]
    s3.tetra[:,n_tetra + (i-1)*8 + 7] = [i23, i20, i12, i13]
    s3.tetra[:,n_tetra + (i-1)*8 + 8] = [i23, i03, i20, i13]
  end
end

s3 = S3Grid(0, [], [], Int64[])
InitS3Grid(s3)
for lvl in 1:5
  println("Subdivide at level $lvl")
  SubdivideOnce(s3)
  println(s3.tetra_levels)
end
