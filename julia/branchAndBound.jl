

# http://julia.readthedocs.org/en/latest/manual/control-flow/
# Vector allows pop! and push! i.e. LIFO
# PriorityQueue allows sorted list according to the key a value is
# stored under.

function ToDeg(angle)
  return angle * 180./pi
end

function BranchAndBound(nodes, objective, lowerBound, upperBound, branch)
  if size(nodes) == 0
    return
  end
  UpperBound = upperBound(nodes[1])
  LowerBound = lowerBound(nodes[1])
  BestNode = nodes[1]
  while abs(UpperBound - LowerBound) > 1e-3
    node = pop!(nodes)
    u = upperBound(node) 
    l = lowerBound(node) 
    if u > UpperBound
      UpperBound = u
      LowerBound = l
      BestNode = node
      println("New best node with $UpperBound")
      println(node)
      println(size(nodes))
      if size(nodes)[1] == 0
        push!(nodes, node)
      end
    else
      for branchedNode in branch(node)
        println((UpperBound, upperBound(branchedNode), size(nodes), branchedNode))
        if upperBound(branchedNode) >= UpperBound
          push!(nodes, branchedNode)
        end
      end
    end
    println((ToDeg(LowerBound), ToDeg(UpperBound), ToDeg(0.5*(UpperBound+LowerBound))))
  end
end

type AngularArea
  alphaStart::Float64
  alphaEnd::Float64
end

function branch(node::AngularArea)
  branchedNodes = Any[]
  push!(branchedNodes, AngularArea(node.alphaStart, (node.alphaEnd+node.alphaStart)*0.5))
  push!(branchedNodes, AngularArea((node.alphaEnd+node.alphaStart)*0.5, node.alphaEnd))
  return branchedNodes
end

function vMFOverlap(node::AngularArea)
  alpha = (node.alphaEnd+node.alphaStart)*0.5
  return vMFOverlap(alpha)
end

function vMFOverlap(alpha::Float64)
  R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)] 
  f = mu1 + R*mu2
  return sqrt(sum(f.^2))
end

#function vMFOverlapDerivative(alpha::Float64)
#  R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)] 
#  dR = [-sin(alpha) -cos(alpha); cos(alpha) -sin(alpha)] 
#  f = mu1 + R*mu2
#  sqrt(sum(f.^2))
#end

function vMFOverlapUpperBound(node::AngularArea)
  alphas = linspace(node.alphaStart, node.alphaEnd, 100)
  fs = map(vMFOverlap, alphas)
  return maximum(fs) - 1.e-4*(node.alphaEnd-node.alphaStart)
end

function vMFOverlapLowerBound(node::AngularArea)
  alphas = linspace(node.alphaStart, node.alphaEnd, 100)
  fs = map(vMFOverlap, alphas)
  return minimum(fs) + 1.e-4/(node.alphaEnd-node.alphaStart)
end

nodes = Any[]
push!(nodes, AngularArea(0., 2.*pi))
mu1 = [1.,0.]
mu2 = [0.,1.]

BranchAndBound(nodes, vMFOverlap, vMFOverlapLowerBound, vMFOverlapUpperBound, branch)
