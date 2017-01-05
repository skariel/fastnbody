using GeometricalPredicates
using Base.Threads

import GeometricalPredicates:getx, gety, getz

immutable Particle <: AbstractPoint3D
    x::Float64
    y::Float64
    z::Float64
    m::Float64
end
getx(p::Particle) = p.x
gety(p::Particle) = p.y
getz(p::Particle) = p.z

type OctTreeNode
    id::Int64
    r::Float64
    midx::Float64
    midy::Float64
    midz::Float64
    is_empty::Bool
    is_divided::Bool
    particle::Particle
    lxlylz::OctTreeNode
    lxhylz::OctTreeNode
    hxlylz::OctTreeNode
    hxhylz::OctTreeNode
    lxlyhz::OctTreeNode
    lxhyhz::OctTreeNode
    hxlyhz::OctTreeNode
    hxhyhz::OctTreeNode
    function OctTreeNode(r::Number, midx::Number, midy::Number, midz::Number)
        n = new(0, r, midx, midy, midz, true, false, Particle(0.0, 0.0, 0.0, 0.0))
        n.lxlylz = n
        n.lxhylz = n
        n.hxlylz = n
        n.hxhylz = n
        n.lxlyhz = n
        n.lxhyhz = n
        n.hxlyhz = n
        n.hxhyhz = n
        n
    end
end

@inline isleaf(q::OctTreeNode) = !q.is_divided
@inline isemptyleaf(q::OctTreeNode) = !q.is_divided && q.is_empty
@inline isfullleaf(q::OctTreeNode) = !q.is_empty

type OctTree
    head::OctTreeNode
    number_of_nodes_used::Int64
    nodes::Vector{OctTreeNode}
    faststack::Vector{OctTreeNode}
    @inbounds function OctTree(r::Number,
        midx::Number, midy::Number, midz::Number, n::Int64=100000)

        n = round(Int64,n*4.5)
        nodes = OctTreeNode[OctTreeNode(r, midx, midy, midz) for i in 1:n]
        new(nodes[1], 1, nodes, [OctTreeNode(0.0, 0.0, 0.0, 0.0) for i in 1:10000])
    end
end

function clear!(h::OctTree)
    h.head.id = 0
    h.head.is_empty = true
    h.head.is_divided = false
    h.number_of_nodes_used = 1
    nothing
end

function insert!(h::OctTree, particle::Particle)
    q = h.head
    while q.is_divided
        modify!(q, particle)
        q = getsubnode(q, particle)
    end
    while !q.is_empty
        const friend = q.particle
        divide!(h, q)
        modify!(q, friend)
        modify!(q, particle)
        q = getsubnode(q, particle)
    end
    q.particle = particle
    q.is_empty = false
    q
end

@inline function initnode!(q::OctTreeNode, r::Number,
    midx::Number, midy::Number, midz::Number)

    q.r = r
    q.midx = midx
    q.midy = midy
    q.midz = midz
    q.is_empty = true
    q.is_divided = false
    q.id = 0
end

@inline function divide!(h::OctTree, q::OctTreeNode)
    # this line doesnt work well when multithreaded...
    @assert length(h.nodes) - h.number_of_nodes_used >= 8

    # populate new nodes
    q.lxlylz = h.nodes[h.number_of_nodes_used+1]
    q.lxhylz = h.nodes[h.number_of_nodes_used+2]
    q.hxlylz = h.nodes[h.number_of_nodes_used+3]
    q.hxhylz = h.nodes[h.number_of_nodes_used+4]
    q.lxlyhz = h.nodes[h.number_of_nodes_used+5]
    q.lxhyhz = h.nodes[h.number_of_nodes_used+6]
    q.hxlyhz = h.nodes[h.number_of_nodes_used+7]
    q.hxhyhz = h.nodes[h.number_of_nodes_used+8]

    # set new nodes properties (dimensions etc.)
    const r2 = q.r/2.0
    initnode!(q.lxlylz, r2, q.midx-r2, q.midy-r2, q.midz-r2)
    initnode!(q.lxhylz, r2, q.midx-r2, q.midy+r2, q.midz-r2)
    initnode!(q.hxlylz, r2, q.midx+r2, q.midy-r2, q.midz-r2)
    initnode!(q.hxhylz, r2, q.midx+r2, q.midy+r2, q.midz-r2)
    initnode!(q.lxlyhz, r2, q.midx-r2, q.midy-r2, q.midz+r2)
    initnode!(q.lxhyhz, r2, q.midx-r2, q.midy+r2, q.midz+r2)
    initnode!(q.hxlyhz, r2, q.midx+r2, q.midy-r2, q.midz+r2)
    initnode!(q.hxhyhz, r2, q.midx+r2, q.midy+r2, q.midz+r2)

    # update tree and parent node
    h.number_of_nodes_used += 8

    q.is_divided = true
    if !q.is_empty
        # move particle in parent node to child node
        const sq = getsubnode(q, q.particle)
        sq.is_empty = false
        q.is_empty = true
        sq.particle = q.particle
        q.particle = Particle(0.0,0.0,0.0,0.0)
    end    
    q
end

@inline function getsubnode(q::OctTreeNode, particle::Particle)
    if particle.x<q.midx
        # lx
        if particle.y<q.midy
            # ly
            particle.z<q.midz && return q.lxlylz
            return q.lxlyhz
        else
            # hy
            particle.z<q.midz && return q.lxhylz
            return q.lxhyhz
        end
    else
        # hx
        if particle.y<q.midy
            # ly
            particle.z<q.midz && return q.hxlylz
            return q.hxlyhz
        else
            # hy
            particle.z<q.midz && return q.hxhylz
            return q.hxhyhz
        end
    end
end

function get_accel(t::OctTree, x::Float64, y::Float64, z::Float64,
                    alpha2::Float64, eps2::Float64)

    curr_stack_ix = 1
    t.faststack[1] = t.head
    ax=0.0
    ay=0.0
    az=0.0
    @inbounds while curr_stack_ix > 0
        q = t.faststack[curr_stack_ix]
        curr_stack_ix -= 1
        should_stop, tax, tay, taz = stop_cond(q, x,y,z, alpha2, eps2)
        if !should_stop
            # of course q should be divided by now
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.lxlylz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.lxlyhz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.lxhylz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.lxhyhz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.hxlylz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.hxlyhz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.hxhylz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.hxhyhz
        end
        ax += tax
        ay += tay
        az += taz
    end
    return ax,ay,az
end

function get_accel_rel(t::OctTree, x::Float64, y::Float64, z::Float64,
                    alpha2::Float64, eps2::Float64, old_acc2)

    curr_stack_ix = 1
    t.faststack[1] = t.head
    ax=0.0
    ay=0.0
    az=0.0
    @inbounds while curr_stack_ix > 0
        q = t.faststack[curr_stack_ix]
        curr_stack_ix -= 1
        should_stop, tax, tay, taz = stop_cond_rel(q, x,y,z, alpha2, eps2, old_acc2)
        if !should_stop
            # of course q should be divided by now
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.lxlylz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.lxlyhz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.lxhylz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.lxhyhz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.hxlylz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.hxlyhz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.hxhylz
            curr_stack_ix += 1
            t.faststack[curr_stack_ix] = q.hxhyhz
        end
        ax += tax
        ay += tay
        az += taz
    end
    return ax,ay,az
end

@inline function modify!(q::OctTreeNode, p::Particle)
    const total_mass = q.particle.m + p.m
    const newx = (q.particle.x*q.particle.m + p.x*p.m)/total_mass
    const newy = (q.particle.y*q.particle.m + p.y*p.m)/total_mass
    const newz = (q.particle.z*q.particle.m + p.z*p.m)/total_mass
    q.particle = Particle(newx, newy, newz, total_mass)
    nothing
end

@inline function stop_cond(q::OctTreeNode,
    x::Float64, y::Float64, z::Float64, alpha2::Float64, eps2::Float64)

    isemptyleaf(q) && return true, 0.0, 0.0, 0.0

    const dx = q.particle.x - x
    const dy = q.particle.y - y
    const dz = q.particle.z - z

    const dx2 = dx*dx
    const dy2 = dy*dy
    const dz2 = dz*dz
    const dr2 = dx2 + dy2 + dz2

    const l = 2.0*q.r
    const l2 = l*l

    if l2/dr2 > alpha2 && q.is_divided
        return false, 0.0, 0.0, 0.0 # we need to further open the node
    end

    const smthdr2 = eps2+dr2
    const smthdr = sqrt(smthdr2)
    const denom = smthdr2*smthdr/q.particle.m

    return true, dx/denom, dy/denom, dz/denom
end;

@inline function stop_cond_rel(q::OctTreeNode,
    x::Float64, y::Float64, z::Float64, alpha2::Float64, eps2::Float64,
    old_acc::Float64)

    isemptyleaf(q) && return true, 0.0, 0.0, 0.0

    const dx = q.particle.x - x
    const dy = q.particle.y - y
    const dz = q.particle.z - z

    const dx2 = dx*dx
    const dy2 = dy*dy
    const dz2 = dz*dz
    const dr2 = dx2 + dy2 + dz2

    const l = 2.0*q.r
    const l2 = l*l

    const smthdr2 = eps2+dr2
    const new_acc = q.particle.m/smthdr2


    if new_acc*l2/dr2 > alpha2*old_acc && q.is_divided
        return false, 0.0, 0.0, 0.0 # we need to further open the node
    end

    const smthdr = sqrt(smthdr2)
    const denom = smthdr2*smthdr/q.particle.m

    return true, dx/denom, dy/denom, dz/denom
end;

type ParallelOctTree
    trees::Vector{OctTree}
    tax::Array{Float64, 2}
    tay::Array{Float64, 2}
    taz::Array{Float64, 2}
    function ParallelOctTree(r::Number,
        midx::Number, midy::Number, midz::Number, n::Int64=100000)
        n4 = div(n,4)
        new([OctTree(r,midx,midy,midz,n4) for i in 1:nthreads()],
            zeros(nthreads(),n), zeros(nthreads(),n), zeros(nthreads(),n))
    end
end

function clear!(pt::ParallelOctTree)
    @threads for t in pt.trees
        clear!(t)
    end
    nothing
end

function insert!(pt::ParallelOctTree, particles::Vector{Particle})
    chunks = Array(UnitRange{Int64}, nthreads())
    li = round(Int64, linspace(1,length(particles),1+nthreads()))
    @inbounds for i in 1:(nthreads()-1)
        chunks[i] = li[i]:(li[i+1]-1)
    end
    @inbounds chunks[end] = li[end-1]:li[end]
    @threads for tid in 1:nthreads()
        const t = pt.trees[tid]
        @inbounds for i in chunks[tid]
            insert!(t, particles[i])
        end
    end
    nothing
end

function get_accel_all_particles!(pt::ParallelOctTree,
    particles::Vector{Particle}, alpha2::Float64, eps2::Float64,
        iax::Vector{Float64}, iay::Vector{Float64}, iaz::Vector{Float64})

    @threads for tid in 1:nthreads()
        const t = pt.trees[tid]
        @inbounds for i in 1:length(particles)
            const p = particles[i]
            ax,ay,az = get_accel(t, p.x, p.y, p.z, alpha2, eps2)
            pt.tax[tid,i] = ax
            pt.tay[tid,i] = ay
            pt.taz[tid,i] = az
        end
    end

    @inbounds for i in 1:length(particles)
        iax[i] = pt.tax[1,i]
        iay[i] = pt.tay[1,i]
        iaz[i] = pt.taz[1,i]
        for tid in 2:nthreads()
            iax[i] += pt.tax[tid,i]
            iay[i] += pt.tay[tid,i]
            iaz[i] += pt.taz[tid,i]
        end
    end
    nothing
end

function get_accel_all_particles_rel!(pt::ParallelOctTree,
    particles::Vector{Particle}, alpha2::Float64, eps2::Float64,
        iax::Vector{Float64}, iay::Vector{Float64}, iaz::Vector{Float64})

    @threads for tid in 1:nthreads()
        const t = pt.trees[tid]
        @inbounds for i in 1:length(particles)
            const p = particles[i]
            const old_acc = sqrt(iax[i]*iax[i]+iay[i]*iay[i]+iaz[i]*iaz[i])
            ax,ay,az = get_accel_rel(t, p.x, p.y, p.z, alpha2, eps2, old_acc)
            pt.tax[tid,i] = ax
            pt.tay[tid,i] = ay
            pt.taz[tid,i] = az
        end
    end

    @inbounds for i in 1:length(particles)
        iax[i] = pt.tax[1,i]
        iay[i] = pt.tay[1,i]
        iaz[i] = pt.taz[1,i]
        for tid in 2:nthreads()
            iax[i] += pt.tax[tid,i]
            iay[i] += pt.tay[tid,i]
            iaz[i] += pt.taz[tid,i]
        end
    end
    nothing
end
