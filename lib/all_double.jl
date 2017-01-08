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

    ax::Float64
    ay::Float64
    az::Float64

    # ax_dx::Float64
    # ax_dy::Float64
    # ax_dz::Float64

    # ay_dx::Float64
    # ay_dy::Float64
    # ay_dz::Float64

    # az_dx::Float64
    # az_dy::Float64
    # az_dz::Float64

    function OctTreeNode(r::Number, midx::Number, midy::Number, midz::Number)
        n = new(r, midx, midy, midz, true, false, Particle(0.0, 0.0, 0.0, 0.0))
        n.lxlylz = n
        n.lxhylz = n
        n.hxlylz = n
        n.hxhylz = n
        n.lxlyhz = n
        n.lxhyhz = n
        n.hxlyhz = n
        n.hxhyhz = n
        n.ax=0.0
        n.ay=0.0
        n.az=0.0
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
    faststack1::Vector{OctTreeNode}
    faststack2::Vector{OctTreeNode}
    function OctTree(r::Number,
        midx::Number, midy::Number, midz::Number, n::Int64=100000)

        n = round(Int64,n*4.5)
        nodes = OctTreeNode[OctTreeNode(r, midx, midy, midz) for i in 1:n]
        new(nodes[1], 1, nodes, [OctTreeNode(0.0, 0.0, 0.0, 0.0) for i in 1:10000],[OctTreeNode(0.0, 0.0, 0.0, 0.0) for i in 1:10000])
    end
end

function clear!(h::OctTree)
    h.head.is_empty = true
    h.head.is_divided = false
    h.number_of_nodes_used = 1
    nothing
end

@inline function modify!(q::OctTreeNode, p::Particle)
    const total_mass = q.particle.m + p.m
    const newx = (q.particle.x*q.particle.m + p.x*p.m)/total_mass
    const newy = (q.particle.y*q.particle.m + p.y*p.m)/total_mass
    const newz = (q.particle.z*q.particle.m + p.z*p.m)/total_mass
    q.particle = Particle(newx, newy, newz, total_mass)
    nothing
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

@inline function initnode!(q::OctTreeNode, r::Number,
    midx::Number, midy::Number, midz::Number)

    q.r = r
    q.midx = midx
    q.midy = midy
    q.midz = midz
    q.is_empty = true
    q.is_divided = false
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
    initnode!(q.lxhylz, r2, q.midx-r2, q.midy+r2, q.midz-r2)
    initnode!(q.lxlylz, r2, q.midx-r2, q.midy-r2, q.midz-r2)
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

function insert!(t::OctTree, particles::Vector{Particle})
    for p in particles
        insert!(t, p)
    end
    nothing
end


immutable stop_cond_data
    should_stop::Bool
    ax1::Float64
    ay1::Float64
    az1::Float64
    ax2::Float64
    ay2::Float64
    az2::Float64
end

function prepare_accel!(t::OctTree, alpha2::Float64, eps2::Float64)
    curr_stack_ix = 1
    t.faststack1[1] = t.head
    t.faststack2[1] = t.head
    @inbounds while curr_stack_ix > 0
        q1 = t.faststack1[curr_stack_ix]
        q2 = t.faststack2[curr_stack_ix]
        curr_stack_ix -= 1
        if q1==q2
            if !q1.is_divided
                continue
            end
            # expand self interactions
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxhyhz
            t.faststack2[curr_stack_ix] = q2.hxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxhylz
            t.faststack2[curr_stack_ix] = q2.hxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxlyhz
            t.faststack2[curr_stack_ix] = q2.hxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxlylz
            t.faststack2[curr_stack_ix] = q2.hxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhyhz
            t.faststack2[curr_stack_ix] = q2.hxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhylz
            t.faststack2[curr_stack_ix] = q2.hxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlyhz
            t.faststack2[curr_stack_ix] = q2.hxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlylz
            t.faststack2[curr_stack_ix] = q2.hxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxhylz
            t.faststack2[curr_stack_ix] = q2.hxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxlyhz
            t.faststack2[curr_stack_ix] = q2.hxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxlylz
            t.faststack2[curr_stack_ix] = q2.hxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhyhz
            t.faststack2[curr_stack_ix] = q2.hxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhylz
            t.faststack2[curr_stack_ix] = q2.hxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlyhz
            t.faststack2[curr_stack_ix] = q2.hxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlylz
            t.faststack2[curr_stack_ix] = q2.hxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxlyhz
            t.faststack2[curr_stack_ix] = q2.hxlyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxlylz
            t.faststack2[curr_stack_ix] = q2.hxlyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhyhz
            t.faststack2[curr_stack_ix] = q2.hxlyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhylz
            t.faststack2[curr_stack_ix] = q2.hxlyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlyhz
            t.faststack2[curr_stack_ix] = q2.hxlyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlylz
            t.faststack2[curr_stack_ix] = q2.hxlyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.hxlylz
            t.faststack2[curr_stack_ix] = q2.hxlylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhyhz
            t.faststack2[curr_stack_ix] = q2.hxlylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhylz
            t.faststack2[curr_stack_ix] = q2.hxlylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlyhz
            t.faststack2[curr_stack_ix] = q2.hxlylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlylz
            t.faststack2[curr_stack_ix] = q2.hxlylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhyhz
            t.faststack2[curr_stack_ix] = q2.lxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhylz
            t.faststack2[curr_stack_ix] = q2.lxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlyhz
            t.faststack2[curr_stack_ix] = q2.lxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlylz
            t.faststack2[curr_stack_ix] = q2.lxhyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxhylz
            t.faststack2[curr_stack_ix] = q2.lxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlyhz
            t.faststack2[curr_stack_ix] = q2.lxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlylz
            t.faststack2[curr_stack_ix] = q2.lxhylz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlyhz
            t.faststack2[curr_stack_ix] = q2.lxlyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlylz
            t.faststack2[curr_stack_ix] = q2.lxlyhz
            curr_stack_ix += 1
            t.faststack1[curr_stack_ix] = q1.lxlylz
            t.faststack2[curr_stack_ix] = q2.lxlylz
            continue
        end
        s = stop_cond(q1, q2, alpha2, eps2)
        if s.should_stop
            q1.ax += s.ax1
            q1.ay += s.ay1
            q1.az += s.az1
            q2.ax += s.ax2
            q2.ay += s.ay2
            q2.az += s.az2
            continue
        end
        if (q1.r > q2.r && q1.is_divided) || !q2.is_divided
            q1,q2 = q2,q1
        end

        # now q1 is the smaller node
        # we split the larger node and push it into the stack
        curr_stack_ix += 1
        t.faststack1[curr_stack_ix] = q1
        t.faststack2[curr_stack_ix] = q2.hxhyhz
        curr_stack_ix += 1
        t.faststack1[curr_stack_ix] = q1
        t.faststack2[curr_stack_ix] = q2.hxhylz
        curr_stack_ix += 1
        t.faststack1[curr_stack_ix] = q1
        t.faststack2[curr_stack_ix] = q2.hxlyhz
        curr_stack_ix += 1
        t.faststack1[curr_stack_ix] = q1
        t.faststack2[curr_stack_ix] = q2.hxlylz
        curr_stack_ix += 1
        t.faststack1[curr_stack_ix] = q1
        t.faststack2[curr_stack_ix] = q2.lxhyhz
        curr_stack_ix += 1
        t.faststack1[curr_stack_ix] = q1
        t.faststack2[curr_stack_ix] = q2.lxhylz
        curr_stack_ix += 1
        t.faststack1[curr_stack_ix] = q1
        t.faststack2[curr_stack_ix] = q2.lxlyhz
        curr_stack_ix += 1
        t.faststack1[curr_stack_ix] = q1
        t.faststack2[curr_stack_ix] = q2.lxlylz
    end
    nothing
end

@inline function stop_cond(q1::OctTreeNode, q2::OctTreeNode, alpha2::Float64, eps2::Float64)
    if isemptyleaf(q1) || isemptyleaf(q2)
        return stop_cond_data(true, 0.0,0.0,0.0, 0.0,0.0,0.0)
    end

    const dx = q2.particle.x - q1.particle.x
    const dy = q2.particle.y - q1.particle.y
    const dz = q2.particle.z - q1.particle.z

    const dx2 = dx*dx
    const dy2 = dy*dy
    const dz2 = dz*dz
    const dr2 = dx2 + dy2 + dz2

    const smthdr2 = eps2+dr2
    const smthdr = sqrt(smthdr2)

    if !q1.is_empty && !q2.is_empty
        const smthdr32 = smthdr2*smthdr
        const denom1 = smthdr32/q2.particle.m
        const denom2 = -smthdr32/q1.particle.m
        return stop_cond_data(true, dx/denom1, dy/denom1, dz/denom1, dx/denom2, dy/denom2, dz/denom2)
    end

    const l = q1.r+q2.r
    if l/smthdr > alpha2
        if q1.is_divided || q2.is_divided
            return stop_cond_data(false, 0.0,0.0,0.0, 0.0,0.0,0.0)
        end
        return stop_cond_data(true, 0.0,0.0,0.0, 0.0,0.0,0.0)
    end

    const smthdr32 = smthdr2*smthdr
    const denom1 = smthdr32/q2.particle.m
    const denom2 = -smthdr32/q1.particle.m
    return stop_cond_data(true, dx/denom1, dy/denom1, dz/denom1, dx/denom2, dy/denom2, dz/denom2)
end;

function get_accel(t::OctTree, particle::Particle)
    q = t.head
    ax=q.ax
    ay=q.ay
    az=q.az
    while q.is_divided
        q = getsubnode(q, particle)
        ax += q.ax
        ay += q.ay
        az += q.az
    end        
    ax,ay,az
end







