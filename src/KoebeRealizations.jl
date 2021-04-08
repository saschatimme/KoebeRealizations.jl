module KoebeRealizations

export EuclideanMedialGraph, layout, plot_circle_packing, koebe_realization

using LinearAlgebra
using ModelingToolkit: @variables, ~, NonlinearSystem, NonlinearProblem
import NonlinearSolve
import Polymake
import Plots

"""
    vertex_facet_to_incidence(J)

Given a vertex-facet incidence matrix `J` where `J[i,j] == true` if the `i`-th face
contains the `j`-th vertex, this function returns a matrix `M` with `M[i,j] == true` if
there is an endge between vertex `i` and `j`.
"""
function vertex_facet_to_incidence(J)
    nfacets, nvertices = size(J)
    incidence = falses(nvertices, nvertices)
    for i = 1:nvertices, j = i+1:nvertices
        found_first = false
        for cont = 1:nfacets
            (J[cont, i] == 1 && J[cont, j] == 1) || continue

            if found_first
                incidence[i, j] = incidence[j, i] = 1
                break
            else
                found_first = true
            end
        end
    end
    incidence
end

"Given a face, e.g., [1,3,6] construct the boundary edge cycle [(6, 1), (1,3), (3,6)]"
function cycle(face)
    E = NTuple{2,Int}[]
    f = face[end]
    for g in face
        push!(E, (f, g))
        f = g
    end
    E
end

"[(6, 1), (1,3), (3,6)] -> [(6,3), (3, 1), (1, 6)]"
rev(cycle) = map(reverse, reverse(cycle))

"""
    EuclideanMedialGraph(vertex_facet_incidence)

Given a vertex-facet incidence matrix it constructs the Euclidean medial graph where the
vertex corresponding to the first edge is removed.
"""
struct EuclideanMedialGraph
    vertex_facet::BitMatrix
    incidence::BitMatrix
    faces::Vector{Vector{Int}}
    vertex_face_indicator::Vector{Bool}
    face_cycles::Vector{Vector{NTuple{2,Int}}}
    polytope_edges::Vector{NTuple{2,Int}}
    oriented_edges::Dict{NTuple{2,Int},Int}
end

function EuclideanMedialGraph(vertex_facet)
    J = vertex_facet
    M = vertex_facet_to_incidence(J)

    # STEP 1: a) Collect all edges of the polytope -> will become vertices
    #         b) Assemble faces of the medial graph
    #
    # a) collect edges
    edges = NTuple{2,Int}[]
    for i = 1:size(M, 1), j = i+1:size(M, 1)
        M[i, j] && push!(edges, (i, j))
    end

    # b) Assemble faces
    # Construct faces for each vertex
    vertex_faces = map(1:size(J, 2)) do t
        edges_containing_vertex = tuple.(t, findall(M[:, t]))
        f = map(edges_containing_vertex) do e
            findfirst(f -> minmax(e...) == f, edges)
        end

        # We collected all vertices in the face but we need to store the vertices in an order
        # such that `cycle(f)` describes the boundary edges of `f`.
        boundary = Int[]
        push!(boundary, popfirst!(f))
        faces_containing_vertex = findall(J[:, t])
        while !isempty(f)
            i, j = edges[boundary[end]]
            ind = popfirst!(f)
            k, l = edges[ind]
            # check if there exists a facet containing that vertex which contains both the edges.
            if any(
                id -> J[id, i] && J[id, j] && J[id, k] && J[id, l],
                faces_containing_vertex,
            )
                push!(boundary, ind)
            else
                push!(f, ind)
            end
        end
        boundary
    end

    # Construct faces for each face
    face_faces = map(1:size(J, 1)) do t
        V = findall(J[t, :])
        f = Int[]
        for i = 1:length(V), j = i+1:length(V)
            k, l = V[i], V[j]
            if M[k, l]
                push!(f, findfirst(e -> e == minmax(k, l), edges))
            end
        end

        # We collected all vertices in the face but we need to store the vertices in an order
        # such that `cycle(f)` describes the boundary edges of `f`.
        boundary = Int[]
        push!(boundary, popfirst!(f))
        while !isempty(f)
            b = edges[boundary[end]]
            ind = popfirst!(f)
            k, l = edges[ind]
            # check if they share a vertex
            if k ∈ b || l ∈ b
                push!(boundary, ind)
            else
                push!(f, ind)
            end
        end

        boundary
    end


    faces = Vector{Int}[]
    vertex_face_indicator = Bool[]
    for f in vertex_faces
        if all(!isone, f)
            push!(faces, f)
            push!(vertex_face_indicator, true)
        end
    end
    for f in face_faces
        if all(!isone, f)
            push!(faces, f)
            push!(vertex_face_indicator, false)
        end
    end

    # STEP 2: The faces of our graph should have a consistent orientation.
    #         Otherwise the layout algorithm becomes a mess.
    #

    # To guarantee consistent orientation we keep a set dict S where the keys are all
    # oriented edges and the values are the corresponding face
    S = Dict{NTuple{2,Int},Int}()
    # The first face is by definition consistent
    conistent_faces = faces[1]
    for e in cycle(faces[1])
        S[e] = 1
    end

    # need to process all others. put them in a queue
    faces_to_process = collect(2:length(faces))
    # empty queue until we are done
    limit_counter = 1
    limit = length(faces_to_process)^3

    while !isempty(faces_to_process)
        k = popfirst!(faces_to_process)
        face = faces[k]

        # now we build the boundary cycle in both orientations
        c = cycle(face)
        rc = rev(c)
        # if the orientation is correct as it is, then an edge of the reverse cycle
        # needs to be in S
        if any(e -> haskey(S, e), rc)
            for e in c
                S[e] = k
            end
            # maybe we need to reverse the orientation?
        elseif any(e -> haskey(S, e), c)
            for e in rc
                S[e] = k
            end
            reverse!(face)
            # If we are here, then no neighboring faces are processed so far
            # -> put to the end of the queue again
        else
            push!(faces_to_process, k)
        end
        limit_counter += 1
        if limit_counter > limit
            @warn("Aborted orientation algorithm.")
            break
        end
    end

    return EuclideanMedialGraph(
        BitMatrix(J),
        BitMatrix(M),
        faces,
        vertex_face_indicator,
        cycle.(faces),
        edges,
        S,
    )
end

function is_boundary_edge(G::EuclideanMedialGraph, e::NTuple{2,Int})
    i, j = e
    !haskey(G.oriented_edges, e) || !haskey(G.oriented_edges, (j, i))
end

is_interior_edge(G, e) = !is_boundary_edge(G, e)

is_boundary_face(G::EuclideanMedialGraph, k) =
    any(e -> is_boundary_edge(G, e), G.face_cycles[k])

is_vertex_face(G, k) = G.vertex_face_indicator[k]

function adjacent_face_indices(G::EuclideanMedialGraph, k)
    faces = Int[]
    for (i, j) in G.face_cycles[k]
        if haskey(G.oriented_edges, (j, i))
            push!(faces, G.oriented_edges[(j, i)])
        end
    end
    faces
end

"Computes the radii of the circles corresponding to the faces of `G`"
function circle_radii(G::EuclideanMedialGraph)
    @variables ρ[1:length(G.faces)]

    eqs = map(enumerate(G.faces)) do (i, f)
        s = is_boundary_face(G, i) ? π : 2π
        for j in adjacent_face_indices(G, i)
            s -= 2 * atan(tanh((-ρ[i] + ρ[j]) / 2)) + π / 2
        end
        0 ~ s
    end
    sys = NonlinearSystem(eqs[2:end], ρ[2:end], ρ[1:1])
    # random guess, doesn't matter since convex
    guess = ρ[2:end] .=> rand.()
    # fix ρ_f[3] to 0
    prob = NonlinearProblem(
        sys,
        guess,
        [0];
        # analytical jacobian
        jac = true,
    )
    res = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(); atol = 1e-12)

    #no radius must be too big
    return 0.5 * exp.([0; res.u])
end


"Rotate `(x,y)` by `θ` counter-clockwise"
rotate((x, y), θ) = (x * cos(θ) - y * sin(θ), x * sin(θ) + y * cos(θ))

function layout(G::EuclideanMedialGraph)
    J = G.vertex_facet
    nedges = length(G.polytope_edges)
    faces = G.faces
    face_cycles = G.face_cycles
    cycle_face_map = G.oriented_edges

    face_coords = Union{NTuple{2,Float64},Nothing}[nothing for _ in faces]
    vertex_coords = Union{NTuple{2,Float64},Nothing}[nothing for _ = 1:nedges]
    face_radii = circle_radii(G)

    face_coords[1] = (0, 0)

    for (i, j) in face_cycles[1]
        if is_interior_edge(G, (i, j))
            vertex_coords[i] = (0, face_radii[1])
            break
        end
    end

    # We perform the layout algorithm in a depth-first way.
    # So we need to recurse
    iterative_layout =
        (face_id, level) -> begin

            for (i, j) in face_cycles[face_id]
                haskey(G.oriented_edges, (j, i)) || continue
                isnothing(vertex_coords[i]) && isnothing(vertex_coords[j]) && continue
                next_edge = (i, j)
                next_face_id = cycle_face_map[(j, i)]
                if !isnothing(vertex_coords[i]) &&
                   !isnothing(vertex_coords[j]) &&
                   !isnothing(face_coords[next_face_id])
                    continue
                end

                #              * p_j
                #
                #    * f                   * g
                #
                #              * p_i

                p_i, p_j = vertex_coords[i], vertex_coords[j]
                f = face_coords[face_id]
                g = face_coords[next_face_id]

                φ_e = let r_l = face_radii[face_id], r_r = face_radii[next_face_id]
                    imag(log((r_l - r_r * cis(pi * (-0.5))) / (r_l - r_r * cis(pi * 0.5))) / 2)
                end

                s = face_radii[next_face_id] / face_radii[face_id]
                if !isnothing(p_i)
                    g = rotate(f .- p_i, -π / 2) .* s .+ p_i
                    vertex_coords[j] = rotate(p_i .- f, 2φ_e) .+ f
                else
                    g = rotate(f .- p_j, π / 2) .* s .+ p_j
                    vertex_coords[i] = rotate(p_j .- f, -2φ_e) .+ f
                end
                face_coords[next_face_id] = g

                iterative_layout(next_face_id, level + 1)
            end
            nothing
        end
    # Start recursion
    iterative_layout(1, 0)


    return (
        vertex_coordinates = map(xy -> isnothing(xy) ? xy : [xy...], vertex_coords),
        face_coords = map(xy -> isnothing(xy) ? xy : [xy...], face_coords),
        face_radii = face_radii,
    )
end

# visualization

function circle_shape(x, y, r)
    θ = range(0, 2π, length = 500)
    x .+ r .* sin.(θ), y .+ r .* cos.(θ)
end

function draw_circle!(p, x, y, r; color = :dodgerblue)
    Plots.plot!(
        p,
        circle_shape(x, y, r),
        series_type = [:shape],
        fillalpha = 0,
        linecolor = color,
    )
end

"""
    plot_circle_packing(vertex_facet)
    plot_circle_packing(::EuclideanMedialGraph)

Visualize the (euclidean) circle packing.
"""
function plot_circle_packing(G::EuclideanMedialGraph; show_vertices = false)
    vertex_coords, face_coords, face_radii = layout(G)

    p = Plots.plot(; aspect_ratio = :equal, grid = false, legend = false)
    ind = 1:length(face_coords)
    for (pt, r, k) in zip(face_coords[ind], face_radii[ind], ind)
        isnothing(pt) && continue
        (x, y) = pt
        color = is_vertex_face(G, k) ? :indianred : :dodgerblue
        draw_circle!(p, x, y, r; color = color)
    end
    if show_vertices
        vc = filter(!isnothing, vertex_coords)
        Plots.scatter!(p, first.(vc), last.(vc); color = :black)
    end
    p
end
function plot_circle_packing(vertex_facet)
    plot_circle_packing(EuclideanMedialGraph(vertex_facet))
end

# 2D -> 3D

function stereo_proj((x, y))
    den = x^2 + y^2 + 1
    return (2 * x / den, 2 * y / den, (den - 2) / den)
end

function tangency_points_on_sphere(vertex_coordinates)
    [[[0, 0, 1]]; stereo_proj.(vertex_coordinates[2:end])]
end

function point_minimal_distance_sum(tangency_points, G::EuclideanMedialGraph)
    M = G.incidence
    nvertices = size(M, 1)
    nedges = length(tangency_points)

    @variables x[1:3]

    Spring_Functional = map(1:3) do k
        s = 0
        for (pt, (i, j)) in zip(tangency_points, G.polytope_edges)
            s += (-pt[k]) / (1 - sum(pt .* x)) + x[k] / (1 - sum(x .^ 2))
        end
        0 ~ s
    end

    sys = NonlinearSystem(Spring_Functional, x[1:3], [])

    while true
        guess = x[1:3] .=> rand.()

        #I want my guess to be inside the ball
        while sum(guess[1][2]^2 + guess[2][2]^2 + guess[3][2]^2) > 1
            guess = x[1:3] .=> rand.()
        end

        prob = NonlinearProblem(
            sys,
            guess,
            # analytical jacobian
            jac = true,
        )
        res = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(); tol = 1e-12)

        if sum(res.u .^ 2) < 1
            return res.u
        end
    end
end

function projective_transformation(tangency_points, G::EuclideanMedialGraph)
    pt_min = [1; point_minimal_distance_sum(tangency_points, G)]
    e₁ = [1, 0, 0, 0]
    e₂ = [0, 1, 0, 0]
    e₃ = [0, 0, 1, 0]
    e₄ = [0, 0, 0, 1]
    L = [-e₁ e₂ e₃ e₄]

    # We perform Gram-Schmidt with respect to the L scalar product
    norm_pt_min = sqrt(-pt_min' * (L * pt_min))
    pt_min ./= norm_pt_min

    b₂ = e₂ .+ (e₂' * (L * pt_min)) .* pt_min
    sign_norm_2 = sign(b₂' * (L * b₂))
    norm_2 = sqrt(abs(b₂' * (L * b₂)))
    b₂ = sign_norm_2 * b₂ / norm_2

    b₃ = e₃ .+ (e₃' * (L * pt_min)) .* pt_min .- (e₃' * (L * b₂)) .* b₂
    norm_3 = sqrt(abs(b₃' * (L * b₃)))
    sign_norm_3 = sign(b₃' * (L * b₃))
    b₃ = sign_norm_3 * b₃ / norm_3

    b₄ =
        e₄ .+ (e₄' * (L * pt_min)) .* pt_min .- (e₄' * (L * b₂)) .* b₂ .-
        (e₄' * (L * b₃)) .* b₃
    norm_4 = sqrt(abs(b₄' * (L * b₄)))
    sign_norm_4 = sign(b₄' * (L * b₄))
    b₄ = sign_norm_4 * b₄ / norm_4

    # We return the inverse of the matrix whose columns
    # are given by the computed orthonormal basis
    return L * ([pt_min b₂ b₃ b₄]' * L)
end

"""

"""
function koebe_realization(G::EuclideanMedialGraph)
    vertex_coordinates, face_coordinates, face_radii = layout(G)
    tangency_points = tangency_points_on_sphere(vertex_coordinates)
    nvertices = size(G.incidence, 1)

    # Let us find the vertices using the tangency points
    # solving linear systems
    E = zeros(3, 3)
    vertices = map(1:nvertices) do i
        indices = Int[]
        for j = 1:nvertices
            G.incidence[i, j] && push!(indices, j)
            length(indices) == 3 && break
        end

        for a = 1:3, b = 1:3
            e = minmax(i, indices[a])
            k = findfirst(x -> x == e, G.polytope_edges)
            E[a, b] = tangency_points[k][b]
        end
        E \ [1, 1, 1]
    end
    A = projective_transformation(tangency_points, G)
    V = (A * [ones(1, nvertices); reduce(hcat, vertices)])'
    V[:,2:4] ./ V[:,1]
end

function koebe_realization(vertex_facet::AbstractMatrix)
    koebe_realization(EuclideanMedialGraph(vertex_facet))
end
function koebe_realization(poly::Polymake.BigObject)
    pts = koebe_realization(poly.VERTICES_IN_FACETS)
    Polymake.polytope.Polytope(; POINTS = [ones(size(pts, 1), 1) pts])
end


end # module
