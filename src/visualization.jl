function RecipesBase.plot(pl::AbstractPointSet; radius::Number = lattice(pl).radius, arrow=true)
    lpts, pts = points(lattice(pl)), points(pl)
    window = repeat(Float64[-radius radius], dim(pl))
    if dim(pl) == 2
        lattice_x = [p[1] for p in lpts]
        lattice_y = [p[2] for p in lpts]
        points_x = [p[1] for p in pts]
        points_y = [p[2] for p in pts]

        p = scatter(lattice_x, lattice_y, color=:red, markersize=3, label="Grid")

        # Add arrows from lattice to perturbed points
        if arrow
            for i in 1:length(pts)
                quiver!([lpts[i][1]], [lpts[i][2]],
                       quiver=([pts[i][1] - lpts[i][1]],
                              [pts[i][2] - lpts[i][2]]),
                       color=:gray, alpha=0.3)
            end
        end

        scatter!(points_x, points_y,
                color=:blue, markersize=3, label="Perturbed")
        xlims!(window[1,1], window[1,2])
        ylims!(window[2,1], window[2,2])
        xlabel!("x")
        ylabel!("y")
        plot!(aspect_ratio=:equal)

        return p
    elseif dim(pl) == 3
        lattice_x = [p[1] for p in lpts]
        lattice_y = [p[2] for p in lpts]
        lattice_z = [p[3] for p in lpts]
        points_x = [p[1] for p in pts]
        points_y = [p[2] for p in pts]
        points_z = [p[3] for p in pts]

        p = plot3d()

        # Connection lines
        for i in 1:length(pts)
            plot3d!([lpts[i][1], pts[i][1]],
                   [lpts[i][2], pts[i][2]],
                   [lpts[i][3], pts[i][3]],
                   color=:gray, alpha=0.3, label="", linewidth=0.5)
        end

        # Grid points
        scatter3d!(lattice_x, lattice_y, lattice_z,
                  color=:red, label="Grid", markersize=2)

        # Perturbed points
        scatter3d!(points_x, points_y, points_z,
                  color=:blue, label="Perturbed", markersize=2)
        xlims!(window[1,1], window[1,2])
        ylims!(window[2,1], window[2,2])
        zlims!(window[3,1], window[3,2])
        xlabel!("x")
        ylabel!("y")
        zlabel!("z")
        plot!(aspect_ratio=:equal)

        return p
    end
end