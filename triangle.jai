Triangle :: struct {
    p0, p1, p2: Vector3; // rename to a b c later probably.
    #overlay (p0) positions : [3] Vector3;
}

bounds :: (t : $T/interface Triangle) -> AABB {
    min := min(min(t.p0, t.p1), t.p2);
    max := max(max(t.p0, t.p1), t.p2);
    return .{min, max};
}

bounds :: (ts : [] $T/interface Triangle) -> AABB {
    assert(ts.count > 0);
    root := bounds(ts[0]);
    for ts merge(*root, bounds(it));
    return root;
}

Projection_Interval :: struct {
    min, max : float32;
}

project_triangle_onto_axis :: (a : Vector3, b : Vector3, c : Vector3, axis : Vector3) -> Projection_Interval {
    projection_a := dot(a, axis);
    projection_b := dot(b, axis);
    projection_c := dot(c, axis);

    projection_min := projection_a;
    projection_max := projection_a;
    if projection_b < projection_min projection_min = projection_b;
    if projection_b > projection_max projection_max = projection_b;
    if projection_c < projection_min projection_min = projection_c;
    if projection_c > projection_max projection_max = projection_c;

    return .{projection_min, projection_max};
}

project_aabb_onto_axis :: (box : AABB, axis : Vector3) -> Projection_Interval {
    center := (box.min + box.max) * 0.5;
    half_size := (box.max - box.min) * 0.5;
    projected_center := dot(center, axis);
    projected_radius :=
        half_size.x * abs(axis.x) +
        half_size.y * abs(axis.y) +
        half_size.z * abs(axis.z);

    return .{
        projected_center - projected_radius,
        projected_center + projected_radius,
    };
}

projection_intervals_overlap :: (a : Projection_Interval, b : Projection_Interval) -> bool {
    return a.min <= b.max && b.min <= a.max;
}

triangle_and_aabb_overlap_on_axis :: (
    a : Vector3,
    b : Vector3,
    c : Vector3,
    box : AABB,
    axis : Vector3,
) -> bool {
    triangle_interval := project_triangle_onto_axis(a, b, c, axis);
    box_interval := project_aabb_onto_axis(box, axis);
    return projection_intervals_overlap(triangle_interval, box_interval);
}

triangle_intersects_aabb :: (a : Vector3, b : Vector3, c : Vector3, box : AABB) -> bool  {
    box_axes : [3] Vector3 = .[
        .{1, 0, 0},
        .{0, 1, 0},
        .{0, 0, 1},
    ];

    triangle_edges : [3] Vector3 = .[
        b - a,
        c - b,
        a - c,
    ];

    // The box face normals also provide the inexpensive triangle-AABB bounds test.
    for axis : box_axes {
        if !triangle_and_aabb_overlap_on_axis(a, b, c, box, axis) return false;
    }

    triangle_normal := cross(triangle_edges[0], triangle_edges[1]);
    if !triangle_and_aabb_overlap_on_axis(a, b, c, box, triangle_normal) return false;

    // These nine axes cover possible edge-edge separation. Parallel edges produce
    // a zero axis, whose two projection intervals both collapse to zero and overlap.
    for triangle_edge : triangle_edges {
        for box_axis : box_axes {
            edge_cross_axis := cross(triangle_edge, box_axis);
            if !triangle_and_aabb_overlap_on_axis(a, b, c, box, edge_cross_axis) return false;
        }
    }

    return true;
}

triangle_intersects_aabb :: (triangle : $T/interface Triangle, box : AABB) -> bool {
    return triangle_intersects_aabb(triangle.p0, triangle.p1, triangle.p2, box);
}
