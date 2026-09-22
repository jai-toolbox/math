#import "Math";

AABB :: Axis_Aligned_Bounding_Box;

Axis_Aligned_Bounding_Box :: struct {
    min, max : Vector3; 
}

size :: (a : AABB) -> Vector3 {
    return a.max - a.min;
}

merge :: (a : AABB, b : AABB) -> AABB {
    min := min(a.min, b.min);
    max := max(a.max, b.max);
    return .{min, max};
}

merge :: (a : *AABB, b : AABB) {
    a.min = min(a.min, b.min);
    a.max = max(a.max, b.max);
}

create_aabb_from_positions :: (positions: [..] Vector3) -> AABB {
    min_so_far := ONE3 * FLOAT32_MAX;
    max_so_far := ONE3 * -FLOAT32_MAX;

    for positions {
        min_so_far = min(min_so_far, it);
        max_so_far = max(max_so_far, it);
    }

    return .{min_so_far, max_so_far};
}

// bad function name, shouldn't be called change should have copy or something in it.
change_aabb_origin :: (aabb : AABB, position: Vector3) -> AABB {
    out : AABB;

    half_extents := (aabb.max - aabb.min) * 0.5;

    out.min = position - half_extents;
    out.max = position + half_extents;

    return out;
}

Ray :: struct {
    // note we didn't call it direction because that usually implies that length doesn't matter.
    origin, reach : Vector3;
}
