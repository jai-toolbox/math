#import "Math";

AABB :: Axis_Aligned_Bounding_Box;

Axis_Aligned_Bounding_Box :: struct {
    min, max : Vector3;
}


/*

in the context of doing minimums and maximums, then the identity element
is no longer zero like it is in usual addition.

for the min binary operation, then min(x, infinity) = x
and for the max we have max(x, -infinity) = x
because computers have a finite fixed greatest value.



*/

MERGE_IDENTITY :: AABB.{
    .{FLOAT32_MAX, FLOAT32_MAX, FLOAT32_MAX},
    .{-FLOAT32_MAX, -FLOAT32_MAX, -FLOAT32_MAX},
};

size :: (a : AABB) -> Vector3 {
    return a.max - a.min;
}

longest_axis :: (a : AABB) -> Axis {
    extent := size(a);
    if extent.x > extent.y {
        if extent.x > extent.z {
            return .X;
        } else {
            return .Z;
        }
    } else { // extent.x <= extent.y
        if extent.y > extent.z {
            return .Y;
        } else {
            return .Z;
        }
    }
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

disjoint :: (a: AABB, b: AABB) -> bool {
    return a.max.x < b.min.x || b.max.x < a.min.x
        || a.max.y < b.min.y || b.max.y < a.min.y
        || a.max.z < b.min.z || b.max.z < a.min.z;
}

intersects :: (a: AABB, b: AABB) -> bool {
    return !disjoint(a, b);
}

inside :: (a: AABB, p: Vector3) -> bool {
    return a.min.x <= p.x && p.x <= a.max.x
        && a.min.y <= p.y && p.y <= a.max.y
        && a.min.z <= p.z && p.z <= a.max.z;
}

/*
expand the given aabb to contain p

definitionally:

min(a.min, p) <= a.min, p
max(a.max, p) >= a.max, p

thus the following aabb contains the
old aabb and p.

*/

merge :: (a : *AABB, p : Vector3) {
    a.min = min(a.min, p);
    a.max = max(a.max, p);
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
