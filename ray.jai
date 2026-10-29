// Slab test for a finite ray segment against an axis-aligned bounding box.
// The tested segment is origin + direction * t, where t is in
// [0, end_fraction].
ray_segment_intersects_aabb :: (
    origin: Vector3,
    direction: Vector3,
    bounds_min: Vector3,
    bounds_max: Vector3,
    end_fraction: float32,
) -> bool {
    near_fraction: float32 = 0;
    far_fraction := end_fraction;
    for axis: 0..2 {
        axis_origin := origin.component[axis];
        axis_direction := direction.component[axis];
        axis_minimum := bounds_min.component[axis];
        axis_maximum := bounds_max.component[axis];

        if abs(axis_direction) <= 0.00000001 {
            if axis_origin < axis_minimum || axis_origin > axis_maximum return false;
            continue;
        }

        inverse_direction := 1.0 / axis_direction;
        axis_near := (axis_minimum - axis_origin) * inverse_direction;
        axis_far := (axis_maximum - axis_origin) * inverse_direction;
        if axis_near > axis_far swap(*axis_near, *axis_far);
        near_fraction = max(near_fraction, axis_near);
        far_fraction = min(far_fraction, axis_far);
        if near_fraction > far_fraction return false;
    }
    return true;
}
