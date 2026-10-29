#import "Sort";

Triangle_Metadata :: struct {
    triangle_index: int;
    bounds: AABB;
    centroid: Vector3;
}

Triangle_BVH_Node :: struct {
    bounds: AABB;

    // this range indexes Triangle_BVH.references.
    first_triangle_reference_index: int;
    triangle_reference_count: int;
    
    // a node is either a leaf or has two children, this is enforced during construction

    // negative indices mean this node has no children and is a leaf.
    left_node_index: int = -1;
    right_node_index: int = -1;
}

MAXIMUM_TRIANGLES_PER_TRIANGLE_BVH_LEAF_NODE :: 8;

Triangle_BVH :: struct {
    // building reorders these references so every node owns a contiguous range.
    // triangle_index still addresses the caller's unchanged triangle array.
    references: [..] Triangle_Metadata;
    nodes: [..] Triangle_BVH_Node;
}

Triangle_Ray_Hit :: struct {
    // This is a world-space distance when direction is normalized.
    distance: float32;
    triangle_index: int = -1;

    // Weights for triangle a, b, and c, respectively.
    barycentric_coordinates: Vector3;
}

reference :: (
    triangle: $T/interface Triangle,
    triangle_index: int,
) -> Triangle_Metadata {
    return .{
        triangle_index = triangle_index,
        bounds = bounds(triangle),
        centroid = (triangle.a + triangle.b + triangle.c) / 3.0,
    };
}

reference_range_bounds :: (
    references: [] Triangle_Metadata,
    first: int,
    count: int,
) -> triangle_bounds: AABB, centroid_bounds: AABB {
    assert(count > 0 && first >= 0 && first + count <= references.count);

    triangle_bounds := MERGE_IDENTITY;
    centroid_bounds := MERGE_IDENTITY;
    for reference_index: first..first + count - 1 {
        reference := references[reference_index];
        merge(*triangle_bounds, reference.bounds);
        merge(*centroid_bounds, reference.centroid);
    }
    return triangle_bounds, centroid_bounds;
}

sort_reference_range_along_axis :: (
    bvh: *Triangle_BVH,
    node: Triangle_BVH_Node,
    axis: Axis,
) {
    reference_range := array_view(
        bvh.references,
        node.first_triangle_reference_index,
        node.triangle_reference_count,
    );

    if axis == {
        case .X; quick_sort(reference_range, reference => reference.centroid.x);
        case .Y; quick_sort(reference_range, reference => reference.centroid.y);
        case .Z; quick_sort(reference_range, reference => reference.centroid.z);
    }
}

try_subdivide :: (bvh: *Triangle_BVH, node_index: int) {
    // copy the node because appending children can reallocate bvh.nodes.
    node := bvh.nodes[node_index];
    if node.triangle_reference_count
        <= MAXIMUM_TRIANGLES_PER_TRIANGLE_BVH_LEAF_NODE return;

    _, centroid_bounds := reference_range_bounds(
        bvh.references,
        node.first_triangle_reference_index,
        node.triangle_reference_count,
    );
    split_axis := longest_axis(centroid_bounds);
    sort_reference_range_along_axis(bvh, node, split_axis);

    left_count := node.triangle_reference_count / 2;
    right_count := node.triangle_reference_count - left_count;
    left_first := node.first_triangle_reference_index;
    right_first := left_first + left_count;

    left_bounds, _ := reference_range_bounds(
        bvh.references, left_first, left_count);
    right_bounds, _ := reference_range_bounds(
        bvh.references, right_first, right_count);

    left_node_index := bvh.nodes.count;
    array_add(*bvh.nodes, .{
        bounds = left_bounds,
        first_triangle_reference_index = left_first,
        triangle_reference_count = left_count,
    });
    right_node_index := bvh.nodes.count;
    array_add(*bvh.nodes, .{
        bounds = right_bounds,
        first_triangle_reference_index = right_first,
        triangle_reference_count = right_count,
    });

    // look the parent up again because either append may have moved the array.
    bvh.nodes[node_index].left_node_index = left_node_index;
    bvh.nodes[node_index].right_node_index = right_node_index;

    try_subdivide(bvh, left_node_index);
    try_subdivide(bvh, right_node_index);
}

create :: (triangles: [] $T/interface Triangle) -> Triangle_BVH {
    result: Triangle_BVH;
    if triangles.count == 0 return result;

    array_reserve(*result.references, triangles.count);
    for triangle, triangle_index: triangles {
        array_add(*result.references,
            reference(triangle, triangle_index));
    }

    root_bounds, _ := reference_range_bounds(
        result.references, 0, result.references.count);
    array_add(*result.nodes, .{
        bounds = root_bounds,
        first_triangle_reference_index = 0,
        triangle_reference_count = result.references.count,
    });
    try_subdivide(*result, 0);
    return result;
}

node_is_leaf :: (node: Triangle_BVH_Node) -> bool {
    has_left_child := node.left_node_index >= 0;
    has_right_child := node.right_node_index >= 0;
    assert(has_left_child == has_right_child);
    // because of how we insert left child first, if left child doesn't exist neither does right
    return !has_left_child;
}

closest_ray_hit :: (
    origin: Vector3,
    direction: Vector3,
    tnear: float32,
    tfar: float32,
    triangles: [] $T/interface Triangle,
    bvh: *Triangle_BVH,
) -> found: bool, hit: Triangle_Ray_Hit {
    if bvh.nodes.count == 0 || tfar < tnear return false, .{};

    found := false;
    hit: Triangle_Ray_Hit;
    closest := tfar;
    stack: [128] int;
    stack_count := 1;
    stack[0] = 0;

    while stack_count > 0 {
        stack_count -= 1;
        node := bvh.nodes[stack[stack_count]];
        if !ray_segment_intersects_aabb(
            origin, direction, node.bounds.min, node.bounds.max, closest) continue;

        if node_is_leaf(node) {
            reference_end := node.first_triangle_reference_index
                + node.triangle_reference_count;
            for reference_index:
                node.first_triangle_reference_index..reference_end - 1 {
                triangle_index := bvh.references[reference_index].triangle_index;
                triangle := triangles[triangle_index];
                distance := ray_triangle_intersection(
                    origin, direction, triangle.a, triangle.b, triangle.c);
                if distance == FLOAT32_MAX
                    || distance < tnear || distance > closest continue;

                closest = distance;
                hit.distance = distance;
                hit.triangle_index = triangle_index;
                hit_position := origin + direction * distance;
                hit.barycentric_coordinates, _ = barycentric_coordinates_3d(
                    hit_position, triangle.a, triangle.b, triangle.c);
                found = true;
            }
        } else {
            assert(stack_count + 2 <= stack.count);
            stack[stack_count] = node.left_node_index;
            stack[stack_count + 1] = node.right_node_index;
            stack_count += 2;
        }
    }

    return found, hit;
}

deinit :: (using bvh: *Triangle_BVH) {
    array_free(references);
    array_free(nodes);
}

// queries

/*

returns true iff a triangle in the bvh intersects the aabb

node bounds provide the broad-phase rejection. Leaf triangles still require
the exact triangle-aabb test because overlapping bounds do not prove that the
enclosed geometry intersects the query.

The average case is log_2(n)

*/
aabb_intersects :: (
    bvh: *Triangle_BVH,
    triangles: [] $T/interface Triangle,
    aabb: AABB,
) -> bool {
    if bvh.nodes.count == 0 return false;
    return aabb_intersects(bvh, triangles, aabb, 0);
}

aabb_intersects :: (
    bvh: *Triangle_BVH,
    triangles: [] $T/interface Triangle,
    aabb: AABB,
    node_index: int,
) -> bool {
    node := bvh.nodes[node_index];
    if disjoint(node.bounds, aabb) return false;

    if node_is_leaf(node) {
        for offset: 0..node.triangle_reference_count - 1 {
            metadata_index := node.first_triangle_reference_index + offset;
            metadata := bvh.references[metadata_index];
            triangle := triangles[metadata.triangle_index];
            if triangle_intersects_aabb(triangle, aabb) return true;
        }
        return false;
    }
    
    // a node is either a leaf or has two children, this is enforced during construction

    return aabb_intersects(bvh, triangles, aabb, node.left_node_index) || 
    aabb_intersects(bvh, triangles, aabb, node.right_node_index);
}
