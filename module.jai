#import "Math";

Axis_Plane :: enum u32 {
    XY;
    XZ;
    YZ;
}

Plane :: struct {
    normal: Vector3;
    // this is bad naming fix this.
    d: float;  // distance from origin
}

Triangle :: struct {
    p0, p1, p2: Vector3;
}

point_in_triangle :: (point: Vector3, tri: Triangle) -> bool {

    // barycentric coordinates
    edge_to_vertex2 := tri.p2 - tri.p0;
    edge_to_vertex1 := tri.p1 - tri.p0;
    edge_to_point := point - tri.p0;
    
    dot_edge2_edge2 := dot(edge_to_vertex2, edge_to_vertex2);
    dot_edge2_edge1 := dot(edge_to_vertex2, edge_to_vertex1);
    dot_edge2_point := dot(edge_to_vertex2, edge_to_point);
    dot_edge1_edge1 := dot(edge_to_vertex1, edge_to_vertex1);
    dot_edge1_point := dot(edge_to_vertex1, edge_to_point);
    
    inverse_denominator := 1.0 / (dot_edge2_edge2 * dot_edge1_edge1 - dot_edge2_edge1 * dot_edge2_edge1);
    barycentric_u := (dot_edge1_edge1 * dot_edge2_point - dot_edge2_edge1 * dot_edge1_point) * inverse_denominator;
    barycentric_v := (dot_edge2_edge2 * dot_edge1_point - dot_edge2_edge1 * dot_edge2_point) * inverse_denominator;
    
    return (barycentric_u >= 0) && (barycentric_v >= 0) && (barycentric_u + barycentric_v <= 1);
}

ZERO3 :: Vector3.{0, 0, 0};
ONE3 :: Vector3.{1, 1, 1};
X3 :: Vector3.{1, 0, 0};
Y3 :: Vector3.{0, 1, 0};
Z3 :: Vector3.{0, 0, 1};
X2 :: Vector2.{1, 0};
Y2 :: Vector2.{0, 1};

non_neg_mod :: (value : int, mod : int) -> int { 
    return (value % mod + mod) % mod;
}

componentwise_mul :: (a : Vector3, b: Vector3 ) -> Vector3 {
    return .{a.x * b.x, a.y * b.y, a.z * b.z};
}

normalized_cross :: (a: Vector3, b: Vector3, epsilon := 0.0, fallback :Vector3 = .{}) -> Vector3 { 
    return normalize(cross_product(a, b), epsilon, fallback);
}

// get lowest root of quadratic equation within range [0, max_r]
get_smallest_non_negative_root_of_quadratic :: (a: float, b: float, c: float, max_r: float) -> (found: bool, root: float) {
    determinant := b * b - 4.0 * a * c;
    
    if determinant < 0.0 { return false, 0; }
    
    sqrt_d := sqrt(determinant);
    r1 := (-b - sqrt_d) / (2.0 * a);
    r2 := (-b + sqrt_d) / (2.0 * a);
    
    // sort so r1 <= r2
    if r1 > r2 {
        temp := r1;
        r1 = r2;
        r2 = temp;
    }
    
    if r1 > 0 && r1 < max_r { return true, r1; }
    if r2 > 0 && r2 < max_r { return true, r2; }
    
    return false, 0;
}
