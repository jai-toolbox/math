#import "Basic";
#import "Math";

#load "axis_aligned_bounding_box.jai";
#load "triangle.jai";
#load "triangle_bounding_volume_hierarchy.jai";

Axis :: enum {
    X;
    Y;
    Z;
}

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

// Cubic interpolation over [0, 1] with smooth endpoints. It maps 0 to 0 and
// 1 to 1, and its first derivative is 0 at both ends, so a transition using it
// starts and finishes without an abrupt change in slope.
smooth_step :: (t: float32) -> float32 {
    return t * t * (3.0 - 2.0 * t);
}

// note this is integer based.
log_2 :: (n : int) -> int {
    assert(n > 0);
    result := 0;
    // 00100101...0110 a non-zero generic binary number
    // always reduces to 1 after some amount of bitshifts.
    while n != 1 {
        result += 1;
        n >>= 1;
    }
    return result;
}



// Moller-Trumbore ray/triangle intersection. The returned fraction is along
// direction, which may be a finite segment; callers clamp it to their reach.
// FLOAT32_MAX denotes a miss.
ray_triangle_intersection :: (origin: Vector3, direction: Vector3, v0: Vector3, v1: Vector3, v2: Vector3) -> float32 {
    epsilon: float32 = 1.0e-12;
    e1 := v1 - v0;
    e2 := v2 - v0;
    p := cross(direction, e2);
    det := dot(e1, p);
    det_near_zero := abs(det) < epsilon;
    if det_near_zero det = 1.0;

    s := origin - v0;
    u := dot(s, p) / det;
    q := cross(s, e1);
    v := dot(direction, q) / det;
    t := dot(e2, q) / det;

    if det_near_zero || u < 0.0 || v < 0.0 || u + v > 1.0 || t < 0.0 return FLOAT32_MAX;
    return t;
}

// proof: https://www.openmath.net/algebra/linear/vectors.html#proposition-barycentric-coordinates-from-dot-products
point_in_triangle :: (point: Vector3, tri: Triangle) -> bool {

    // barycentric coordinates
    edge_to_vertex2 := tri.c - tri.a;
    edge_to_vertex1 := tri.b - tri.a;
    edge_to_point := point - tri.a;
    
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

// The weights correspond to a, b, and c. They are also defined outside the
// triangle; inside is false there or when the triangle has zero UV area.
barycentric_coordinates_2d :: (point: Vector2, a: Vector2, b: Vector2, c: Vector2) -> weights: Vector3, inside: bool {
    ab := b - a;
    ac := c - a;
    ap := point - a;
    denominator := ab.x * ac.y - ab.y * ac.x;
    if abs(denominator) <= 0.000000000001 return .{}, false;

    weight_b := (ap.x * ac.y - ap.y * ac.x) / denominator;
    weight_c := (ab.x * ap.y - ab.y * ap.x) / denominator;
    weights := Vector3.{1 - weight_b - weight_c, weight_b, weight_c};

    epsilon :: 0.000001;
    inside := weights.x >= -epsilon && weights.y >= -epsilon && weights.z >= -epsilon;
    return weights, inside;
}

point_in_triangle :: (point: Vector2, a: Vector2, b: Vector2, c: Vector2) -> bool {
    _, inside := barycentric_coordinates_2d(point, a, b, c);
    return inside;
}

ZERO3 :: Vector3.{0, 0, 0};
ONE3 :: Vector3.{1, 1, 1};
NEGATIVE_ONE3 :: Vector3.{-1, -1, -1};

X3 :: Vector3.{1, 0, 0};
Y3 :: Vector3.{0, 1, 0};
Z3 :: Vector3.{0, 0, 1};
X2 :: Vector2.{1, 0};
Y2 :: Vector2.{0, 1};

/*
Branchfree saturating addition for unsigned 32-bit integers.

Recall this fact from number theory:
https://www.openmath.net/number_theory/quotient_and_remainder.html#proposition-remainder-of-sum-less-than-first-summand

When adding two unsigned 32 bit integers if overflow happens, it is equvalent to just doing % (U32_MAX + 1), thus the mapping from code to math is valid.

Let a, b be two unsigned 32-bit integers, so a, b in [0, ..., U32_MAX], and use n = U32_MAX + 1.

If the addition overflows, res wraps around and becomes less than a by the proposition.

That means that when overflow occurs then res < a evalutes to 1, 0 otherwise.

-1 in twos complement is 0xFFFFFFFF
-0 in twos complement is 0x00000000

if you bitwise or 0xFFFFFFFF with any number you get 0xFFFFFFFF back.
if you bitwise or 0x00000000 with any number you get back the original.

This is how we "bake" the conditional logic into one thing.

 */
saturating_add :: (a: u32, b: u32) -> u32 {
    res := a + b;

    res |= cast,no_check(u32)(-(cast(s32)(res < a)));
    return res;
}

/*
Branchfree saturating subtraction for unsigned 32-bit integers.

If the subtraction underflows, res wraps around and becomes greater than a.
We exploit this: (res <= a) evaluates to 1 when no underflow, 0 on underflow.
Negating that: -(1) in two's complement is 0xFFFFFFFF and -(0) is 0x00000000.

So on underflow, res &= 0x00000000 clamps to 0.
On no underflow, res &= 0xFFFFFFFF doesn't change res.
 */
saturating_sub :: (a: u32, b: u32) -> u32 {
    res := a - b;
    res &= cast,no_check(u32)(-(cast(s32)(res <= a)));
    return res;
}

round :: (x: float) -> s64 {
    if x >= 0 {
        return cast(s64)(x + 0.5);
    } else {
        return cast(s64)(x - 0.5);
    }
}


non_neg_mod :: (value : int, mod : int) -> int { 
    return (value % mod + mod) % mod;
}

componentwise_mul :: (a : Vector3, b: Vector3 ) -> Vector3 {
    return .{a.x * b.x, a.y * b.y, a.z * b.z};
}

normalized_cross :: (a: Vector3, b: Vector3, epsilon := 0.0, fallback :Vector3 = .{}) -> Vector3 { 
    return normalize(cross_product(a, b), epsilon, fallback);
}

cos_of_angle_between :: (a: Vector3, b: Vector3) -> float {
    return dot(normalize(a), normalize(b));
}

angle_between :: (a: Vector3, b: Vector3) -> float {
    return acos(clamp(cos_of_angle_between(a, b), -1.0, 1.0));
}

DIRECTION_EPSILON :: 0.00001;

points_same_way :: (a: Vector3, b: Vector3) -> bool {
    if length_squared(a) <= DIRECTION_EPSILON return false;
    if length_squared(b) <= DIRECTION_EPSILON return false;

    an := normalize(a);
    bn := normalize(b);

    return dot(an, bn) > DIRECTION_EPSILON;
}

points_opposite_way :: (a: Vector3, b: Vector3) -> bool {
    if length_squared(a) <= DIRECTION_EPSILON return false;
    if length_squared(b) <= DIRECTION_EPSILON return false;

    an := normalize(a);
    bn := normalize(b);

    return dot(an, bn) < -DIRECTION_EPSILON;
}

is_perpendicular :: (a: Vector3, b: Vector3) -> bool {
    if length_squared(a) <= DIRECTION_EPSILON return false;
    if length_squared(b) <= DIRECTION_EPSILON return false;

    an := normalize(a);
    bn := normalize(b);

    return abs(dot(an, bn)) <= DIRECTION_EPSILON;
}

// get lowest root of quadratic equation within range [0, max_r]
get_smallest_non_negative_root_of_quadratic :: (a: float, b: float, c: float, max_r: float) -> (found: bool, root: float) {
    EPSILON :: 0.000001;

    // degenerate into a linear equation when the quadratic term is effectively zero, we solve b*t + c = 0 directly.
    if abs(a) <= EPSILON {
        if abs(b) <= EPSILON return false, 0;


        r := -c / b;
        if r > 0 && r < max_r return true, r;
        return false, 0;
    }

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
