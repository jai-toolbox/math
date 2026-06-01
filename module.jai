#import "Math";

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
