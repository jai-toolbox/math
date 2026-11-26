#import "Basic";
#import "Math";

#load "noise/perlin.jai";
#load "noise/simplex.jai";
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

linear_interpolate :: lerp;

/*
Unnormalized Gaussian with value 1 at zero: e^(-x^2*squared_input_scale).
A positive squared-input scale controls narrowing; the integral is not normalized to 1.
This multiplies x^2; a horizontal input scale c would correspond to c^2 here.
*/
gaussian :: inline (x: float32, squared_input_scale: float32 = 1) -> float32 {
    return exp(-x*x*squared_input_scale);
}

/*
Approximate blackbody chromaticity as linear sRGB, with brightness kept separate.
Tanner Helland's temperature-to-RGB fit (2012), decoded from sRGB to linear:
https://tannerhelland.com/2012/09/18/convert-temperature-rgb-algorithm-code.html
Finite positive inputs clamp to the fit's supported 1000..40000 kelvin range.
Nonpositive temperatures return black. This is a color approximation, not radiance.
*/
black_body_radiation_color :: (temperature_kelvin: float32) -> Vector3 {
    if temperature_kelvin <= 0 return .{0, 0, 0};
    temperature: float32 = clamp(temperature_kelvin, 1000.0, 40000.0) / 100;
    color: Vector3;
    if temperature <= 66 {
        color.x = 255;
        color.y = cast(float32) (99.4708025861 * log(temperature) - 161.1195681661);
    } else {
        color.x = cast(float32) (329.698727446 * pow(temperature - 60, -0.1332047592));
        color.y = cast(float32) (288.1221695283 * pow(temperature - 60, -0.0755148492));
    }
    if temperature >= 66 color.z = 255;
    else if temperature <= 19 color.z = 0;
    else color.z = cast(float32) (138.5177312231 * log(temperature - 10) - 305.0447927307);
    for *channel: color.component {
        encoded: float32 = clamp(channel.*, 0.0, 255.0) / 255;
        if encoded <= 0.04045 channel.* = encoded / 12.92;
        else channel.* = pow((encoded + 0.055) / 1.055, 2.4);
    }
    return color;
}

/*
Rotate a 32-bit unsigned value left, wrapping the high bits into the low bits.
Counts wrap modulo 32; masking both shifts also handles a count of zero.
*/
rotate_bits_left :: inline (value: u32, count: u32) -> u32 {
    shift := count & 31;
    return (value << shift) | (value >> ((32 - shift) & 31));
}

/*
Slab test for a finite ray segment against an axis-aligned bounding box.
The tested segment is origin + direction * t, where t is in
[0, end_fraction].
*/
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

/*
Evaluate a scalar cubic Bezier curve. a and d are endpoints; b and c are
control points. t runs from 0 to 1 along the curve and is not clamped.
*/
cubic_bezier :: (a: float32, b: float32, c: float32, d: float32, t: float32) -> float32 {
    u := 1 - t;
    return u*u*u*a + 3*u*u*t*b + 3*u*t*t*c + t*t*t*d;
}

/*
Convert a value's position in an ascending range to a fraction in [0, 1].
Values outside the range are clamped. Equal endpoints form a step:
values at or below the endpoint return 0, and values above it return 1.
*/
normalize_to_unit_interval :: (value: float32, start: float32, end: float32) -> float32 {
    if value <= start return 0;
    if value >= end return 1;
    return (value - start) / (end - start);
}

/*
Cubic interpolation over [0, 1] with smooth endpoints. It maps 0 to 0 and
1 to 1, and its first derivative is 0 at both ends, so a transition using it
starts and finishes without an abrupt change in slope.
Input is clamped to [0, 1].
*/
smoothstep :: (t: float32) -> float32 {
    t = clamp(t, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

/*
Quintic interpolation over [0, 1], also called smootherstep.
The unique polynomial of degree at most five satisfying
f(0) = 0, f(1) = 1, f'(0) = f'(1) = 0, and f''(0) = f''(1) = 0
is f(t) = 6t^5 - 15t^4 + 10t^3. These endpoint conditions make a
transition start and finish with zero slope and zero curvature.
The caller supplies t in [0, 1]; this procedure does not clamp it.
*/
smootherstep :: (t: float32) -> float32 {
    return t * t * t * (t * (t * 6 - 15) + 10);
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



/*
Moller-Trumbore ray/triangle intersection. The returned fraction is along
direction, which may be a finite segment; callers clamp it to their reach.
FLOAT32_MAX denotes a miss.
*/
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

/*
The weights correspond to a, b, and c. They are also defined outside the
triangle; inside is false there or when the triangle has zero area.
Proof: https://www.openmath.net/algebra/linear/vectors.html#proposition-barycentric-coordinates-from-dot-products
*/
barycentric_coordinates_3d :: (
    point: Vector3,
    a: Vector3,
    b: Vector3,
    c: Vector3,
) -> weights: Vector3, inside: bool {
    ab := b - a;
    ac := c - a;
    ap := point - a;

    ab_ab := dot(ab, ab);
    ab_ac := dot(ab, ac);
    ac_ac := dot(ac, ac);
    ap_ab := dot(ap, ab);
    ap_ac := dot(ap, ac);
    denominator := ab_ab * ac_ac - ab_ac * ab_ac;
    if abs(denominator) <= 0.000000000001 return .{}, false;

    weight_b := (ac_ac * ap_ab - ab_ac * ap_ac) / denominator;
    weight_c := (ab_ab * ap_ac - ab_ac * ap_ab) / denominator;
    weights := Vector3.{1.0 - weight_b - weight_c, weight_b, weight_c};

    epsilon :: 0.000001;
    inside := weights.x >= -epsilon && weights.y >= -epsilon && weights.z >= -epsilon;
    return weights, inside;
}

point_in_triangle :: (point: Vector3, tri: Triangle) -> bool {
    _, inside := barycentric_coordinates_3d(point, tri.a, tri.b, tri.c);
    return inside;
}

/*
The weights correspond to a, b, and c. They are also defined outside the
triangle; inside is false there or when the triangle has zero UV area.
*/
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
ZERO2 :: Vector2.{0, 0};

ONE3 :: Vector3.{1, 1, 1};
ONE2 :: Vector2.{1, 1};

NEGATIVE_ONE3 :: Vector3.{-1, -1, -1};

Vector3i :: struct {
    x, y, z: int;
}

// vector comparisons succeed only when every component satisfies the comparison.
operator < :: inline (a: Vector3i, b: Vector3i) -> bool {
    return a.x < b.x && a.y < b.y && a.z < b.z;
}

operator <= :: inline (a: Vector3i, b: Vector3i) -> bool {
    return a.x <= b.x && a.y <= b.y && a.z <= b.z;
}

operator > :: inline (a: Vector3i, b: Vector3i) -> bool {
    return a.x > b.x && a.y > b.y && a.z > b.z;
}

operator >= :: inline (a: Vector3i, b: Vector3i) -> bool {
    return a.x >= b.x && a.y >= b.y && a.z >= b.z;
}

operator < :: inline (a: Vector3i, b: int) -> bool {
    return a.x < b && a.y < b && a.z < b;
}

operator <= :: inline (a: Vector3i, b: int) -> bool {
    return a.x <= b && a.y <= b && a.z <= b;
}

operator > :: inline (a: Vector3i, b: int) -> bool {
    return a.x > b && a.y > b && a.z > b;
}

operator >= :: inline (a: Vector3i, b: int) -> bool {
    return a.x >= b && a.y >= b && a.z >= b;
}

X3 :: Vector3.{1, 0, 0};
Y3 :: Vector3.{0, 1, 0};
Z3 :: Vector3.{0, 0, 1};

project :: inline (vector_to_project: Vector3, onto_direction: Vector3) -> float {
    assert(onto_direction.x != 0 || onto_direction.y != 0 || onto_direction.z != 0, "Projection direction must be nonzero");
    return dot(vector_to_project, normalize(onto_direction));
}

project_onto_line :: inline (vector_to_project: Vector3, line_direction: Vector3) -> Vector3 {
    assert(line_direction.x != 0 || line_direction.y != 0 || line_direction.z != 0, "Projection direction must be nonzero");
    unit_direction := normalize(line_direction);
    return dot(vector_to_project, unit_direction) * unit_direction;
}

project_onto_plane :: inline (vector_to_project: Vector3, plane_normal: Vector3) -> Vector3 {
    return vector_to_project - project_onto_line(vector_to_project, plane_normal);
}

// Engine world convention: right-handed XYZ, +X forward, +Z up, -Y right.
// These describe world directions; asset-local axes and screen axes are separate.
WORLD_FORWARD  :: X3;
WORLD_BACKWARD :: #run -WORLD_FORWARD;
WORLD_RIGHT    :: #run -Y3;
WORLD_LEFT     :: #run -WORLD_RIGHT;
WORLD_UP       :: Z3;
WORLD_DOWN     :: #run -WORLD_UP;

// Conversions maintain handedness.
world_to_y_up :: inline (v: Vector3) -> Vector3 {
    return .{dot(v, WORLD_FORWARD), dot(v, WORLD_UP), dot(v, WORLD_RIGHT)};
}
world_from_y_up :: inline (v: Vector3) -> Vector3 {
    return v.x * WORLD_FORWARD + v.y * WORLD_UP + v.z * WORLD_RIGHT;
}

X2 :: Vector2.{1, 0};
Y2 :: Vector2.{0, 1};


// we inline these things to open them up for optimization

is_finite :: inline (v: Vector2) -> bool {
    return is_finite(v.x) && is_finite(v.y);
}

is_finite :: inline (v: Vector3) -> bool {
    return is_finite(v.x) && is_finite(v.y) && is_finite(v.z);
}

is_finite :: inline (v: Vector4) -> bool {
    return is_finite(v.x) && is_finite(v.y) && is_finite(v.z) && is_finite(v.w);
}

component_product :: inline (using v : Vector3) -> float {
    return x * y * z;
}

component_product :: inline (using v: Vector3i) -> int {
    return x * y * z;
}

pow :: inline (x: Vector3, power: float32) -> Vector3 {
    return .{pow(x.x, power), pow(x.y, power), pow(x.z, power)};
}

floor :: inline (using v : Vector3) -> Vector3 {
    return .{floor(x), floor(y), floor(z)};
}

floor :: inline (using v : *Vector3) {
    x = floor(x); y = floor(y); z = floor(z);
}

to_vector3i :: inline (v: Vector3) -> Vector3i {
    return .{cast(int) v.x, cast(int) v.y, cast(int) v.z};
}

to_vector3 :: inline (v: Vector3i) -> Vector3 {
    return .{cast(float32) v.x, cast(float32) v.y, cast(float32) v.z};
}

round :: inline (v: Vector3) -> Vector3i {
    return .{cast(int) round(v.x), cast(int) round(v.y), cast(int) round(v.z)};
}

clamp :: inline (value: Vector3, minimum: Vector3, maximum: Vector3) -> Vector3 {
    return .{
        clamp(value.x, minimum.x, maximum.x),
        clamp(value.y, minimum.y, maximum.y),
        clamp(value.z, minimum.z, maximum.z),
    };
}

clamp :: inline (value: Vector3, minimum: float32, maximum: float32) -> Vector3 {
    return .{
        clamp(value.x, minimum, maximum),
        clamp(value.y, minimum, maximum),
        clamp(value.z, minimum, maximum),
    };
}

clamp :: inline (value: Vector4, minimum: float32, maximum: float32) -> Vector4 {
    return .{
        clamp(value.x, minimum, maximum),
        clamp(value.y, minimum, maximum),
        clamp(value.z, minimum, maximum),
        clamp(value.w, minimum, maximum),
    };
}

/*

(0, 0), (1, 0)    ... (w-1, 0)
(0, 1), (1, 1)    ... (w-1, 1)
  .       .              .
  .       .              .
(0,h-1), (1,h-1)  ... (w-1,h-1)

there is a canonical bijection between this and {0, ..., w * h - 1}

  0      1             w-1      w                   w + (w-1)      h*w                   h*w + (w-1)
(0, 0), (1, 0),  ... (w-1, 0) (0, 1), (1, 1),  ... (w-1, 1) ...  (0,h-1), (1,h-1),  ... (w-1,h-1)

so that

f(x, y) = w * y + x

*/


flatten_index :: inline (coordinate: Vector2, dimensions: Vector2) -> int {
    x := cast(int) coordinate.x;
    y := cast(int) coordinate.y;
    width := cast(int) dimensions.x;
    return y * width + x;
}

// the same idea works for 3d but here we iterate over all the x's, then iterate over the y's and then the z's

flatten_index :: inline (
    coordinate: Vector3i,
    dimensions: Vector3i,
) -> int {
    return (coordinate.z * dimensions.y + coordinate.y)
        * dimensions.x + coordinate.x;
}


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

round :: inline (x: float) -> s64 {
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

/*
Shortest rotation taking from's direction onto to's direction. Input lengths
need not match; both vectors must be finite and nonzero. Opposite directions
have no unique rotation axis, so choose a deterministic perpendicular one.
*/
rotation_between_vectors :: (from: Vector3, to: Vector3) -> Quaternion {
    assert(is_finite(from) && is_finite(to), "Rotation directions must be finite");
    
    from_scale := max(abs(from.x), max(abs(from.y), abs(from.z)));
    to_scale := max(abs(to.x), max(abs(to.y), abs(to.z)));
    
    assert(from_scale > 0 && to_scale > 0, "Rotation directions must be nonzero");
    
    // Scaling first keeps normalization valid for very large/small vectors.
    a := normalize(from / from_scale);
    b := normalize(to / to_scale);
    alignment := clamp(dot(a, b), -1.0, 1.0);
    axis := cross(a, b);
    if alignment >= 0 return normalize(Quaternion.{axis.x, axis.y, axis.z, 1 + alignment});

    // Cancellation near opposite directions can leave a small component along
    // a. So we remove it before normalizing the very short cross product.
    axis -= a * dot(a, axis);
    axis_length := length(axis);
    if axis_length > 0 {
        // equivalent half-angle construction, stable near 180 degrees where
        // 1 + alignment can round to zero. 
        // The vector part is a unit axis.
        axis /= axis_length;
        return normalize(Quaternion.{axis.x, axis.y, axis.z, axis_length / (1 - alignment)});
    }

    basis := Z3;
    least_alignment := abs(a.z);
    if abs(a.y) < least_alignment { basis = Y3; least_alignment = abs(a.y); }
    if abs(a.x) < least_alignment basis = X3;
    axis = normalize(cross(a, basis));
    return .{axis.x, axis.y, axis.z, 0};
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



/*
Linearly map between ranges, including reversed endpoints. Without clamping,
values outside the input range extrapolate. Equal input endpoints produce
zero, restricted to the output range when clamp_result is enabled.
*/
map_range :: (value: float32, from_min: float32, from_max: float32, to_min: float32, to_max: float32, clamp_result: bool) -> float32 {
    result: float32 = 0;
    if from_min != from_max {
        factor: float32 = (value - from_min) / (from_max - from_min);
        result = to_min + factor * (to_max - to_min);
    }
    if clamp_result { result = clamp(result, min(to_min, to_max), max(to_min, to_max)); }
    return result;
}
