#import "Math";

Axis_Plane :: enum u32 {
    XY;
    XZ;
    YZ;
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
