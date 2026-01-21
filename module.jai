#import "Math";

Axis_Plane :: enum u32 {
    XY;
    XZ;
    YZ;
}

x3 := Vector3.{1, 0, 0};
y3 := Vector3.{0, 1, 0};
z3 := Vector3.{0, 0, 1};
x2 := Vector2.{1, 0};
y2 := Vector2.{0, 1};

non_neg_mod :: (value : int, mod : int) -> int { 
    return (value % mod + mod) % mod;
}


normalized_cross :: (a: Vector3, b: Vector3, epsilon := 0.0, fallback :Vector3 = .{}) -> Vector3 { 
    return normalize(cross_product(a, b), epsilon, fallback);
}
