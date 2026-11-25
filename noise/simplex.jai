#import "Basic";
#import "Math";
/*
3D simplex noise following Stefan Gustavson, Simplex noise demystified (2005),
pages 9, 11-12: https://itn-web.it.liu.se/~stegu76/TNM084-2019/simplexnoise.pdf
Fixed reference permutation and twelve unnormalized cube-edge gradients; scale 32.
Pure, deterministic, allocation-free and safe for concurrent sampling. No seed or
texture is required. Finite coordinates with skewed cells inside s32 are expected;
float32 precision limits useful detail long before that bound. Permutation wraps
in lattice space, not independently along Cartesian input axes.
*/
simplex_noise :: (position: Vector3) -> float32 {
    skew: float32 = (position.x + position.y + position.z) / 3;
    i: s32 = simplex_floor(position.x + skew);
    j: s32 = simplex_floor(position.y + skew);
    k: s32 = simplex_floor(position.z + skew);
    unskew: float32 = (cast(float32)i + cast(float32)j + cast(float32)k) / 6;
    origin: Vector3 = position - Vector3.{cast(float32)i, cast(float32)j, cast(float32)k} + Vector3.{unskew, unskew, unskew};
    first: Vector3;
    second: Vector3;
    // Reference tie-breaking matters at equal coordinates.
    if origin.x >= origin.y {
        if origin.y >= origin.z { first = .{1,0,0}; second = .{1,1,0}; }
        else if origin.x >= origin.z { first = .{1,0,0}; second = .{1,0,1}; }
        else { first = .{0,0,1}; second = .{1,0,1}; }
    } else {
        if origin.y < origin.z { first = .{0,0,1}; second = .{0,1,1}; }
        else if origin.x < origin.z { first = .{0,1,0}; second = .{0,1,1}; }
        else { first = .{0,1,0}; second = .{1,1,0}; }
    }
    d1: Vector3 = origin - first + Vector3.{1.0/6,1.0/6,1.0/6};
    d2: Vector3 = origin - second + Vector3.{1.0/3,1.0/3,1.0/3};
    d3: Vector3 = origin - Vector3.{0.5,0.5,0.5};
    // Wrap cells before adding corner offsets, avoiding integer overflow.
    i &= 255; j &= 255; k &= 255;
    n0: float32 = simplex_corner(i, j, k, origin);
    n1: float32 = simplex_corner(i + cast(s32)first.x, j + cast(s32)first.y, k + cast(s32)first.z, d1);
    n2: float32 = simplex_corner(i + cast(s32)second.x, j + cast(s32)second.y, k + cast(s32)second.z, d2);
    n3: float32 = simplex_corner(i+1, j+1, k+1, d3);
    return 32 * (n0+n1+n2+n3);
}

#scope_module
simplex_floor :: (value: float32) -> s32 {
    integer: s32 = cast(s32)value;
    // Unlike the PDF's fastfloor, exact zero/negative integers use true floor.
    if cast(float32)integer > value integer -= 1;
    return integer;
}

simplex_permutation :: (index: s32) -> s32 {
    // A local fixed table also translates directly to a GLSL array constructor.
    // Native optimizing compilers keep this immutable data in read-only storage.
    permutation: [256] s32 = s32.[
151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,
140,36,103,30,69,142,8,99,37,240,21,10,23,190,6,148,
247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,
57,177,33,88,237,149,56,87,174,20,125,136,171,168,68,175,
74,165,71,134,139,48,27,166,77,146,158,231,83,111,229,122,
60,211,133,230,220,105,92,41,55,46,245,40,244,102,143,54,
65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,
200,196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,
52,217,226,250,124,123,5,202,38,147,118,126,255,82,85,212,
207,206,59,227,47,16,58,17,182,189,28,42,223,183,170,213,
119,248,152,2,44,154,163,70,221,153,101,155,167,43,172,9,
129,22,39,253,19,98,108,110,79,113,224,232,178,185,112,104,
218,246,97,228,251,34,242,193,238,210,144,12,191,179,162,241,
81,51,145,235,249,14,239,107,49,192,214,31,181,199,106,157,
184,84,204,176,115,121,50,45,127,4,150,254,138,236,205,93,
222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180];
    return permutation[index & 255];
}

simplex_corner :: (i: s32, j: s32, k: s32, displacement: Vector3) -> float32 {
    weight: float32 = 0.6 - displacement.x*displacement.x - displacement.y*displacement.y - displacement.z*displacement.z;
    if weight < 0 return 0;
    gradient: s32 = simplex_permutation(i + simplex_permutation(j + simplex_permutation(k))) % 12;
    contribution: float32 = 0;
    // Reference grad3 table, expressed as dot products to avoid a second array.
    if gradient == 0 contribution = displacement.x + displacement.y;
    if gradient == 1 contribution = -displacement.x + displacement.y;
    if gradient == 2 contribution = displacement.x - displacement.y;
    if gradient == 3 contribution = -displacement.x - displacement.y;
    if gradient == 4 contribution = displacement.x + displacement.z;
    if gradient == 5 contribution = -displacement.x + displacement.z;
    if gradient == 6 contribution = displacement.x - displacement.z;
    if gradient == 7 contribution = -displacement.x - displacement.z;
    if gradient == 8 contribution = displacement.y + displacement.z;
    if gradient == 9 contribution = -displacement.y + displacement.z;
    if gradient == 10 contribution = displacement.y - displacement.z;
    if gradient == 11 contribution = -displacement.y - displacement.z;
    weight *= weight;
    return weight * weight * contribution;
}

#scope_export
/*
Normalized octave composition of the reference 3D simplex primitive.
Returns signed noise in the nominal range [-1, 1], not unsigned [0, 1].
Normalization divides by the sum of layer weights; the result is not clamped.
Same count/persistence/lacunarity conventions as the existing Perlin fBm API.
A distinct API keeps existing fractal_brownian_motion callers on Perlin.
*/
simplex_fractal_brownian_motion :: (position: Vector3, octave_count: float32 = 4, persistence: float32 = 0.5, lacunarity: float32 = 2) -> float32 {
    octaves: float32 = octave_count;
    if octaves < 1 octaves = 1;
    if octaves > 16 octaves = 16;
    whole_octaves: s32 = cast(s32)octaves;
    fractional_octave: float32 = octaves - cast(float32)whole_octaves;
    amplitude_multiplier: float32 = persistence;
    if amplitude_multiplier < 0 amplitude_multiplier = 0;
    amplitude: float32 = 1;
    frequency: float32 = 1;
    total: float32 = 0;
    total_weight: float32 = 0;
    for octave: 0..whole_octaves-1 {
        total += simplex_noise(position * frequency) * amplitude;
        total_weight += amplitude;
        amplitude *= amplitude_multiplier;
        frequency *= lacunarity;
    }
    result: float32 = total / total_weight;
    if fractional_octave == 0 || amplitude == 0 return result;
    next_total: float32 = total + simplex_noise(position * frequency) * amplitude;
    return result + (next_total / (total_weight + amplitude) - result) * fractional_octave;
}
