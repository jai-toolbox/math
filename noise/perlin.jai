#import "Basic";
#import "Math";
#scope_module

// Fixed Fisher-Yates permutation generated with LCG seed 0x5045524c.
// Read-only data makes sampling deterministic and safe across bake workers.
PERLIN_PERMUTATION :: s32.[
    66, 2, 8, 242, 72, 1, 150, 96, 6, 14, 140, 18, 126, 161, 182, 156,
    176, 136, 20, 135, 60, 53, 98, 33, 110, 84, 181, 179, 146, 170, 7, 80,
    148, 34, 36, 50, 208, 153, 92, 83, 252, 185, 188, 166, 206, 220, 62, 55,
    42, 194, 214, 167, 52, 29, 0, 100, 5, 216, 236, 149, 184, 41, 40, 77,
    183, 111, 4, 247, 122, 248, 232, 97, 54, 226, 65, 217, 28, 212, 22, 114,
    88, 128, 125, 177, 26, 9, 16, 230, 32, 89, 241, 99, 70, 103, 108, 127,
    121, 235, 244, 129, 74, 172, 37, 51, 64, 253, 78, 87, 254, 171, 10, 231,
    222, 213, 211, 130, 209, 120, 246, 223, 152, 124, 192, 137, 250, 143, 106, 159,
    224, 218, 204, 71, 58, 101, 186, 189, 190, 151, 145, 196, 102, 119, 200, 210,
    107, 17, 38, 67, 160, 207, 164, 197, 162, 69, 82, 35, 118, 81, 165, 240,
    154, 221, 13, 11, 24, 178, 195, 43, 68, 113, 105, 47, 169, 45, 141, 227,
    56, 19, 109, 23, 116, 25, 132, 147, 138, 203, 48, 15, 134, 201, 30, 199,
    76, 191, 174, 219, 175, 75, 104, 255, 168, 3, 46, 215, 205, 131, 112, 187,
    173, 133, 144, 95, 180, 245, 57, 91, 90, 157, 158, 238, 198, 21, 44, 228,
    93, 155, 94, 251, 225, 49, 139, 234, 63, 123, 233, 79, 86, 237, 39, 27,
    249, 31, 12, 163, 202, 115, 229, 243, 73, 193, 61, 239, 142, 85, 117, 59,
];

#scope_export
/*
Improved Perlin noise in two, three or four dimensions. The vector type selects
an implementation specialized at compile time, sampling 4, 8 or 16 cell corners.
Gradient directions have unit length: eight compass directions in 2D, cube-edge
directions in 3D, and permutations of (0, +/-1, +/-1, +/-1) in 4D.
The signed result is in [-1, 1]; integer lattice points return zero. The fixed
permutation repeats every 256 units per axis. Coordinates must be finite;
float32 precision limits detail at large values.
*/
perlin_noise :: (position: Vector2) -> float32 {
    return sample_perlin_noise(2, position.component);
}

perlin_noise :: (position: Vector3) -> float32 {
    return sample_perlin_noise(3, position.component);
}

perlin_noise :: (position: Vector4) -> float32 {
    return sample_perlin_noise(4, position.component);
}

#scope_module
sample_perlin_noise :: ($dimensions: s32, coordinates: [dimensions] float32) -> float32 {
    #assert dimensions >= 2 && dimensions <= 4;
    CORNER_COUNT :: 1 << dimensions;
    cell: [dimensions] s32;
    fraction: [dimensions] float32;
    blend: [dimensions] float32;
    for axis: 0..dimensions-1 {
        coordinate := coordinates[axis];
        lower := floor(coordinate);
        // Wrap before converting to an integer, including for negative cells.
        cell[axis] = cast(s32)(lower - 256 * floor(lower / 256));
        fraction[axis] = coordinate - lower;
        blend[axis] = smootherstep(fraction[axis]);
    }

    corner_values: [CORNER_COUNT] float32;
    for corner: 0..CORNER_COUNT-1 {
        hash: s32;
        displacement: [dimensions] float32;
        for axis: 0..dimensions-1 {
            offset := (corner >> axis) & 1;
            hash = PERLIN_PERMUTATION[(hash + cell[axis] + offset) & 255];
            displacement[axis] = fraction[axis] - cast(float32)offset;
        }

        #if dimensions == 2 {
            direction := hash & 7;
            if direction < 4 {
                x := displacement[0]; y := displacement[1];
                if (direction & 1) != 0 x = -x;
                if (direction & 2) != 0 y = -y;
                corner_values[corner] = (x + y) / sqrt(2.0);
            } else {
                component := displacement[(direction >> 1) & 1];
                if (direction & 1) != 0 component = -component;
                corner_values[corner] = component;
            }
        } else #if dimensions == 3 {
            // Improved Perlin's 16 hash cases select the 12 cube-edge directions.
            direction := hash & 15;
            first := ifx direction < 8 then displacement[0] else displacement[1];
            second := displacement[2];
            if direction < 4 second = displacement[1];
            else if direction == 12 || direction == 14 second = displacement[0];
            if (direction & 1) != 0 first = -first;
            if (direction & 2) != 0 second = -second;
            corner_values[corner] = first + second;
        } else {
            omitted_axis := (hash & 31) >> 3;
            sign_index: s32;
            gradient_dot: float32;
            for axis: 0..dimensions-1 {
                if axis == omitted_axis continue;
                component := displacement[axis];
                if (hash & (1 << sign_index)) != 0 component = -component;
                gradient_dot += component;
                sign_index += 1;
            }
            corner_values[corner] = gradient_dot;
        }
    }

    remaining := CORNER_COUNT;
    for axis: 0..dimensions-1 {
        remaining /= 2;
        for i: 0..remaining-1 {
            corner_values[i] = linear_interpolate(corner_values[2*i], corner_values[2*i+1], blend[axis]);
        }
    }
    #if dimensions == 2 return corner_values[0];
    else return corner_values[0] / sqrt(cast(float32)(dimensions - 1));
}

#scope_export
/*
Normalized fractal Brownian motion built from Perlin noise. Each octave scales
frequency by lacunarity and amplitude by persistence. The result is a weighted
average in [-1, 1], with the dimension selected by the vector type at compile time.
Octave counts are clamped to [1, 16]. A fractional count crossfades the normalized
results for the two adjacent whole counts. Negative persistence is treated as
zero; finite parameters and positive lacunarity are expected.
*/
fractal_brownian_motion :: (position: Vector2, octave_count: float32 = 4, persistence: float32 = 0.5, lacunarity: float32 = 2) -> float32 {
    return sample_fractal_brownian_motion(position, octave_count, persistence, lacunarity);
}

fractal_brownian_motion :: (position: Vector3, octave_count: float32 = 4, persistence: float32 = 0.5, lacunarity: float32 = 2) -> float32 {
    return sample_fractal_brownian_motion(position, octave_count, persistence, lacunarity);
}

fractal_brownian_motion :: (position: Vector4, octave_count: float32 = 4, persistence: float32 = 0.5, lacunarity: float32 = 2) -> float32 {
    return sample_fractal_brownian_motion(position, octave_count, persistence, lacunarity);
}

#scope_module
sample_fractal_brownian_motion :: (position: $Vector, octave_count: float32, persistence: float32, lacunarity: float32) -> float32 {
    octaves := clamp(octave_count, 1.0, 16.0);
    whole_octaves := cast(s32)octaves;
    fractional_octave := octaves - cast(float32)whole_octaves;
    amplitude: float32 = 1;
    frequency: float32 = 1;
    total: float32;
    total_weight: float32;
    amplitude_multiplier := max(persistence, 0.0);
    for octave: 0..whole_octaves-1 {
        total += perlin_noise(position * frequency) * amplitude;
        total_weight += amplitude;
        amplitude *= amplitude_multiplier;
        frequency *= lacunarity;
    }
    result := total / total_weight;
    if fractional_octave == 0 || amplitude == 0 return result;
    next_total := total + perlin_noise(position * frequency) * amplitude;
    return linear_interpolate(result, next_total / (total_weight + amplitude), fractional_octave);
}
