struct Surface {
    float sd; // signed distance value
    vec3 col; // color
};

Surface sdSphere(vec3 p, float s, vec3 col) {
    float d = length(p) - s;
    return Surface(d, col);
}

//rotate around point
vec3 rotate(vec3 p, vec3 axis, float angle) {
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    return vec3(
        p.x * (axis.x * axis.x * oc + c) +
        p.y * (axis.x * axis.y * oc - axis.z * s) +
        p.z * (axis.x * axis.z * oc + axis.y * s),
        p.x * (axis.y * axis.x * oc + axis.z * s) +
        p.y * (axis.y * axis.y * oc + c) +
        p.z * (axis.y * axis.z * oc - axis.x * s),
        p.x * (axis.z * axis.x * oc - axis.y * s) +
        p.y * (axis.z * axis.y * oc + axis.x * s) +
        p.z * (axis.z * axis.z * oc + c)
    );
}


//Noise

vec3 mod289(vec3 x) {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 mod289(vec4 x) {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x) {
    return mod289(((x * 34.0) + 1.0) * x);
}

vec4 taylorInvSqrt(vec4 r) {
    return 1.79284291400159 - 0.85373472095314 * r;
}

float simplexNoise(vec3 v) {
    const vec2 C = vec2(1.0 / 6.0, 1.0 / 3.0);
    const vec4 D = vec4(0.0, 0.5, 1.0, 2.0);

    // First corner
    vec3 i = floor(v + dot(v, C.yyy));
    vec3 x0 = v - i + dot(i, C.xxx);

    // Other corners
    vec3 g = step(x0.yzx, x0.xyz);
    vec3 l = 1.0 - g;
    vec3 i1 = min(g.xyz, l.zxy);
    vec3 i2 = max(g.xyz, l.zxy);

    // x1, x2, x3
    vec3 x1 = x0 - i1 + C.xxx;
    vec3 x2 = x0 - i2 + C.yyy;
    vec3 x3 = x0 - D.yyy;

     // Permutations
    i = mod289(i);
    vec4 p1 = permute(i.z + vec4(0.0, i1.z, i2.z, 1.0));
    vec4 p2 = permute(p1 + i.y + vec4(0.0, i1.y, i2.y, 1.0));
    vec4 p = permute(p2 + i.x + vec4(0.0, i1.x, i2.x, 1.0));

    // Gradients
    vec4 t = vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3));
    vec4 norm = taylorInvSqrt(t);
    float grad1 = dot(x0, p.xyz);
    float grad2 = dot(x1, p.yzw);
    float grad3 = dot(x2, p.zww);
    float grad4 = dot(x3, p.www);
    vec4 grad_combine = vec4(grad1, grad2, grad3, grad4);
    vec4 grads = 70.0 * grad_combine * norm;

    // Perlin noise
    float noiseValue = dot(grads, vec4(x0.x, x1.x, x2.x, x3.x)); // Take only x component for each vec3

    return 2.2 * noiseValue;
}


float hash(float n) {
    return fract(sin(n) * 43758.5453);
}

float noise(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = smoothstep(0.0, 1., f);
    float n = p.x + p.y * 157.0 + 113.0 * p.z;
    return mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
                   mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
               mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
                   mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
}


Surface apply_octave_noise(in Surface s, vec3 x, int octaves, float factor) {
    float n = 0.;
    float a = factor < 1. ? 1. - factor : 1.;
    for(int i = 0; i < octaves; i++) {
        n += a * noise(x);
        x *= 2.;
        a *= 0.5;
    }
    s.sd -= factor * n;
    return s;
}

float apply_simplex_noise(in Surface s, vec3 pos, float scale) {
    // Apply simplex noise to the position
    float noiseValue = simplexNoise(pos);

    // Modify the signed distance value using noise
    float modifiedSd = s.sd - scale * noiseValue;

    // Return the modified signed distance value
    return modifiedSd;
}

Surface apply_sky(in Surface s, vec3 x, float scale) {
    s.sd = apply_simplex_noise(s, x, scale);
    return s;
}

float mod289(float x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
//vec4 mod289(vec4 x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
vec4 perm(vec4 x){return mod289(((x * 34.0) + 1.0) * x);}

float noises(vec3 p){
    vec3 a = floor(p);
    vec3 d = p - a;
    d = d * d * (3.0 - 2.0 * d);

    vec4 b = a.xxyy + vec4(0.0, 1.0, 0.0, 1.0);
    vec4 k1 = perm(b.xyxy);
    vec4 k2 = perm(k1.xyxy + b.zzww);

    vec4 c = k2 + a.zzzz;
    vec4 k3 = perm(c);
    vec4 k4 = perm(c + 1.0);

    vec4 o1 = fract(k3 * (1.0 / 41.0));
    vec4 o2 = fract(k4 * (1.0 / 41.0));

    vec4 o3 = o2 * d.z + o1 * (1.0 - d.z);
    vec2 o4 = o3.yw * d.x + o3.xz * (1.0 - d.x);

    return o4.y * d.y + o4.x * (1.0 - d.y);
}

vec3 apply_coloring(in Surface s) {
        // Coloring depending on the distance
        const vec3 high_color = vec3(0.50, 0.35, 0.15);
        const float high_level = 0.8;
        const vec3 mid_color = vec3(0.05, 1.15, 0.10);
        const float mid_level = 0.5;
        const vec3 low_color = vec3(0.0, 0.0, 1.0);
        const float low_level = 0.0;

        // Coloring
        float d = s.sd * 1.2;

        if (d >= high_level) {
            return high_color;
        } else if (d >= mid_level) {
            return mix(mid_color, high_color, (d - mid_level) / (high_level - mid_level));
        } else if (d >= low_level) {
            return mix(low_color, mid_color, (d - low_level) / (mid_level - low_level));
        } else {
            return low_color;
        }
    }

//---------------------------------
//includes Coloring
Surface map(in vec3 pos) {
    vec3 k = pos;
    Surface s = sdSphere(k, 0.5, vec3(0.));

    // Add some noise
    Surface s1 = apply_octave_noise(s, k, 4, .5);
    Surface s2 = apply_octave_noise(s1, k * 2.0, 2, .5);
    s2.col = apply_coloring(s);

    // Surface sky = apply_sky(sdSphere(k, 0.8, vec3(0.)), k, 1.);
    // sky = apply_sky(sky, k*1.1, .2);
    // sky.col = apply_coloring(sky);

    return s2;
}

vec3 calcNormal(in vec3 pos) {
    const float ep = 0.0001;
    vec2 e = vec2(1.0, -1.0) * 0.5773;
    return normalize(e.xyy * (map(pos + e.xyy * ep).sd) +
        e.yyx * (map(pos + e.yyx * ep).sd) +
        e.yxy * (map(pos + e.yxy * ep).sd) +
        e.xxx * (map(pos + e.xxx * ep).sd));
}

float calcSoftshadow(in vec3 ro, in vec3 rd, float tmin, float tmax, const float k) {
    float res = 1.0;
    float t = tmin;
    for(int i = 0; i < 50; i++) {
        Surface s = map(ro + rd * t);
        float h = s.sd;
        res = min(res, k * h / t);
        t += clamp(h, 0.02, 0.20);
        if(res < 0.005 || t > tmax)
            break;
    }
    return clamp(res, 0.0, 1.0);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 tot = vec3(0.0);

   // 1) Set up the Camera (primary ray for this pixel) 
    vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;

    vec3 ro = vec3(0., 3.0, 10.4);
    vec3 rd = normalize(vec3(p - vec2(0.1, 1.9), -5));

    Surface s; 
   // 2) Raymarching 
    float t = 0.1;
    for(int i = 0; i < 256; i++) {
        vec3 p = ro + t * rd;
        s = map(p);
        float h = s.sd;
        if(abs(h) < 0.0001 || t > 50.0)
            break;
        t += h;
    }

    vec3 col = vec3(0.0);

   // 3) Shading
    if(t < 50.0) {
        vec3 pos = ro + t * rd;
        vec3 nor = calcNormal(pos);
        vec3 lig = normalize(vec3(1.0, 1., 0.2));
        lig = rotate(lig, vec3(0.0, 1.0, 0.5), 1.5 * iTime);
        float dif = clamp(dot(nor, lig), 0.0, 1.0);
        float sha = calcSoftshadow(pos, lig, 0.001, 1.0, 16.0);
        float amb = 0.5 + 0.5 * nor.y;
        col = vec3(0.05, 0.1, 0.15) * amb + vec3(1.00, 0.9, 0.80) * dif * sha * s.col;
    }

    col = sqrt(col);
    tot += col;

    fragColor = vec4(tot, 1.0);
}