struct Surface
{
    float sd; // signed distance value
    vec3 col; // color
};

struct Noise
{
    float factor; // factor to multiply the noise by   
    int octaves; // number of octaves that represent how many times the noise is repeated
    float lacunarity; // factor that determines how much the frequency of the noise increases with each octave
    float gain; // factor that determines how much the amplitude of the noise decreases with each octave
};

float opSmoothUnion(float d1, float d2, float k)
{
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

Surface sdSphere(vec3 p, float s, vec3 col)
{
    float d = length(p) - s;
    return Surface(d, col);
}

float rand(float n)
{
    return fract(sin(n) * 43758.5453123);
}

vec3 get_random_star_color(float n, vec3 k)
{
    float r = rand(n + 1.0);
    float g = rand(n + 2.0);
    float b = rand(n + 3.0);

    return vec3(r, g, b);
}

vec3 apply_coloring_by_height(in Surface s)
{
        // Coloring depending on the distance
    const vec3 high_color = vec3(0.39, 0.17, 0.04);
    float high_level = 0.2;
    const vec3 mid_color = vec3(0.13, 0.52, 0.11);
    float mid_level = 0.1;
    const vec3 low_color = vec3(0.04, 0.16, 0.53);
    float low_level = 0.0;

        // Coloring
    float d = s.sd * 1.2;

    if(d >= high_level)
    {
        return high_color;
    }
    else if(d >= mid_level && d < high_level)
    {
        return mix(mid_color, high_color, (d - mid_level) / (high_level - mid_level));
    }
    else if(d >= low_level && d < mid_level)
    {
        return mix(low_color, mid_color, (d - low_level) / (mid_level - low_level));
    }
    else if(d < low_level)
    {
        return low_color;
    }
}

//rotate around point
vec3 rotate(vec3 p, vec3 axis, float angle)
{
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    return vec3(p.x * (axis.x * axis.x * oc + c) +
        p.y * (axis.x * axis.y * oc - axis.z * s) +
        p.z * (axis.x * axis.z * oc + axis.y * s), p.x * (axis.y * axis.x * oc + axis.z * s) +
        p.y * (axis.y * axis.y * oc + c) +
        p.z * (axis.y * axis.z * oc - axis.x * s), p.x * (axis.z * axis.x * oc - axis.y * s) +
        p.y * (axis.z * axis.y * oc + axis.x * s) +
        p.z * (axis.z * axis.z * oc + c));
}

// Simplex noise
vec3 mod289(vec3 x)
{
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 mod289(vec4 x)
{
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x)
{
    return mod289(((x * 34.0) + 1.0) * x);
}

vec4 taylorInvSqrt(vec4 r)
{
    return 1.79284291400159 - 0.85373472095314 * r;
}

float snoise(vec3 v)
{
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

  //   x0 = x0 - 0.0 + 0.0 * C.xxx;
  //   x1 = x0 - i1  + 1.0 * C.xxx;
  //   x2 = x0 - i2  + 2.0 * C.xxx;
  //   x3 = x0 - 1.0 + 3.0 * C.xxx;
    vec3 x1 = x0 - i1 + C.xxx;
    vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
    vec3 x3 = x0 - D.yyy;      // -1.0+3.0*C.x = -0.5 = -D.y

// Permutations
    i = mod289(i);
    vec4 p = permute(permute(permute(i.z + vec4(0.0, i1.z, i2.z, 1.0)) + i.y + vec4(0.0, i1.y, i2.y, 1.0)) + i.x + vec4(0.0, i1.x, i2.x, 1.0));

// Gradients: 7x7 points over a square, mapped onto an octahedron.
// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
    float n_ = 0.142857142857; // 1.0/7.0
    vec3 ns = n_ * D.wyz - D.xzx;

    vec4 j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)

    vec4 x_ = floor(j * ns.z);
    vec4 y_ = floor(j - 7.0 * x_);    // mod(j,N)

    vec4 x = x_ * ns.x + ns.yyyy;
    vec4 y = y_ * ns.x + ns.yyyy;
    vec4 h = 1.0 - abs(x) - abs(y);

    vec4 b0 = vec4(x.xy, y.xy);
    vec4 b1 = vec4(x.zw, y.zw);

  //vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
  //vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
    vec4 s0 = floor(b0) * 2.0 + 1.0;
    vec4 s1 = floor(b1) * 2.0 + 1.0;
    vec4 sh = -step(h, vec4(0.0));

    vec4 a0 = b0.xzyw + s0.xzyw * sh.xxyy;
    vec4 a1 = b1.xzyw + s1.xzyw * sh.zzww;

    vec3 p0 = vec3(a0.xy, h.x);
    vec3 p1 = vec3(a0.zw, h.y);
    vec3 p2 = vec3(a1.xy, h.z);
    vec3 p3 = vec3(a1.zw, h.w);

//Normalise gradients
    vec4 norm = taylorInvSqrt(vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3)));
    p0 *= norm.x;
    p1 *= norm.y;
    p2 *= norm.z;
    p3 *= norm.w;

// Mix final noise value
    vec4 m = max(0.6 - vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3)), 0.0);
    m = m * m;
    return 42.0 * dot(m * m, vec4(dot(p0, x0), dot(p1, x1), dot(p2, x2), dot(p3, x3)));
}

float accumulateNoises(Noise[10] noises, vec3 pos, float height)
{
    float result = 0.0;

    for(int i = 0; i < noises.length(); i++)
    {
        Noise noise = noises[i];

        if(noise.factor == 0.0)
            continue;

        float acc = 0.0;
        for(int j = 0; j < noise.octaves; j++)
        {
            acc += noise.gain * snoise(vec3(pos * noise.lacunarity));
        }

        result += noise.factor * acc;
    }

    return result;
}

// Used to create mountains and valleys
Surface applyUpperNoises(Surface s, vec3 pos, Noise[10] noises)
{
    Surface result = s;
    result.sd += accumulateNoises(noises, pos, s.sd);
    return result;
}

// Used to create oceans and lakes reliefs
Surface applyLowerNoises(Surface s, vec3 pos, Noise[10] noises)
{
    Surface result = s;
    result.sd -= abs(accumulateNoises(noises, pos, s.sd));
    return result;
}

//---------------------------------

Surface map(in vec3 pos)
{
    vec3 k = pos;
    Surface s = sdSphere(k, 2.5, vec3(0.0));

    Noise[10] planetUpperNoises;
    planetUpperNoises[0] = Noise(0.4, 1, 0.65, 0.7);
    planetUpperNoises[1] = Noise(0.3, 1, 1.4, 0.2);
    planetUpperNoises[2] = Noise(0.3, 3, 3., 0.01);
    planetUpperNoises[3] = Noise(0.2, 5, 9., 0.005);

    Surface s2 = applyUpperNoises(s, pos, planetUpperNoises);
    s2.sd = min(s.sd + 0.03, s2.sd);
    s2.col = apply_coloring_by_height(s);

    return s2;
}

vec3 calcNormal(in vec3 pos)
{
    const float ep = 0.0001;
    vec2 e = vec2(1.0, -1.0) * 0.5773;
    return normalize(e.xyy * (map(pos + e.xyy * ep).sd) +
        e.yyx * (map(pos + e.yyx * ep).sd) +
        e.yxy * (map(pos + e.yxy * ep).sd) +
        e.xxx * (map(pos + e.xxx * ep).sd));
}

float calcSoftshadow(in vec3 ro, in vec3 rd, float tmin, float tmax, const float k)
{
    float res = 1.0;
    float t = tmin;
    for(int i = 0; i < 50; i++)
    {
        Surface s = map(ro + rd * t);
        float h = s.sd;
        res = min(res, k * h / t);
        t += clamp(h, 0.02, 0.20);
        if(res < 0.005 || t > tmax)
            break;
    }
    return clamp(res, 0.0, 1.0);
}

// Camera

// Set up a camera looking at the scene.
// origin - camera is positioned relative to, and looking at, this point
// dist(ance) - how far camera is from origin
// rotation - about x & y axes, by left-hand screw rule, relative to camera looking along +z
// zoom- the relative length of the lens
vec3 localRay;
void handleCamera(out vec3 pos, out vec3 ray, in vec3 origin, in vec2 rotation, in float dist, in float zoom, in vec2 fragCoord)
{
	// get rotation coefficients
    vec2 c = vec2(cos(rotation.x), cos(rotation.y));
    vec4 s;
    s.xy = vec2(sin(rotation.x), sin(rotation.y)); // worth testing if this is faster as sin or sqrt(1.0-cos);
    s.zw = -s.xy;

	// ray in view space
    ray.xy = fragCoord.xy - iResolution.xy * .5;
    ray.z = iResolution.y * zoom;
    ray = normalize(ray);
    localRay = ray;

	// rotate ray
    ray.yz = ray.yz * c.xx + ray.zy * s.zx;
    ray.xz = ray.xz * c.yy + ray.zx * s.yw;

	// position camera
    pos = origin - dist * vec3(c.x * s.y, s.z, c.x * c.y);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec3 tot = vec3(0.0);

   // 1) Set up the Camera (primary ray for this pixel) 
    vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;

   // Camera Handling
    vec2 cameraRotation = vec2(.5, .5) + vec2(-.35, 4.5) * (iMouse.yx / iResolution.yx);
    vec3 ro, rd;
    float minRenderDistance = 0.1, maxRenderDistance = 20.0;

    handleCamera(ro, rd, vec3(0.), cameraRotation, 10.0, 1., fragCoord);

    Surface s;
   // 2) Raymarching 
    for(int i = 0; i < 256; i++)
    {
        vec3 p = ro + minRenderDistance * rd;
        s = map(p);
        float h = s.sd;
        if(abs(h) < 0.0001 || minRenderDistance > maxRenderDistance)
            break;
        minRenderDistance += h;
    }

    vec3 col = vec3(0.0);

   // 3) Shading
    vec3 ambiant_color = vec3(0.0, 0.02, 0.03);
    vec3 light_color = vec3(0.97, 0.85, 0.78);

    if(minRenderDistance < maxRenderDistance)
    {
        vec3 pos = ro + minRenderDistance * rd;
        vec3 nor = calcNormal(pos);
        vec3 lig = normalize(vec3(0.8, 0.8,-1.));
        lig = rotate(lig, vec3(0.1, 0.1, 0.), 1.3 * iTime); // Rotate the light around the planet
        float dif = clamp(dot(nor, lig), 0.0, 1.0);
        float sha = calcSoftshadow(pos, lig, 0.001, 1.0, 16.0);
        float amb = 0.5 + 0.5 * nor.y;
        col = ambiant_color * amb + light_color * dif * sha * s.col;
    }
    else
    {
        vec3 k = vec3(p, 0.);
        vec3 atmColor = vec3(0.12, 0.53, 0.88);
        Surface atm = sdSphere(k, 0.6, atmColor);

        // Coloring depending on the distance to create a atmosphere effect
        if(atm.sd < 0.0)
        {
            float d = atm.sd * 1.2;
            col = atm.col * smoothstep(0.0, 0.1, -d);
        }
        else
        {
            // use the noise to add some colored stars in the sky
            Noise[10] starNoise1;
            starNoise1[0] = Noise(1., 3, 30.0, 0.43);

            float n = accumulateNoises(starNoise1, k, 0.0);

            if(n > 0.98)
            {
                vec3 starColor = get_random_star_color(n, k);

                // Add random twinkle effect
                float twinkle = 0.5 + 0.5 * sin(iTime * 3.0 + n * 100.0);

                vec3 star = starColor * n * twinkle;

                col = star;
            }
        }
    }

    // Gamma correction
    col = sqrt(col);

    tot += col;

    fragColor = vec4(tot, 1.0);
}