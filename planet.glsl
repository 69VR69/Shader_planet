struct Surface {
    float sd; // signed distance value
    vec3 col; // color
};

struct Noise {
    float factor; // factor to multiply the noise by   
    int octaves; // number of octaves that represent how many times the noise is repeated
    float lacunarity; // factor that determines how much the frequency of the noise increases with each octave
    float gain; // factor that determines how much the amplitude of the noise decreases with each octave
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

float fbm(vec3 x, int octaves, float lacunarity, float gain) {
    float sum = 0.0;
    float amp = 1.0;
    float freq = 1.0;
    for(int i = 0; i < octaves; i++) {
        sum += amp * noise(x * freq);
        amp *= gain;
        freq *= lacunarity;
    }
    return sum;
}

Surface apply_noises(in Surface s, vec3 x, Noise[10] n) {
    float accumulated_noise = 0.;
    for(int i = 0; i < n.length(); i++) {
        if(n[i].factor == 0.) continue;
        accumulated_noise += n[i].factor * fbm(x, n[i].octaves, n[i].lacunarity, n[i].gain);
    }
    s.sd -= accumulated_noise;
    return s;
}

float getFactorSum(Noise[10] noises) {
    float sum = 0.;
    for(int i = 0; i < noises.length(); i++) {
        sum += noises[i].factor;
    }
    return sum;
}
 
vec3 apply_coloring_by_height(in Surface s, in float base) {
        // Coloring depending on the distance
        const vec3 high_color = vec3(0.89, 0.42, 0.04);
        float high_level = base + 0.7;
        const vec3 mid_color = vec3(0.12, 0.47, 0.04);
        float mid_level = base + 0.5;
        const vec3 low_color = vec3(0.13, 0.33, 0.88);
        float low_level = base + 0.0;

        // Coloring
        float d = s.sd*1.3;

        if (d >= high_level) {
            return high_color;
        } else if (d >= mid_level && d < high_level) {
            return mix(mid_color, high_color, (d - mid_level) / (high_level - mid_level));
        } else if (d >= low_level && d < mid_level) {
            return mix(low_color, mid_color, (d - low_level) / (mid_level - low_level));
        } else if (d < low_level){
            return low_color;
        }
    }

//---------------------------------

Surface map(in vec3 pos) {
    vec3 k = pos;
    Surface s = sdSphere(k, 0.5, vec3(0.));

    Noise[10] planetNoises;
    planetNoises[0] = Noise(0.52, 6, 2.2, 0.35);
    planetNoises[1] = Noise(0.6, 6, 1.4, 0.6);
    planetNoises[2] = Noise(0.3, 10, 5., 0.1);
    planetNoises[3] = Noise(0.2, 10, 9., 0.1);

    Surface s2 = apply_noises(s, k, planetNoises);

    float base = getFactorSum(planetNoises);

    s2.col = apply_coloring_by_height(s, base);

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

// Camera

// Set up a camera looking at the scene.
// origin - camera is positioned relative to, and looking at, this point
// dist(ance) - how far camera is from origin
// rotation - about x & y axes, by left-hand screw rule, relative to camera looking along +z
// zoom- the relative length of the lens

vec3 localRay;
void handleCamera(out vec3 pos, out vec3 ray, in vec3 origin, in vec2 rotation, in float dist, in float zoom, in vec2 fragCoord) {
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

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 tot = vec3(0.0);

   // 1) Set up the Camera (primary ray for this pixel) 
    vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;

   // Camera Handling
    vec2 cameraRotation = vec2(.5, .5) + vec2(-.35, 4.5) * (iMouse.yx / iResolution.yx);
    vec3 ro, rd;
    
    handleCamera(ro, rd, vec3(0.), cameraRotation, 10.0, 1., fragCoord);

    Surface s; 
   // 2) Raymarching 
    float t = 0.1;
    for(int i = 0; i < 256; i++) {
        vec3 p = ro + t * rd;
        s = map(p);
        float h = s.sd;
        if(abs(h) < 0.0001 || t > 20.0)
            break;
        t += h;
    }

    vec3 col = vec3(0.0);

   // 3) Shading
    vec3 ambiant_color = vec3(0.05, 0.1, 0.15);
    vec3 light_color = vec3(0.97, 0.85, 0.78);

    if(t < 20.0) {
        vec3 pos = ro + t * rd;
        vec3 nor = calcNormal(pos);
        vec3 lig = normalize(vec3(1.0, 1., 0.2));
        lig = rotate(lig, vec3(1.0, 1., 1.), 1.5 * iTime);
        float dif = clamp(dot(nor, lig), 0.0, 1.0);
        float sha = calcSoftshadow(pos, lig, 0.001, 1.0, 16.0);
        float amb = 0.5 + 0.5 * nor.y;
        col = vec3(0.05, 0.1, 0.15) * amb + light_color * dif * sha * s.col;
    }

    col = sqrt(col);
    tot += col;

    fragColor = vec4(tot, 1.0);
}