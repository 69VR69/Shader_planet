struct Surface {
    float sd; // signed distance value
    vec3 col; // color
};

Surface sdSphere(vec3 p, float s, vec3 col) {
    float d = length(p) - s;
    return Surface(d, col);
}

//Noise
float hash(float n) {
    return fract(sin(n) * 43758.5453);
}

float noise(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    float n = p.x + p.y * 157.0 + 113.0 * p.z;
    return mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
                   mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
               mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
                   mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
}

Surface apply_octave_noise(Surface s, vec3 x, int octaves) {
    float n = 0.0;
    float a = 1.0;
    for(int i = 0; i < octaves; i++) {
        n += a * noise(x);
        x *= 2.0;
        a *= 0.5;
    }
    s.sd *= n;
    return s;
}

//---------------------------------
//includes Coloring
Surface map(in vec3 pos) {
    vec3 k = pos;
    Surface s = sdSphere(k, 0.5, vec3(0.0, 0.0, 1.0));

    // Add some noise
    s = apply_octave_noise(s, k, 4);
    s = apply_octave_noise(s, k * 2.0, 2);

    return s;
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
    for(int i = 0; i < 128; i++) {
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
        float dif = clamp(dot(nor, lig), 0.0, 1.0);
        float sha = calcSoftshadow(pos, lig, 0.001, 1.0, 16.0);
        float amb = 0.5 + 0.5 * nor.y;
        col = vec3(0.05, 0.1, 0.15) * amb + vec3(1.00, 0.9, 0.80) * dif * sha * s.col;
    }

    col = sqrt(col);
    tot += col;

    fragColor = vec4(tot, 1.0);
}