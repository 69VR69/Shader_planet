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

//---------------------------------
//includes Coloring
Surface map(in vec3 pos) {
    vec3 k = pos;
    Surface s = sdSphere(k, 0.5, vec3(0.));

    // Add some noise
    Surface s1 = apply_octave_noise(s, k, 4, .5);
    Surface s2 = apply_octave_noise(s1, k * 2.0, 2, .5);
    s = apply_octave_noise(s, k * 4.0, 3, .5);

    // Coloring depinding on the distance
    const vec3 high_color = vec3(1.0, 0.0, 0.0);
    const float high_level = 0.8;
    const vec3 mid_color = vec3(0.0, 1.0, 0.0);
    const float mid_level = 0.5;
    const vec3 low_color = vec3(0.0, 0.0, 1.0);
    const float low_level = 0.0;

    // Coloring
    float d = s.sd*2.;
    
    if(d>=high_level){
        s2.col = high_color;
    }else if(d>=mid_level){
        s2.col = mix(mid_color, high_color, (d-mid_level)/(high_level-mid_level));
    }else if(d>=low_level){
        s2.col = mix(low_color, mid_color, (d-low_level)/(mid_level-low_level));
    }else{
        s2.col = low_color;
    }

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