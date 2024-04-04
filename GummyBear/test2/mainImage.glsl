#iChannel0 "file://bufferA.glsl"
#iChannel1 "file://bufferB.glsl"

#define numboids 20.		// number of boids (must be integer value represented as float)

float opUnion(float d1, float d2)
{
    return min(d1, d2);
}

float sdBox(vec3 p, vec3 b)
{
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

vec3 getBoidPosition(float id)
{
    return texture(iChannel1, vec2(id + .5, .5) / iResolution.xy).xyz;
}

vec3 getBoidVelocity(float id)
{
    return texture(iChannel0, vec2(id + .5, .5) / iResolution.xy).xyz;
}

float map(in vec3 pos)
{

    float res = 1000.0;

    for(int i = 0; i < int(numboids); i++)
    {

        vec3 boidPosition = getBoidPosition(float(i));

        vec3 p = pos - boidPosition;
        float d = sdBox(p, vec3(0.05));

        res = opUnion(res, d);
    }
    return res;
}

vec3 calcNormal(in vec3 pos)
{
    const float ep = 0.0001;
    vec2 e = vec2(1.0, -1.0) * 0.5773;
    return normalize(e.xyy * map(pos + e.xyy * ep) +
        e.yyx * map(pos + e.yyx * ep) +
        e.yxy * map(pos + e.yxy * ep) +
        e.xxx * map(pos + e.xxx * ep));
}

float calcSoftshadow(in vec3 ro, in vec3 rd, float tmin, float tmax, const float k)
{
    float res = 1.0;
    float t = tmin;
    for(int i = 0; i < 50; i++)
    {
        float h = map(ro + rd * t);
        res = min(res, k * h / t);
        t += clamp(h, 0.02, 0.20);
        if(res < 0.005 || t > tmax)
            break;
    }
    return clamp(res, 0.0, 1.0);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec3 tot = vec3(0.0);

    vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;

    vec3 ro = vec3(0.0, 3.0, 8.4);
    vec3 rd = normalize(vec3(p - vec2(0.1, 1.9), -6.));

    float t = 0.01;
    for(int i = 0; i < 128; i++)
    {
        vec3 p = ro + t * rd;
        float h = map(p);
        if(abs(h) < 0.0001 || t > 11.0)
            break;
        t += h;
    }

    vec3 col = vec3(0.0);

    if(t < 11.0)
    {
        vec3 pos = ro + t * rd;
        vec3 nor = calcNormal(pos);
        vec3 lig = normalize(vec3(1.0, 0.8, -0.2));
        float dif = clamp(dot(nor, lig), 0.0, 1.0);
        float sha = calcSoftshadow(pos, lig, 0.001, 1.0, 16.0);
        float amb = 0.5 + 0.5 * nor.y;
        col = vec3(0.05, 0.1, 0.15) * amb + vec3(1.00, 0.9, 0.80) * dif * sha;
    }
    else
        col = vec3(0.0);

    col = sqrt(col);
    tot += col;

    fragColor = vec4(tot, 1.0);
}