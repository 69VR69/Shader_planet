// "Boidies" by dr2 - 2018
// License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License

/**
 * Calculates the distance field value for a box shape.
 * 
 * @param p The position in 3D space.
 * @param b The size of the box in each dimension.
 * @return The distance field value.
 */
float PrBoxDf(vec3 p, vec3 b);

/**
 * Calculates the distance field value for a 2D box shape.
 * 
 * @param p The position in 2D space.
 * @param b The size of the box in each dimension.
 * @return The distance field value.
 */
float PrBox2Df(vec2 p, vec2 b);

/**
 * Calculates the distance field value for a cylindrical shape.
 * 
 * @param p The position in 3D space.
 * @param r The radius of the cylinder.
 * @param h The height of the cylinder.
 * @return The distance field value.
 */
float PrCylDf(vec3 p, float r, float h);

/**
 * Calculates the smooth minimum of two values.
 * 
 * @param a The first value.
 * @param b The second value.
 * @param r The smoothing factor.
 * @return The smoothed minimum value.
 */
float SmoothMin(float a, float b, float r);

/**
 * Calculates the smooth bump function.
 * 
 * @param lo The lower bound of the function.
 * @param hi The upper bound of the function.
 * @param w The width of the bump.
 * @param x The input value.
 * @return The smoothed bump value.
 */
float SmoothBump(float lo, float hi, float w, float x);

/**
 * Rotates a 2D vector by a given angle.
 * 
 * @param q The input vector.
 * @param a The rotation angle.
 * @return The rotated vector.
 */
vec2 Rot2D(vec2 q, float a);

/**
 * Converts a quaternion to a rotation matrix.
 * 
 * @param q The input quaternion.
 * @return The rotation matrix.
 */
mat3 QtToRMat(vec4 q);

/**
 * Generates a noise value for a 2D position.
 * 
 * @param p The input position.
 * @return The noise value.
 */
float Noisefv2(vec2 p);

/**
 * Generates a fractal brownian motion value for a 2D position.
 * 
 * @param p The input position.
 * @return The fractal brownian motion value.
 */
float Fbm2(vec2 p);

/**
 * Varies the normal vector based on a position, normal vector, and factor.
 * 
 * @param p The position in 3D space.
 * @param n The normal vector.
 * @param f The factor.
 * @return The varied normal vector.
 */
vec3 VaryNf(vec3 p, vec3 n, float f);

/**
 * Loads a vector from a specific index in a variable array.
 * 
 * @param idVar The index of the variable array.
 * @return The loaded vector.
 */
vec4 Loadv4(int idVar);

/**
 * Calculates the distance field value for a boid shape.
 * 
 * @param p The position in 3D space.
 * @return The distance field value.
 */
float BoidDf(vec3 p);

/**
 * Calculates the distance to the nearest boid along a ray.
 * 
 * @param ro The ray origin.
 * @param rd The ray direction.
 * @return The distance to the nearest boid.
 */
float BoidRay(vec3 ro, vec3 rd);

/**
 * Calculates the normal vector for a boid shape.
 * 
 * @param p The position in 3D space.
 * @return The normal vector.
 */
vec3 BoidNf(vec3 p);

/**
 * Calculates the color for a boid shape.
 * 
 * @param ro The ray origin.
 * @return The color.
 */
vec4 BoidCol(vec3 ro);

/**
 * Calculates the distance field value for an object shape.
 * 
 * @param p The position in 3D space.
 * @return The distance field value.
 */
float ObjDf(vec3 p);

/**
 * Calculates the distance to the nearest object along a ray.
 * 
 * @param ro The ray origin.
 * @param rd The ray direction.
 * @return The distance to the nearest object.
 */
float ObjRay(vec3 ro, vec3 rd);

/**
 * Calculates the normal vector for an object shape.
 * 
 * @param p The position in 3D space.
 * @return The normal vector.
 */
vec3 ObjNf(vec3 p);

/**
 * Calculates the background color.
 * 
 * @param ro The ray origin.
 * @param rd The ray direction.
 * @return The background color.
 */
vec3 BgCol(vec3 ro, vec3 rd);

/**
 * Shows the scene based on the ray origin and direction.
 * 
 * @param ro The ray origin.
 * @param rd The ray direction.
 * @return The color of the scene.
 */
vec3 ShowScene(vec3 ro, vec3 rd);

/**
 * Shows the widgets based on the UV coordinates, canvas size, color, zoom variable, flight variable, and view mode.
 * 
 * @param uv The UV coordinates.
 * @param canvas The size of the canvas.
 * @param col The color.
 * @param zmVar The zoom variable.
 * @param flVar The flight variable.
 * @param vuMode The view mode.
 * @return The color of the widgets.
 */
vec3 ShowWg(vec2 uv, vec2 canvas, vec3 col, float zmVar, float flVar, float vuMode);

/**
 * The main image function.
 * 
 * @param fragColor The output color.
 * @param fragCoord The fragment coordinates.
 */
void mainImage(out vec4 fragColor, in vec2 fragCoord);

#iChannel0 "file://b.glsl"

/*
 Real boids (everybody knows what they are) in 3D. The red boid is the leader; birds of
 the other two colors try to follow the leader and also prefer to group with others of the same
 color. The leader flies a fixed path which sometimes takes it through the square hoop. 
 The region is enclosed by a hemispherical dome, and boids bounce off both the dome and the 
 ground (dome size and projected path are marked on ground).
 
 The two sliders control zoom and the various interaction parameters (all lumped together for
 simplicity, red is strongest). Click on the circle to select the view (fixed, tracking, or
 attached to the red boid); mouse can also be used to look around.

 Explore...
*/

float PrBoxDf(vec3 p, vec3 b);
float PrBox2Df(vec2 p, vec2 b);
float PrCylDf(vec3 p, float r, float h);
float SmoothMin(float a, float b, float r);
float SmoothBump(float lo, float hi, float w, float x);
vec2 Rot2D(vec2 q, float a);
mat3 QtToRMat(vec4 q);
float Noisefv2(vec2 p);
float Fbm2(vec2 p);
vec3 VaryNf(vec3 p, vec3 n, float f);
vec4 Loadv4(int idVar);

const int nBoid = 64;
vec3 qHit, sunDir;
float dstFar, tCur, regSz, hoopSz, hoopThk, hoopHt, vuMode;
int idBoid;
const float pi = 3.14159;

float BoidDf(vec3 p)
{
    vec3 q;
    q = p;
    q.z += 0.3;
    return 0.8 * SmoothMin(PrCylDf(p, 0.16 - 0.97 * p.z * p.z, 0.4), max(PrCylDf(q.xzy, 0.5, 0.02 - 0.072 * min(dot(q.xz, q.xz) +
        (0.5 - q.z) * (0.5 - q.z), 0.25)), -q.z), 0.05);
}

float BoidRay(vec3 ro, vec3 rd)
{
    mat3 bMat;
    vec3 q, qm, rdm;
    float dHit, dMin, d, tol;
    tol = 0.0005;
    dMin = dstFar;
    for(int n = 0; n < nBoid; n++)
    {
        bMat = QtToRMat(Loadv4(3 * n + 2));
        qm = bMat * (ro - Loadv4(3 * n).xyz);
        rdm = bMat * rd;
        dHit = 0.;
        for(int j = 0; j < 120; j++)
        {
            q = qm + dHit * rdm;
            d = BoidDf(q);
            dHit += d;
            if(d < tol || dHit > dstFar)
                break;
        }
        if(d < tol && dHit < dMin)
        {
            dMin = dHit;
            idBoid = n;
            qHit = q;
        }
    }
    return dMin;
}

vec3 BoidNf(vec3 p)
{
    mat3 bMat;
    vec4 v;
    vec3 vn;
    vec2 e = vec2(0.0002, -0.0002);
    bMat = QtToRMat(Loadv4(3 * idBoid + 2));
    p = bMat * (p - Loadv4(3 * idBoid).xyz);
    v = vec4(BoidDf(p + e.xxx), BoidDf(p + e.xyy), BoidDf(p + e.yxy), BoidDf(p + e.yyx));
    vn = normalize(vec3(v.x - v.y - v.z - v.w) + 2. * v.yzw);
    if(vuMode == 2. || vuMode == 3.)
        vn.xy = Rot2D(vn.xy, -0.3 * sin(32. * 2. * pi * qHit.x / 0.5) * smoothstep(0., 0.2, -qHit.z) *
            step(-0.3, qHit.z) * (1. - smoothstep(0.01, 0.02, abs(qHit.y))));
    return vn * bMat;
}

vec4 BoidCol(vec3 ro)
{
    vec4 objCol;
    objCol.a = 0.1;
    if(qHit.z > 0.35)
        objCol.rgb = mix(vec3(0.9), vec3(0.2), (1. - smoothstep(0.004, 0.005, abs(qHit.y))) *
            SmoothBump(0.2, 0.8, 0.05, mod(128. * abs(qHit.x), 1.)));
    else if(length(qHit.yz - vec2(0.08, 0.18)) < 0.01)
        objCol.rgb = vec3(1., 1., 0.3);
    else if(length(qHit.yz - vec2(0.08, 0.18)) < 0.03)
        objCol.rgb = vec3(0.1, 0.1, 0.4);
    else if(idBoid == 0)
        objCol.rgb = vec3(1., 0., 0.) * (0.8 + 0.2 * sign(qHit.y));
    else
        objCol.rgb = ((Loadv4(3 * idBoid).w == 0.) ? vec3(1., 1., 0.) : vec3(0., 1., 1.)) * (0.8 + 0.2 * sign(qHit.y));
    return objCol;
}

float ObjDf(vec3 p)
{
    vec3 q;
    float d;
    q = p;
    q.y -= hoopHt;
    d = max(PrBoxDf(q, vec3(hoopSz, hoopSz, hoopThk)), -PrBox2Df(q.xy, vec2(hoopSz - 2. * hoopThk)));
    return d;
}

float ObjRay(vec3 ro, vec3 rd)
{
    vec3 p;
    float dHit, d;
    dHit = 0.;
    for(int j = 0; j < 120; j++)
    {
        p = ro + dHit * rd;
        d = ObjDf(ro + dHit * rd);
        dHit += d;
        if(d < 0.0005 || dHit > dstFar)
            break;
    }
    return dHit;
}

vec3 ObjNf(vec3 p)
{
    vec4 v;
    vec2 e = vec2(0.0002, -0.0002);
    v = vec4(ObjDf(p + e.xxx), ObjDf(p + e.xyy), ObjDf(p + e.yxy), ObjDf(p + e.yyx));
    return normalize(vec3(v.x - v.y - v.z - v.w) + 2. * v.yzw);
}

vec3 BgCol(vec3 ro, vec3 rd)
{
    vec3 vn, col;
    if(rd.y >= 0.)
    {
        col = vec3(0.1, 0.2, 0.4) + 0.2 * pow(1. - rd.y, 8.) +
            0.35 * pow(max(dot(rd, sunDir), 0.), 6.);
        col = mix(col, vec3(1.), clamp(0.1 + 0.8 * rd.y *
            Fbm2(5. * rd.xz / max(rd.y, 0.001)), 0., 1.));
    }
    else
    {
        ro -= rd * ro.y / rd.y;
        col = mix(mix(vec3(0.3, 0.4, 0.1), vec3(0.4, 0.5, 0.2), Fbm2(ro.xz)) *
            (1. - 0.2 * Noisefv2(32. * ro.xz)), vec3(0.35, 0.45, 0.65), pow(1. + rd.y, 32.));
        col *= (0.1 + 0.9 * max(dot(VaryNf(ro, vec3(0., 1., 0.), 2. * (1. - smoothstep(0.1, 1., length(ro.xz) / dstFar))), sunDir), 0.));
        col = mix(col, vec3(1., 0.5, 0.5), 0.2 * (1. - smoothstep(0.2, 0.3, abs(length(ro.xz) - regSz))));
        col = mix(col, vec3(1., 1., 0.5), 0.2 * (1. - smoothstep(0.2, 0.3, abs(length(vec2(abs(ro.x) - 0.35 * regSz, ro.z)) - 0.35 * regSz))));
        col = mix(col, vec3(0.7), 0.2 * (1. - smoothstep(0., 0.1, max(abs(ro.x) - hoopSz, abs(ro.z) - hoopThk))));
    }
    return col;
}

vec3 ShowScene(vec3 ro, vec3 rd)
{
    vec4 objCol;
    vec3 col, vn;
    float dstBoid, dstObj;
    dstBoid = BoidRay(ro, rd);
    dstObj = ObjRay(ro, rd);
    if(min(dstObj, dstBoid) < dstFar)
    {
        if(dstBoid < min(dstObj, dstFar))
        {
            ro += rd * dstBoid;
            objCol = BoidCol(ro);
            vn = BoidNf(ro);
        }
        else
        {
            ro += rd * dstObj;
            objCol = vec4(0.9, 0.6, 0.2, 0.2);
            vn = VaryNf(16. * ro, ObjNf(ro), 1.);
        }
        col = objCol.rgb * (0.3 + 0.7 * max(dot(vn, sunDir), 0.)) +
            objCol.a * pow(max(dot(normalize(sunDir - rd), vn), 0.), 128.);
    }
    else
        col = BgCol(ro, rd);
    return clamp(col, 0., 1.);
}

vec3 ShowWg(vec2 uv, vec2 canvas, vec3 col, float zmVar, float flVar, float vuMode)
{
    vec4 wgBx[3];
    vec2 ust;
    float asp;
    asp = canvas.x / canvas.y;
    wgBx[0] = vec4(0.45 * asp, -0.2, 0.023, 0.);
    wgBx[1] = vec4(0.43 * asp, 0.05, 0.01 * asp, 0.15);
    wgBx[2] = vec4(0.48 * asp, 0.05, 0.01 * asp, 0.15);
    if(abs(length(0.5 * uv - wgBx[0].xy) - wgBx[0].z) * canvas.y < 1.5)
        col = ((vuMode == 0.) ? vec3(0., 1., 0.) : (vuMode == 1.) ? vec3(1., 1., 0.) : vec3(1., 0., 0.));
    ust = abs(0.5 * uv - wgBx[1].xy) - wgBx[1].zw;
    if(abs(max(ust.x, ust.y)) * canvas.y < 1.)
        col = vec3(0.8);
    ust = 0.5 * uv - wgBx[1].xy;
    ust.y -= (zmVar - 0.5) * 2. * wgBx[1].w;
    if(length(ust) < 0.9 * wgBx[1].z)
        col = vec3(0.7, 0.1, 0.7);
    ust = abs(0.5 * uv - wgBx[2].xy) - wgBx[2].zw;
    if(abs(max(ust.x, ust.y)) * canvas.y < 1.)
        col = vec3(0.8);
    ust = 0.5 * uv - wgBx[2].xy;
    ust.y -= (flVar - 0.5) * 2. * wgBx[2].w;
    if(abs(length(ust) - 0.7 * wgBx[2].z) * canvas.y < 3.)
        col = (flVar > 0.7) ? vec3(1., 0.1, 0.1) : ((flVar > 0.3) ? vec3(1., 1., 0.1) : vec3(0.1, 1., 0.1));
    return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    mat3 vuMat, flMat;
    vec4 stDat, mPtr;
    vec3 rd, ro, vd, rLd, vLd, col;
    vec2 mMid[4], ut[4], mSize, canvas, uv, ori, ca, sa;
    float tCur, el, az, zmFac, pDist, zmVar, flVar, vuCorn;
    int wgSel;
    canvas = iResolution.xy;
    uv = 2. * fragCoord.xy / canvas - 1.;
    uv.x *= canvas.x / canvas.y;
    tCur = iTime;
    mSize = 0.22 * vec2(canvas.x / canvas.y, 1.);
    mMid[0] = (1. / mSize.y - 1.) * mSize;
    mMid[1] = mMid[0] * vec2(1., -1.);
    mMid[2] = mMid[0] * vec2(-1., -1.);
    mMid[3] = mMid[0] * vec2(-1., 1.);
    for(int k = 0; k < 4; k++) ut[k] = abs(uv - mMid[k]) - mSize;
    vuCorn = 0.;
    for(int k = 0; k < 4; k++)
    {
        if(max(ut[k].x, ut[k].y) < 0.)
        {
            uv = (uv - mMid[k]) / mSize.y;
            vuCorn = float(k + 1);
            break;
        }
    }
    stDat = Loadv4(3 * nBoid + 0);
    vuMode = mod(stDat.y + vuCorn, 5.);
    zmVar = (vuCorn == 0.) ? stDat.z : 0.5;
    regSz = stDat.w;
    stDat = Loadv4(3 * nBoid + 1);
    hoopSz = stDat.x;
    hoopThk = stDat.y;
    hoopHt = stDat.z;
    stDat = Loadv4(3 * nBoid + 2);
    el = (vuCorn == 0.) ? stDat.x : 0.;
    az = (vuCorn == 0.) ? stDat.y : 0.;
    flVar = stDat.z;
    wgSel = int(stDat.w);
    mPtr = Loadv4(3 * nBoid + 3);
    rLd = Loadv4(0).xyz;
    dstFar = 4. * regSz;
    if(vuMode == 0.)
    {
        zmFac = 2. + 8. * zmVar;
        el = clamp(el, -0.2 * pi, 0.3 * pi);
    }
    else if(vuMode == 1.)
    {
        ro = vec3(0., 2., -regSz);
        vd = rLd - ro;
        pDist = length(vd);
        vd = normalize(vd);
        zmFac = 0.7 + 4. * pDist / 50. + 25. * zmVar;
        az = clamp(0.25 * az, -0.2 * pi, 0.2 * pi);
        el = clamp(0.25 * el, -0.2 * pi, 0.2 * pi);
        az += 0.5 * pi + atan(-vd.z, vd.x);
        el += asin(vd.y);
    }
    else if(vuMode >= 2.)
    {
        az = clamp(az, -pi, pi);
        el = clamp(el, -0.3 * pi, 0.3 * pi);
        zmFac = 0.5 + 4.5 * zmVar;
    }
    ori = vec2(el, az);
    ca = cos(ori);
    sa = sin(ori);
    vuMat = mat3(ca.y, 0., -sa.y, 0., 1., 0., sa.y, 0., ca.y) *
        mat3(1., 0., 0., 0., ca.x, -sa.x, 0., sa.x, ca.x);
    rd = vuMat * normalize(vec3(uv, zmFac));
    if(vuMode == 0.)
    {
        ro = vuMat * vec3(0., 2., -2. * regSz);
        ro.y = max(ro.y, 0.4);
    }
    else if(vuMode >= 2.)
    {
        flMat = QtToRMat(Loadv4(2));
        ro = rLd + ((vuMode == 2.) ? vec3(0., 0.25, 0.) : ((vuMode == 3.) ? vec3(0., 0.3, -1.7) : vec3(0., 0.25, 1.5))) * flMat;
        if(vuMode == 4.)
            rd.z *= -1.;
        rd = rd * flMat;
    }
    sunDir = normalize(vec3(1., 3., -1.));
    col = ShowScene(ro, rd);
    if(vuCorn == 0.)
        col = ShowWg(uv, canvas, col, zmVar, flVar, vuMode);
    if(mPtr.z > 0. && wgSel < 0 && vuMode > 0. &&
        max(abs(uv.x), abs(uv.y)) < 0.03 &&
        min(abs(uv.x), abs(uv.y)) < 0.003)
        col = vec3(0.8, 0.8, 0.1);
    for(int k = 0; k < 4; k++)
    {
        if(max(ut[k].x, ut[k].y) < 0. && min(abs(ut[k].x), abs(ut[k].y)) * canvas.y < 2.)
            col = vec3(0.8, 0.8, 0.2);
    }
    fragColor = vec4(col, 1.);
}

float PrBoxDf(vec3 p, vec3 b)
{
    vec3 d;
    d = abs(p) - b;
    return min(max(d.x, max(d.y, d.z)), 0.) + length(max(d, 0.));
}

float PrBox2Df(vec2 p, vec2 b)
{
    vec2 d;
    d = abs(p) - b;
    return min(max(d.x, d.y), 0.) + length(max(d, 0.));
}

float PrCylDf(vec3 p, float r, float h)
{
    return max(length(p.xy) - r, abs(p.z) - h);
}

float SmoothMin(float a, float b, float r)
{
    float h;
    h = clamp(0.5 + 0.5 * (b - a) / r, 0., 1.);
    return mix(b, a, h) - r * h * (1. - h);
}

float SmoothBump(float lo, float hi, float w, float x)
{
    return (1. - smoothstep(hi - w, hi + w, x)) * smoothstep(lo - w, lo + w, x);
}

vec2 Rot2D(vec2 q, float a)
{
    return q * cos(a) + q.yx * sin(a) * vec2(-1., 1.);
}

mat3 QtToRMat(vec4 q)
{
    mat3 m;
    float a1, a2, s;
    q = normalize(q);
    s = q.w * q.w - 0.5;
    m[0][0] = q.x * q.x + s;
    m[1][1] = q.y * q.y + s;
    m[2][2] = q.z * q.z + s;
    a1 = q.x * q.y;
    a2 = q.z * q.w;
    m[0][1] = a1 + a2;
    m[1][0] = a1 - a2;
    a1 = q.x * q.z;
    a2 = q.y * q.w;
    m[2][0] = a1 + a2;
    m[0][2] = a1 - a2;
    a1 = q.y * q.z;
    a2 = q.x * q.w;
    m[1][2] = a1 + a2;
    m[2][1] = a1 - a2;
    return 2. * m;
}

const float cHashM = 43758.54;

vec2 Hashv2v2(vec2 p)
{
    vec2 cHashVA2 = vec2(37., 39.);
    return fract(sin(vec2(dot(p, cHashVA2), dot(p + vec2(1., 0.), cHashVA2))) * cHashM);
}

float Noisefv2(vec2 p)
{
    vec2 t, ip, fp;
    ip = floor(p);
    fp = fract(p);
    fp = fp * fp * (3. - 2. * fp);
    t = mix(Hashv2v2(ip), Hashv2v2(ip + vec2(0., 1.)), fp.y);
    return mix(t.x, t.y, fp.x);
}

float Fbm2(vec2 p)
{
    float f, a;
    f = 0.;
    a = 1.;
    for(int i = 0; i < 5; i++)
    {
        f += a * Noisefv2(p);
        a *= 0.5;
        p *= 2.;
    }
    return f * (1. / 1.9375);
}

float Fbmn(vec3 p, vec3 n)
{
    vec3 s;
    float a;
    s = vec3(0.);
    a = 1.;
    for(int i = 0; i < 5; i++)
    {
        s += a * vec3(Noisefv2(p.yz), Noisefv2(p.zx), Noisefv2(p.xy));
        a *= 0.5;
        p *= 2.;
    }
    return dot(s, abs(n));
}

vec3 VaryNf(vec3 p, vec3 n, float f)
{
    vec3 g;
    vec2 e = vec2(0.1, 0.);
    g = vec3(Fbmn(p + e.xyy, n), Fbmn(p + e.yxy, n), Fbmn(p + e.yyx, n)) - Fbmn(p, n);
    return normalize(n + f * (g - n * dot(n, g)));
}

#define txBuf iChannel0
#define txSize iChannelResolution[0].xy

const float txRow = 128.;

vec4 Loadv4(int idVar)
{
    float fi;
    fi = float(idVar);
    return texture(txBuf, (vec2(mod(fi, txRow), floor(fi / txRow)) + 0.5) /
        txSize);
}
