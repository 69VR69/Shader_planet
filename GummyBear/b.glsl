#iChannel0 "file://b.glsl"

vec4 RMatToQt(mat3 m);
float SmoothBump(float lo, float hi, float w, float x);
vec2 Rot2D(vec2 q, float a);
float Hashff(float p);
vec4 Loadv4(int idVar);
void Savev4(int idVar, vec4 val, inout vec4 fCol, vec2 fCoord);

const int nBoid = 64;
vec3 rLd, vLd, aLd;
float vFly, regSz, fSep, rFlok, fFlok, fAln, fLead, rAttr, nStep, dt, hoopSz, hoopThk, hoopHt;
const float pi = 3.14159;
const float txRow = 128.;

/**
 * Performs a step in the simulation for a boid
 * @param mId The ID of the boid
 * @param r The position of the boid (output)
 * @param v The velocity of the boid (output)
 * @param a The acceleration of the boid (output)
 * @param grp The group ID of the boid (output)
 */
void Step(int mId, out vec3 r, out vec3 v, out vec3 a, out float grp)
{
    vec4 p;
    vec3 dr, rSum, vSum;
    float nNeb, rLen, vMag, rMarg;
    p = Loadv4(3 * mId);
    r = p.xyz;
    grp = p.w;
    v = Loadv4(3 * mId + 1).xyz;
    a = vec3(0.);
    vSum = vec3(0.);
    rSum = vec3(0.);
    nNeb = 0.;
    for(int n = 0; n < nBoid; n++)
    {
        if(n != mId)
        {
            p = Loadv4(3 * n);
            dr = r - p.xyz;
            rLen = length(dr);
            if(rLen < 1.)
                a += fSep * (1. / rLen - 1.) * dr;
            if(rLen < rFlok && grp == p.w)
            {
                rSum += p.xyz;
                vSum += Loadv4(3 * n + 1).xyz;
                ++nNeb;
            }
        }
    }
    if(nNeb > 0.)
        a -= fFlok * (r - rSum / nNeb) + fAln * (v - vSum / nNeb);
    dr = r - rLd;
    rLen = length(dr);
    if(rLen < rAttr)
    {
        a += ((1. - 2. * smoothstep(2., 3., rLen)) * fLead / max(rLen * rLen, 0.001)) * dr;
    }
    rMarg = 1.;
    dr = r;
    dr.xy -= vec2((hoopSz - hoopThk) * sign(r.x), hoopHt);
    dr = max(abs(dr) - vec3(hoopThk, hoopSz, hoopThk), 0.) * sign(dr);
    rLen = length(dr);
    if(rLen < hoopThk + rMarg)
        a += fSep * ((hoopThk + rMarg) / rLen - 1.) * dr;
    dr = r;
    dr.y -= hoopHt + (hoopSz - hoopThk) * sign(r.y);
    dr = max(abs(dr) - vec3(hoopSz, hoopThk, hoopThk), 0.) * sign(dr);
    rLen = length(dr);
    if(rLen < hoopThk + rMarg)
        a += fSep * ((hoopThk + rMarg) / rLen - 1.) * dr;
    a += 0.05 * (vFly - length(v)) * normalize(v);
    v += dt * a;
    r += dt * v;
    rLen = length(r);
    if(rLen > regSz)
    {
        if(dot(r, v) > 0.)
            v = 0.9 * reflect(v, r / rLen);
        r *= (regSz - 0.05) / rLen;
    }
    if(r.y < 0.)
    {
        r.y = 0.05;
        if(v.y < 0.)
            v = 0.9 * reflect(v, vec3(0., 1., 0.));
    }
}

/**
 * Computes the position of a boid on the track
 * @param t The time parameter
 * @return The position of the boid on the track
 */
vec3 TrackPos(float t)
{
    vec3 r;
    float tt = mod(t, 1.);
    r.xz = 0.35 * regSz * ((mod(t, 2.) < 1.) ? vec2(-cos(2. * pi * tt) + 1., sin(2. * pi * tt)) : vec2(cos(2. * pi * tt) - 1., sin(2. * pi * tt)));
    r.y = hoopHt + (hoopSz + 7. * hoopThk) * ((mod(floor(t / 2. - 0.25), 4.) > 1.) ? SmoothBump(0.3, 0.7, 0.15, tt) : 1.);
    return r;
}

/**
 * Initializes the position, velocity, acceleration, and group ID of a boid
 * @param mId The ID of the boid
 * @param r The position of the boid (output)
 * @param v The velocity of the boid (output)
 * @param a The acceleration of the boid (output)
 * @param grp The group ID of the boid (output)
 */
void Init(int mId, out vec3 r, out vec3 v, out vec3 a, out float grp)
{
    float mIdf, w;
    mIdf = float(mId);
    w = mIdf / float(nBoid);
    r = TrackPos(2. * w);
    r.y += 1.5 * (w - 0.5);
    v = vFly * (2. * normalize(vec3(Hashff(mIdf), Hashff(mIdf + 0.3), Hashff(mIdf + 0.6)) - 0.5) + 1.);
    a = vec3(0.);
    grp = floor(2. * Hashff(mIdf));
}

/**
 * Sets the position, velocity, and acceleration of the lead boid
 * @param r The position of the lead boid (output)
 * @param v The velocity of the lead boid (output)
 * @param a The acceleration of the lead boid (output)
 */
void SetLead(out vec3 r, out vec3 v, out vec3 a)
{
    vec3 rf, rb;
    float s, ds;
    s = 5.5 * vFly * nStep * dt / (2. * pi * regSz);
    ds = 0.1;
    r = TrackPos(s);
    rf = TrackPos(s + ds);
    rb = TrackPos(s - ds);
    v = (rf - rb) / (2. * ds);
    a = (rf - 2. * r + rb) / (ds * ds);
}

/**
 * Computes the orientation of a boid
 * @param mId The ID of the boid
 * @param v The velocity of the boid
 * @param a The acceleration of the boid
 * @return The quaternion representing the orientation of the boid
 */
vec4 EvalOri(int mId, vec3 v, vec3 a)
{
    vec3 va, ort, ca, sa;
    float el, az, rl;
    v = normalize(v);
    va = cross(a, v);
    el = -0.7 * asin(v.y);
    az = atan(v.z, v.x) - 0.5 * pi;
    rl = 0.001 * length(va) * sign(va.y);
    ort = vec3(el, az, rl);
    ca = cos(ort);
    sa = sin(ort);
    return RMatToQt(mat3(ca.z, -sa.z, 0., sa.z, ca.z, 0., 0., 0., 1.) *
        mat3(1., 0., 0., 0., ca.x, -sa.x, 0., sa.x, ca.x) *
        mat3(ca.y, 0., -sa.y, 0., 1., 0., sa.y, 0., ca.y));
}

/**
 * The main image rendering function
 * @param fragColor The output fragment color
 * @param fragCoord The coordinates of the current fragment
 */
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec4 mPtr, mPtrP;
    vec4 wgBx[3], stDat, p;
    vec3 r, v, a;
    vec2 iFrag, canvas, ust;
    float tCur, grp, vuMode, asp, az, el, zmVar, flVar;
    int mId, pxId, wgSel, wgReg, kp;
    bool doInit;
    iFrag = floor(fragCoord);
    pxId = int(iFrag.x + txRow * iFrag.y);
    if(iFrag.x >= txRow || pxId >= 3 * nBoid + 4)
        discard;
    canvas = iResolution.xy;
    tCur = iTime;
    mPtr = iMouse;
    mPtr.xy = mPtr.xy / canvas - 0.5;
    mId = (pxId < 3 * nBoid) ? pxId / 3 : -1;
    vFly = 1.2;
    regSz = 40.;
    hoopSz = 2.5;
    hoopThk = 0.5;
    hoopHt = 5.;
    fSep = 10.;
    rFlok = 6.;
    dt = 0.05;
    wgReg = -2;
    doInit = false;
    if(iFrame <= 5)
    {
        mPtrP = mPtr;
        vuMode = 0.;
        zmVar = 0.5;
        flVar = 0.5;
        az = 0.;
        el = 0.;
        wgSel = -1;
        doInit = true;
    }
    else
    {
        stDat = Loadv4(3 * nBoid + 0);
        nStep = stDat.x;
        vuMode = stDat.y;
        zmVar = stDat.z;
        stDat = Loadv4(3 * nBoid + 2);
        el = stDat.x;
        az = stDat.y;
        flVar = stDat.z;
        wgSel = int(stDat.w);
        mPtrP = Loadv4(3 * nBoid + 3);
        fLead = mix(1., 10., flVar);
        rAttr = mix(3., 12., flVar);
        fAln = mix(0., 0.04, flVar);
        fFlok = mix(0., 0.04, flVar);
    }
    if(doInit)
    {
        nStep = 0.;
        if(mId >= 0)
            SetLead(rLd, vLd, aLd);
        if(mId > 0)
            Init(mId, r, v, a, grp);
    }
    else
    {
        ++nStep;
        if(mId >= 0)
            SetLead(rLd, vLd, aLd);
        if(mId > 0)
            Step(mId, r, v, a, grp);
    }
    if(mId == 0)
    {
        r = rLd;
        v = vLd;
        a = aLd;
        grp = 2.;
    }
    asp = canvas.x / canvas.y;
    if(mPtr.z > 0.)
    {
        wgBx[0] = vec4(0.45 * asp, -0.2, 0.023, 0.);
        wgBx[1] = vec4(0.43 * asp, 0.05, 0.01 * asp, 0.15);
        wgBx[2] = vec4(0.48 * asp, 0.05, 0.01 * asp, 0.15);
        if(length(mPtr.xy * vec2(asp, 1.) - wgBx[0].xy) < wgBx[0].z)
            wgReg = 0;
        ust = abs(mPtr.xy * vec2(asp, 1.) - wgBx[1].xy) - wgBx[1].zw;
        if(max(ust.x, ust.y) < 0.)
            wgReg = 1;
        ust = abs(mPtr.xy * vec2(asp, 1.) - wgBx[2].xy) - wgBx[2].zw;
        if(max(ust.x, ust.y) < 0.)
            wgReg = 2;
        if(mPtrP.z <= 0.)
            wgSel = wgReg;
    }
    else
    {
        wgSel = -1;
        wgReg = -2;
        az = (vuMode == 0.) ? 0.1 * pi : 0.;
        el = (vuMode == 0.) ? -0.04 * pi : 0.;
    }
    if(wgSel < 0)
    {
        if(mPtr.z > 0.)
        {
            az = 2. * pi * mPtr.x;
            el = pi * mPtr.y;
        }
    }
    else
    {
        if(wgSel == 0)
        {
            if(mPtrP.z <= 0.)
            {
                vuMode = mod(++vuMode, 5.);
                zmVar = 0.5;
                az = (vuMode == 0.) ? 0.1 * pi : 0.;
                el = (vuMode == 0.) ? -0.04 * pi : 0.;
            }
        }
        else if(wgSel == 1)
        {
            zmVar = clamp(0.5 + 0.5 * (mPtr.y - wgBx[1].y) / wgBx[1].w, 0., 1.);
        }
        else if(wgSel == 2)
        {
            flVar = clamp(0.5 + 0.5 * (mPtr.y - wgBx[2].y) / wgBx[2].w, 0., 1.);
        }
    }
    if(pxId < 3 * nBoid)
    {
        kp = 3 * mId;
        if(pxId == kp + 0)
            stDat = vec4(r, grp);
        else if(pxId == kp + 1)
            stDat = vec4(v, 0.);
        else if(pxId == kp + 2)
            stDat = EvalOri(mId, v, a);
    }
    else
    {
        kp = 3 * nBoid;
        if(pxId == kp + 0)
            stDat = vec4(nStep, vuMode, zmVar, regSz);
        else if(pxId == kp + 1)
            stDat = vec4(hoopSz, hoopThk, hoopHt, 0.);
        else if(pxId == kp + 2)
            stDat = vec4(el, az, flVar, float(wgSel));
        else if(pxId == kp + 3)
            stDat = mPtr;
    }
    Savev4(pxId, stDat, fragColor, fragCoord);
}

vec4 RMatToQt(mat3 m)
{
    vec4 q;
    const float tol = 1e-6;
    q.w = 0.5 * sqrt(max(1. + m[0][0] + m[1][1] + m[2][2], 0.));
    if(abs(q.w) > tol)
        q.xyz = vec3(m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0]) / (4. * q.w);
    else
    {
        q.x = sqrt(max(0.5 * (1. + m[0][0]), 0.));
        if(abs(q.x) > tol)
            q.yz = vec2(m[0][1], m[0][2]) / q.x;
        else
        {
            q.y = sqrt(max(0.5 * (1. + m[1][1]), 0.));
            if(abs(q.y) > tol)
                q.z = m[1][2] / q.y;
            else
                q.z = 1.;
        }
    }
    return normalize(q);
}

float SmoothBump(float lo, float hi, float w, float x)
{
    return (1. - smoothstep(hi - w, hi + w, x)) * smoothstep(lo - w, lo + w, x);
}

vec2 Rot2D(vec2 q, float a)
{
    return q * cos(a) + q.yx * sin(a) * vec2(-1., 1.);
}

const float cHashM = 43758.54;

float Hashff(float p)
{
    return fract(sin(p) * cHashM);
}

#define txBuf iChannel0
#define txSize iChannelResolution[0].xy

vec4 Loadv4(int idVar)
{
    float fi;
    fi = float(idVar);
    return texture(txBuf, (vec2(mod(fi, txRow), floor(fi / txRow)) + 0.5) /
        txSize);
}

void Savev4(int idVar, vec4 val, inout vec4 fCol, vec2 fCoord)
{
    vec2 d;
    float fi;
    fi = float(idVar);
    d = abs(fCoord - vec2(mod(fi, txRow), floor(fi / txRow)) - 0.5);
    if(max(d.x, d.y) < 0.5)
        fCol = val;
}
