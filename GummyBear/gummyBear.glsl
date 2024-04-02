#iChannel0 "file://background.glsl"

struct Surface
{
  float sd; // signed distance value
  vec3 col; // color
  vec3 emission; // emission color
  float roughness; // surface roughness
  float metallic; // metallic property
};

///////////////////////////////////////////////////////////

float rand(vec2 co)
{
  return fract(sin(dot(co.xy, vec2(12.9898, 78.233)) * 43758.5453));
}

float opUnion(float d1, float d2)
{
  return min(d1, d2);
}

Surface opUnion(Surface s1, Surface s2)
{
  if(s1.sd < s2.sd)
    return s1;
  else
    return s2;
}

float opSmoothUnion(float d1, float d2, float k)
{
  float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
  return mix(d2, d1, h) - k * h * (1.0 - h);
}

float sdRoundBox(vec3 p, vec3 b, float r)
{
  vec3 q = abs(p) - b + r;
  return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0) - r;
}

float sdSphere(vec3 p, float s)
{
  return length(p) - s;
}

float sdCutSphere(vec3 p, float r, float h)
{
  // sampling independent computations (only depend on shape)
  float w = sqrt(r * r - h * h);

  // sampling dependant computations
  vec2 q = vec2(length(p.xz), p.y);
  float s = max((h - r) * q.x * q.x + w * w * (h + r - 2.0 * q.y), h * q.x - w * q.y);
  return (s < 0.0) ? length(q) - r : (q.x < w) ? h - q.y : length(q - vec2(w, h));
}

vec3 rotate(vec3 p, vec3 angle)
{
  float a = angle.x;
  float b = angle.y;
  float c = angle.z;
  float sina = sin(a);
  float cosa = cos(a);
  float sinb = sin(b);
  float cosb = cos(b);
  float sinc = sin(c);
  float cosc = cos(c);
  mat3 ma = mat3(cosa, 0, sina, 0, 1, 0, -sina, 0, cosa);
  mat3 mb = mat3(1, 0, 0, 0, cosb, -sinb, 0, sinb, cosb);
  mat3 mc = mat3(cosc, -sinc, 0, sinc, cosc, 0, 0, 0, 1);
  return mc * mb * ma * p;
}

float sdRotatedCutSphere(vec3 p, float r, float h, vec3 angle)
{
  p = rotate(p, angle); // Rotate the point
  return sdCutSphere(p, r, h);
}

float sdPlane(vec3 p, vec3 n, float h)
{
  // n must be normalized
  return dot(p, n) + h;
}

Surface sdPlane(vec3 p, vec3 n, float h, vec3 col, vec3 emission, float roughness, float metallic)
{
  Surface s;
  s.sd = sdPlane(p, n, h);
  s.col = col;
  s.emission = emission;
  s.roughness = roughness;
  s.metallic = metallic;
  return s;
}

Surface sdGummyBear(vec3 pos)
{
  vec3 p = pos;
  float body = sdRoundBox(pos, vec3(.7, 2., 1.), 0.2);
  float bear = body;

    // Hands
  p = pos + vec3(-.7, -0.6, -.8);
  float rHand = sdSphere(p, .5);
  bear = opSmoothUnion(bear, rHand, 0.5);

  p = pos + vec3(-.7, -0.6, .8);
  float lHand = sdSphere(p, .5);
  bear = opSmoothUnion(bear, lHand, 0.5);

    // Foots
  p = pos + vec3(-.7, 2.1, -.6);
  float rFoot = sdCutSphere(p, .5, 0.1);
  bear = opSmoothUnion(bear, rFoot, 0.5);

  p = pos + vec3(-.7, 2.1, .6);
  float lFoot = sdCutSphere(p, .5, 0.1);
  bear = opSmoothUnion(bear, lFoot, 0.5);

  // Head
  p = pos + vec3(0., -2.4, 0.);
  float head = sdRoundBox(p, vec3(0.77, 0.7, 1.), 0.3);
  bear = opSmoothUnion(bear, head, 0.2);

  // Ears
  p = pos + vec3(-0.5, -3., -.8);
  float rEar = sdRotatedCutSphere(p, .5, 0.01, vec3(0., 0.5, 01.5));
  bear = opSmoothUnion(bear, rEar, 0.4);

  p = pos + vec3(-0.5, -3., .8);
  float lEar = sdRotatedCutSphere(p, .5, 0.01, vec3(0., -0.5, 01.5));
  bear = opSmoothUnion(bear, lEar, 0.4);

  vec3 color = vec3(1.0, 0.0, 0.0);
  vec3 emissionColor = color + vec3(0.5, 0.5, 0.5);

  return Surface(bear, color, emissionColor, 0.5, 0.2);
}

///////////////////////////////////////////////////////////////

Surface map(in vec3 pos)
{
  vec3 p = pos + vec3(0.0, -2.0, 0.0);
  Surface bear = sdGummyBear(p);

  Surface plane = sdPlane(pos, vec3(0.0, 1.0, 0.0), 0.0, vec3(1.0), vec3(0.0), 0., 0.1);

  bear = opUnion(bear, plane);

  return bear;
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

vec3 localRay;

void CamPolar(out vec3 pos, out vec3 ray, in vec3 origin, in vec2 rotation, in float distance, in float zoom, in vec2 fragCoord)
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
  pos = origin - distance * vec3(c.x * s.y, s.z, c.x * c.y);
}

//////////////////////////////////////////////////////////////

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
  vec3 tot = vec3(0.0);

  vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y; 
   //center of screen: (0,0)
   //dimensions: +/- 0.5

   //Better Camera 
  vec2 camRot = vec2(.5, .5) + vec2(-.35, 4.5) * (iMouse.yx / iResolution.yx);
  vec3 ro, rd;
  CamPolar(ro, rd, vec3(0), camRot, 20.0, 0.8, fragCoord);

  float t = 0.1;
  float tend = 25.;
  Surface s;
  for(int i = 0; i < 128; i++)
  {
    vec3 p = ro + t * rd;
    s = map(p);
    float h = s.sd;
    if(abs(h) < 0.0001 || t > tend)
      break;
    t += h;
  }

  vec3 col = texture(iChannel0, p / iResolution.xy).rgb;

  vec3 ambiant_color = vec3(0.0, 0.02, 0.03);
  vec3 light_color = vec3(0.97, 0.8, 0.8);
  if(t < tend)
  {
    vec3 pos = ro + t * rd;
    vec3 nor = calcNormal(pos);
    vec3 lig = normalize(vec3(0.8, 0.8, -1.));
    float dif = clamp(dot(nor, lig), 0.0, 1.0);
    float sha = calcSoftshadow(pos, lig, 0.001, 1.0, 16.0);
    float amb = 0.5 + 0.5 * nor.y;
    col = ambiant_color * amb + light_color * dif * sha * s.col;

    if(s.metallic > 0.0 || s.roughness > 0.0)
    {
    // Glass effect
      vec3 refl = reflect(rd, nor);
      vec3 refr = refract(rd, nor, 1.0 / 1.5);
      float fresnel = 0.1 + 0.9 * pow(1.0 - dot(-rd, nor), 5.0);
      col = ambiant_color * amb + light_color * dif * sha * s.col * (1.0 - fresnel) + s.emission * fresnel * texture(iChannel0, p + refr.xy * s.roughness).rgb + 0.1 * texture(iChannel0, p + refl.xy * s.metallic).rgb;
    }
  }

  col = sqrt(col);
  tot += col;

  fragColor = vec4(tot, 1.0);
}