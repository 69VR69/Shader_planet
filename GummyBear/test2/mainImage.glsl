#iChannel0 "file://bufferA.glsl"
#define numboids 20.		// number of boids (must be integer value represented as float)

float sdBox(vec3 p, vec3 b)
{
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

// A function that defines whether the screen coordinate uv is inside a boid.
// boid.xyz = boid screenspace coordinate
// boid.w = boid id
bool drawBoid(vec3 uv, vec4 boid)
{
    // Box
    vec3 position = boid.xyz;
    vec3 bounds = vec3(10.0, 10.0, 10.0);
    return sdBox(uv - position, bounds) < 0.0;
}

vec4 getBoid(float id)
{
    return texture(iChannel0, vec2(id + .5, .5) / iResolution.xy);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    //Background Color
    fragColor = vec4(0.1, 0.1, 0.1, 1.);

    for(float i = 0.; i < numboids; i++)
    {

        vec4 boid = getBoid(i);

        //Boid color
        vec4 col = vec4(boid.x / iResolution.x, 0., boid.y / iResolution.y, 1.);

        if(drawBoid(vec3(fragCoord, 0.), boid))
            fragColor = col;

    }

}