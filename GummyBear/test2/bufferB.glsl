#iChannel0 "file://bufferA.glsl"
#iChannel1 "file://bufferB.glsl"

#define numboids 5.		// Number of boids (must be integer value represented as float)
#define speed 10.0f			// Boid speed
#define a1 0.001f			// Collision factor
#define a2 0.01f			// Cohesion factor	
#define a3 0.05f			// Alignment factor

vec3 getBoidPosition(float id)
{
    return texture(iChannel1, vec2(id + .5, .5) / iResolution.xy).xyz;
}

vec3 getBoidVelocity(float id)
{
    return texture(iChannel0, vec2(id + .5, .5) / iResolution.xy).xyz;
}

vec3 nrand(float id)
{
    return normalize(vec3(2.0 * fract(sin(id) * 43758.5453) - 1.0, 2.0 * fract(sin(id + 1.0) * 43758.5453) - 1.0, 2.0 * fract(sin(id + 2.0) * 43758.5453) - 1.0));

}

// Then the main image program queries the boid data from this buffer.
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{

    fragColor = vec4(0.0, 0.0, 0.0, 1.0);			// default data

    if (fragCoord.y > 1. || fragCoord.x > numboids) discard;

    float id = floor(fragCoord.x);

    // initialize random boid positions and velocities
    if(iFrame < 5)
    {
        vec3 init_p = nrand(id);
        fragColor = vec4(init_p, 1.);
    }
    else
    {        

        // boid of interest
        vec3 p = getBoidPosition(id);

        // loop over neighboring boids and build update vectors accordingly
        for(float nid = 0.0; nid < numboids; nid++)
        {
            p += iTimeDelta * getBoidVelocity(id);			// update position

            // Clamp p between 0 and 1
            p.xy = clamp(p.xy, -2.1, 2.1);
            p.z = clamp(p.z, -5., 1.);
        }

        fragColor = vec4(p.xyz, 1.);

    }
}