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

// Then the main image program queries the boid data from this buffer.
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{

    fragColor = vec4(0.0, 0.0, 0.0, 1.0);			// default data
    if (fragCoord.y > 1. || fragCoord.x > numboids) discard;

    float id = floor(fragCoord.x);

    // initialize random boid positions and velocities
    if(iFrame < 5)
    {
        vec3 init_v = vec3(sin(id), cos(id), 0.);
        init_v /= speed * length(init_v);
        fragColor = vec4(init_v, 1.);
    }
    else
    {        

        // boid of interest
        vec3 v = getBoidVelocity(id);
        vec3 p = getBoidPosition(id);

        // initialize velocity update vectors
        vec3 v_collision, v_cohesion, v_alignment = vec3(0.);

        // loop over neighboring boids and build update vectors accordingly
        for(float nid = 0.0; nid < numboids; nid++)
        {

            if(nid != id)
            {
                vec3 p_neighbor = getBoidPosition(id);	// nieghboring boid position
                vec3 v_neighbor = getBoidVelocity(id);	// neighboring boid velocity

                vec3 sep = p_neighbor - p;
                float ls = length(sep);
                float r = 50.;				// separation radius

            	// Collision Update
                if(ls < r)
                {
                    v_collision -= sep;
                }

                // Cohesion Update
                v_cohesion += sep;

                // Alignment Update
                v_alignment += v_neighbor;
            }

            // Perceived values
            v_cohesion = v_cohesion / (numboids - 1.0);
            v_alignment = v_alignment / (numboids - 1.0) - v;

            // Apply position and velocity updates
            v += a1 * v_collision + a2 * v_cohesion + a3 * v_alignment;
            v /= length(v) * speed;		// constrain velocity
        }

        fragColor = vec4(v.xyz, 1.);

    }
}