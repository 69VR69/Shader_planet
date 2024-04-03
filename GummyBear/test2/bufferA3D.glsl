#iChannel0 "file://bufferA.glsl"

#define numboids 20.0 // Number of boids (must be integer value represented as float)
#define speed 10.0    // Boid speed
#define a1 0.001      // Collision factor
#define a2 0.01       // Cohesion factor
#define a3 0.05       // Alignment factor

vec4 getBoid(float id) {
    return texture(iChannel0, vec2((id + 0.5) / numboids, 0.5));
}

vec3 hash(float n) {
    return fract(sin(vec3(n, n * 7.0, n * 13.0)) * 43758.5453);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    fragColor = vec4(0.0); // Default data

    if (fragCoord.y > 0.5 || fragCoord.x > numboids)
        discard;

    float id = floor(fragCoord.x);

    // Initialize random boid positions and velocities
    if (iFrame < 5) {
        vec3 init_p = hash(id) * iResolution.xyz;
        vec3 init_v = normalize(hash(id + 0.1) * 2.0 - 1.0); // Random velocity direction
        fragColor = vec4(init_p, init_v);
    } else {
        // Boid of interest
        vec4 boid = getBoid(id);
        vec3 p = boid.xyz; // Boid position
        vec3 v = vec3(boid.w * speed); // Boid velocity

        // Initialize velocity update vectors
        vec3 v_collision = vec3(0.0);
        vec3 v_cohesion = vec3(0.0);
        vec3 v_alignment = vec3(0.0);

        // Loop over neighboring boids and build update vectors accordingly
        for (float nid = 0.0; nid < numboids; nid++) {
            if (nid != id) {
                vec4 neighbor = getBoid(nid);
                vec3 p_neighbor = neighbor.xyz; // Neighbor boid position
                vec3 v_neighbor = vec3(neighbor.w * speed); // Neighbor boid velocity

                vec3 sep = p_neighbor - p;
                float ls = length(sep);
                float r = 50.0; // Separation radius

                // Collision Update
                if (ls < r) {
                    v_collision -= sep;
                }

                // Cohesion Update
                v_cohesion += sep;

                // Alignment Update
                v_alignment += v_neighbor;
            }
        }

        // Perceived values
        v_cohesion /= (numboids - 1.0);
        v_alignment = (v_alignment / (numboids - 1.0)) - v;

        // Apply position and velocity updates
        v += a1 * v_collision + a2 * v_cohesion + a3 * v_alignment;
        v = normalize(v) * speed; // Constrain velocity
        p += iTimeDelta * v; // Update position

        // Impose wrapped boundary conditions
        p = mod(p, iResolution.xyz);

        fragColor = vec4(p.xy, v.xy);
    }
}
