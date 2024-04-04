//#iChannel0 "file://bufferC.glsl"

#define NB_PARTICLES 100
#define PARTICLE_SIZE 0.01
#define OPACITY 0.5

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    // Normalized pixel coordinates (from 0 to 1)
    //vec2 uv = (2.0 * fragCoord.xy - iResolution.xy) / iResolution.y;

    vec2 pixelSize = .9 / iResolution.xy;
    vec2 uv = fragCoord.xy * pixelSize;
    vec3 col = vec3(uv, 1.0);

    for(int x = -1; x < NB_PARTICLES; x++)
    {
        for(int y = -1; y < NB_PARTICLES; y++)
        {
             // This is the bottleneck of the shader, there might be a
             // better way to read the particle textures.
            vec4 currentParticle = texture(iChannel0, vec2(x, y) * pixelSize);
            vec2 particlePixelVector = currentParticle.xy - uv;

            // If a particle is close to this pixel, add its color to the final color.
            if(particlePixelVector.x * particlePixelVector.x + particlePixelVector.y * particlePixelVector.y < pixelSize.x * PARTICLE_SIZE)
            {
                vec3 velocityColor = vec3(0.882, 0.875, 0.875) * OPACITY;
                col += velocityColor.xyz;
            }
        }
    }

    // Output to screen
    fragColor = vec4(col, 1.0);
}