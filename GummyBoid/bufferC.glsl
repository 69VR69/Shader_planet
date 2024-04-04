void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = (2.0 * fragCoord.xy - iResolution.xy) / iResolution.y;

    // Time varying pixel color
    vec3 col = vec3(uv, 1.0);

    // Output to screen
    fragColor = vec4(col, 1.0);
}