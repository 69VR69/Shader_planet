void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = (2.0 * fragCoord.xy - iResolution.xy) / iResolution.y;

    // Output to screen
    fragColor = vec4(uv, 1., 1.);
}