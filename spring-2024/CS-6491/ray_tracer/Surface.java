class Surface {
    final Color diffuse;
    final Color specular; // Optional specular color
    final float specPower; // Optional: Specular power (higher values result in tighter highlights)
    final float reflectivity; // Optional: Reflectivity coefficient (0 to 1, how reflective the surface is)
    final float glossRadius; // Optional: Gloss radius (0 for perfect mirror, higher for glossier appearance)

    // Glossy surface
    Surface(Color diffuse, Color specular, float specPower, float reflectivity, float glossRadius) {
        this.diffuse = diffuse;
        this.specular = specular;
        this.specPower = specPower;
        this.reflectivity = reflectivity;
        this.glossRadius = glossRadius;
    }

    // Diffuse-only
    Surface(Color diffuse) {
        this(diffuse, null, 0, 0, 0);
    }
}
