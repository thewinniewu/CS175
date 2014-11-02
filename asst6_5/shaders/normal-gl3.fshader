#version 150

uniform sampler2D uTexColor;
uniform sampler2D uTexNormal;

// lights in eye space
uniform vec3 uLight;
uniform vec3 uLight2;

in vec2 vTexCoord;
in mat3 vNTMat;
in vec3 vEyePos;

out vec4 fragColor;

void main() {
  // get normals from texture image 
  vec4 normal_vec4 = texture(uTexNormal, vTexCoord);
  
  // scale and shift
  vec3 shifted_normal = vec3(2, 2, 2) * vec3(normal_vec4.x, normal_vec4.y, normal_vec4.z) - vec3(1, 1, 1);  
  
  // get normal 
  vec3 normal = normalize(vNTMat * shifted_normal);

  vec3 viewDir = normalize(-vEyePos);
  vec3 lightDir = normalize(uLight - vEyePos);
  vec3 lightDir2 = normalize(uLight2 - vEyePos);

  float nDotL = dot(normal, lightDir);
  vec3 reflection = normalize( 2.0 * normal * nDotL - lightDir);
  float rDotV = max(0.0, dot(reflection, viewDir));
  float specular = pow(rDotV, 32.0);
  float diffuse = max(nDotL, 0.0);

  nDotL = dot(normal, lightDir2);
  reflection = normalize( 2.0 * normal * nDotL - lightDir2);
  rDotV = max(0.0, dot(reflection, viewDir));
  specular += pow(rDotV, 32.0);
  diffuse += max(nDotL, 0.0);

  vec3 color = texture(uTexColor, vTexCoord).xyz * diffuse + specular * vec3(0.6, 0.6, 0.6);

  fragColor = vec4(color, 1);
}
