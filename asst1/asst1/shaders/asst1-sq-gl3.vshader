#version 130

uniform float uVertexScale;
uniform float uXScale;
uniform float uYScale;

float xPos;
float yPos;

in vec2 aPosition;
in vec2 aTexCoord;

out vec2 vTexCoord;
out vec2 vTemp;

void main() {
  yScale = uVertexScale;
  gl_Position = vec4(aPosition.x * uVertexScale * yScale, aPosition.y * yScale, 0, 1);
  xPos = aPosition.x;
  yPos = aPosition.y;
  vTexCoord = aTexCoord;
  vTemp = vec2(1, 1);
}
