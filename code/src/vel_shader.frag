#version 330 core
// no need to specify layout when there is only one output
out vec4 diffuseColor;

uniform vec4 color;

void main() {
    diffuseColor = color;
}
