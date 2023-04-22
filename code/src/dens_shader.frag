#version 330 core
// no need to specify layout when there is only one output
out vec4 diffuseColor;

in float color;

void main() {
    diffuseColor = vec4(color, color, color, 1.0);
}
