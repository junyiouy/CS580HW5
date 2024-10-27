/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include    <algorithm>
#include     <cmath>


#define PI (float) 3.14159265358979323846

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif


static short matlevelNormal;

/* Helper Functions */
void normalize(GzCoord& vec) {
    float length = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (length == 0) return;
    vec[0] /= length;
    vec[1] /= length;
    vec[2] /= length;
}

float dotProduct(const GzCoord& vec1, const GzCoord& vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

void crossProduct(const GzCoord& vec1, const GzCoord& vec2, GzCoord& result) {
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

void multiplyMatrix(const GzMatrix mat1, const GzMatrix mat2, GzMatrix result) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < 4; ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

#define ARRAY(x, y) ((x) + (y) * xres)
#define ctoi(c) ((int)((c) * 4095))


class EdgeDDA {
    // Class to represent an edge with DDA algorithm
public:
    float x = 0.0f, z = 0.0f;
    float u = 0.0f, v = 0.0f; // Add UV interpolation
    float dx = 0.0f, dz = 0.0f;
    float du = 0.0f, dv = 0.0f; // Add slopes for UV interpolation
    int yStart = 0, yEnd = 0;

    static EdgeDDA setupEdgeDDA(const GzCoord& v1, const GzCoord& v2, const GzTextureIndex& uv1, const GzTextureIndex& uv2) {
        EdgeDDA dda;
        dda.yStart = std::ceil(v1[1]);
        dda.yEnd = std::ceil(v2[1]);

        float dy = v2[1] - v1[1];
        dda.dx = (v2[0] - v1[0]) / dy;
        dda.dz = (v2[2] - v1[2]) / dy;

        // Set up UV interpolation slopes
        dda.du = (uv2[0] - uv1[0]) / dy;
        dda.dv = (uv2[1] - uv1[1]) / dy;

        // Set initial values for x, z, u, and v
        dda.x = v1[0] + dda.dx * (dda.yStart - v1[1]);  // x position for ceiled yStart
        dda.z = v1[2] + dda.dz * (dda.yStart - v1[1]);  // z position for ceiled yStart
        dda.u = uv1[0] + dda.du * (dda.yStart - v1[1]); // u position for ceiled yStart
        dda.v = uv1[1] + dda.dv * (dda.yStart - v1[1]); // v position for ceiled yStart

        return dda;
    }
};



int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
    /* HW 3.1
    // Create rotate matrix : rotate along x axis
    // Pass back the matrix using mat value
    */
    float radians = degree * PI / 180.0f;
    mat[0][0] = 1.0f; mat[0][1] = 0.0f;          mat[0][2] = 0.0f;           mat[0][3] = 0.0f;
    mat[1][0] = 0.0f; mat[1][1] = cos(radians);  mat[1][2] = -sin(radians);  mat[1][3] = 0.0f;
    mat[2][0] = 0.0f; mat[2][1] = sin(radians);  mat[2][2] = cos(radians);   mat[2][3] = 0.0f;
    mat[3][0] = 0.0f; mat[3][1] = 0.0f;          mat[3][2] = 0.0f;           mat[3][3] = 1.0f;
    return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
    /* HW 3.2
    // Create rotate matrix : rotate along y axis
    // Pass back the matrix using mat value
    */
    float radians = degree * PI / 180.0f;
    mat[0][0] = cos(radians);  mat[0][1] = 0.0f; mat[0][2] = sin(radians); mat[0][3] = 0.0f;
    mat[1][0] = 0.0f;          mat[1][1] = 1.0f; mat[1][2] = 0.0f;         mat[1][3] = 0.0f;
    mat[2][0] = -sin(radians); mat[2][1] = 0.0f; mat[2][2] = cos(radians); mat[2][3] = 0.0f;
    mat[3][0] = 0.0f;          mat[3][1] = 0.0f; mat[3][2] = 0.0f;         mat[3][3] = 1.0f;
    return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
    /* HW 3.3
    // Create rotate matrix : rotate along z axis
    // Pass back the matrix using mat value
    */
    float radians = degree * PI / 180.0f;
    mat[0][0] = cos(radians); mat[0][1] = -sin(radians); mat[0][2] = 0.0f; mat[0][3] = 0.0f;
    mat[1][0] = sin(radians); mat[1][1] = cos(radians);  mat[1][2] = 0.0f; mat[1][3] = 0.0f;
    mat[2][0] = 0.0f;         mat[2][1] = 0.0f;          mat[2][2] = 1.0f; mat[2][3] = 0.0f;
    mat[3][0] = 0.0f;         mat[3][1] = 0.0f;          mat[3][2] = 0.0f; mat[3][3] = 1.0f;
    return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
    /* HW 3.4
    // Create translation matrix
    // Pass back the matrix using mat value
    */
    mat[0][0] = 1.0f; mat[0][1] = 0.0f; mat[0][2] = 0.0f; mat[0][3] = translate[0];
    mat[1][0] = 0.0f; mat[1][1] = 1.0f; mat[1][2] = 0.0f; mat[1][3] = translate[1];
    mat[2][0] = 0.0f; mat[2][1] = 0.0f; mat[2][2] = 1.0f; mat[2][3] = translate[2];
    mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;
    return GZ_SUCCESS;
}

int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
    /* HW 3.5
    // Create scaling matrix
    // Pass back the matrix using mat value
    */
    mat[0][0] = scale[0]; mat[0][1] = 0.0f;     mat[0][2] = 0.0f;     mat[0][3] = 0.0f;
    mat[1][0] = 0.0f;     mat[1][1] = scale[1]; mat[1][2] = 0.0f;     mat[1][3] = 0.0f;
    mat[2][0] = 0.0f;     mat[2][1] = 0.0f;     mat[2][2] = scale[2]; mat[2][3] = 0.0f;
    mat[3][0] = 0.0f;     mat[3][1] = 0.0f;     mat[3][2] = 0.0f;     mat[3][3] = 1.0f;
    return GZ_SUCCESS;
}

GzRender::GzRender(int xRes, int yRes)
{
    /* HW1.1 create a framebuffer for MS Windows display:
     -- set display resolution
     -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
     -- allocate memory for pixel buffer
     */
     // Set display resolution
    xres = xRes;
    yres = yRes;

    // Allocate memory for framebuffer
    framebuffer = new char[xRes * yRes * 3];

    // Allocate memory for pixel buffer
    pixelbuffer = new GzPixel[xRes * yRes];

    // Initialize pixel buffer
    GzDefault();

    /* HW 3.6
    - setup Xsp and anything only done once
    - init default camera
    */

    // setup Xsp
    Xsp[0][0] = (float)xres / 2.0f;
    Xsp[0][1] = 0.0f;
    Xsp[0][2] = 0.0f;
    Xsp[0][3] = (float)xres / 2.0f;

    Xsp[1][0] = 0.0f;
    Xsp[1][1] = -(float)yres / 2.0f;
    Xsp[1][2] = 0.0f;
    Xsp[1][3] = (float)yres / 2.0f;

    Xsp[2][0] = 0.0f;
    Xsp[2][1] = 0.0f;
    Xsp[2][2] = INT_MAX;
    Xsp[2][3] = 0.0f;

    Xsp[3][0] = 0.0f;
    Xsp[3][1] = 0.0f;
    Xsp[3][2] = 0.0f;
    Xsp[3][3] = 1.0f;

    // Initialize default camera
    m_camera.position[0] = DEFAULT_IM_X;
    m_camera.position[1] = DEFAULT_IM_Y;
    m_camera.position[2] = DEFAULT_IM_Z;

    m_camera.lookat[0] = 0.0f;
    m_camera.lookat[1] = 0.0f;
    m_camera.lookat[2] = 0.0f;

    m_camera.worldup[0] = 0.0f;
    m_camera.worldup[1] = 1.0f;
    m_camera.worldup[2] = 0.0f;

    m_camera.FOV = 90.0f;

    matlevel = 0; // Initialize the matrix stack level
    numlights = 0; // Initialize the number of lights
    interp_mode = GZ_FLAT; // Initialize the interpolation mode




}

GzRender::~GzRender()
{
    // Free buffer memory
    delete[] framebuffer;
    delete[] pixelbuffer;
}

int GzRender::GzDefault()
{
    // Set pixel buffer to default values
    for (int i = 0; i < xres * yres; ++i) {
        pixelbuffer[i] = { 4095, 4095, 4095, 1, INT_MAX }; // White for 12-bit color
    }
    return GZ_SUCCESS;
}
int GzRender::GzBeginRender()
{
    /* HW 3.7
    - setup for start of each frame - init frame buffer color,alpha,z
    - compute Xiw and projection xform Xpi from camera definition
    - init Ximage - put Xsp at base of stack, push on Xpi and Xiw
    - now stack contains Xsw and app can push model Xforms when needed
    */
    // Initialize the frame buffer for a new frame
    GzDefault();

    GzMatrix Xpi;
    GzMatrix Xiw;

    // Compute Xpi
    float d = 1.0f / tanf((m_camera.FOV / 2.0f) * (PI / 180.0f)); // Convert FOV to radians
    Xpi[0][0] = 1.0f; Xpi[0][1] = 0.0f; Xpi[0][2] = 0.0f;    Xpi[0][3] = 0.0f;
    Xpi[1][0] = 0.0f; Xpi[1][1] = 1.0f; Xpi[1][2] = 0.0f;    Xpi[1][3] = 0.0f;
    Xpi[2][0] = 0.0f; Xpi[2][1] = 0.0f; Xpi[2][2] = 1.0f / d;    Xpi[2][3] = 0.0f;
    Xpi[3][0] = 0.0f; Xpi[3][1] = 0.0f; Xpi[3][2] = 1.0f / d; Xpi[3][3] = 1.0f;

    // Compute Xiw
    GzCoord cl = {
        m_camera.lookat[0] - m_camera.position[0],
        m_camera.lookat[1] - m_camera.position[1],
        m_camera.lookat[2] - m_camera.position[2]
    };
    normalize(cl);

    GzCoord upPrime;
    float dotUpCl = dotProduct(m_camera.worldup, cl);
    upPrime[0] = m_camera.worldup[0] - dotUpCl * cl[0];
    upPrime[1] = m_camera.worldup[1] - dotUpCl * cl[1];
    upPrime[2] = m_camera.worldup[2] - dotUpCl * cl[2];
    normalize(upPrime);

    GzCoord cameraX;
    crossProduct(upPrime, cl, cameraX);

    Xiw[0][0] = cameraX[0]; Xiw[0][1] = cameraX[1]; Xiw[0][2] = cameraX[2]; Xiw[0][3] = -dotProduct(cameraX, m_camera.position);
    Xiw[1][0] = upPrime[0]; Xiw[1][1] = upPrime[1]; Xiw[1][2] = upPrime[2]; Xiw[1][3] = -dotProduct(upPrime, m_camera.position);
    Xiw[2][0] = cl[0];      Xiw[2][1] = cl[1];      Xiw[2][2] = cl[2];      Xiw[2][3] = -dotProduct(cl, m_camera.position);
    Xiw[3][0] = 0.0f;       Xiw[3][1] = 0.0f;       Xiw[3][2] = 0.0f;       Xiw[3][3] = 1.0f;

    // Initialize the image matrix stack
    int status = 0;

    matlevel = 0;
    matlevelNormal = 0;
    status |= GzPushMatrix(Xsp);
    status |= GzPushMatrix(Xpi);
    status |= GzPushMatrix(Xiw);


    if (status) return GZ_FAILURE;
    else return GZ_SUCCESS;
}


int GzRender::GzPutCamera(GzCamera camera)
{
    /* HW 3.8
    /*- overwrite renderer camera structure with new camera definition
    */

    // Overwrite the current camera structure with the new camera definition
    m_camera = camera;
    return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix matrix) {
    /* HW 3.9
    - push a matrix onto the Ximage stack
    - check for stack overflow
    */

    // Step 1: Push matrix onto Ximage stack
    if (matlevel >= MATLEVELS) {
        return GZ_FAILURE; // Stack overflow
    }

    if (matlevel == 0) {
        // Base case: first matrix on the stack
        memcpy(Ximage[matlevel], matrix, sizeof(GzMatrix));
    }
    else {
        // Multiply top of the stack with the new matrix
        GzMatrix result;
        multiplyMatrix(Ximage[matlevel - 1], matrix, result);
        memcpy(Ximage[matlevel], result, sizeof(GzMatrix));
    }
    matlevel++;

    // Step 2: Push matrix onto Xnorm stack
    GzMatrix identityMatrix;
    memset(identityMatrix, 0, sizeof(GzMatrix));
    for (int i = 0; i < 4; ++i) {
        identityMatrix[i][i] = 1.0f;
    }

    if (matlevelNormal < 2) {
        // For Xsp and Xpi, push identity matrix
        memcpy(Xnorm[matlevelNormal], identityMatrix, sizeof(GzMatrix));
    }

    else {
        // Normalize each row of the rotation matrix
        for (int i = 0; i < 3; ++i) {
            float length = sqrt(matrix[i][0] * matrix[i][0] +
                matrix[i][1] * matrix[i][1] +
                matrix[i][2] * matrix[i][2]);

            // Normalize the row if length is not zero
            if (length != 0) {
                matrix[i][0] /= length;
                matrix[i][1] /= length;
                matrix[i][2] /= length;
            }
        }
        // Zero out the translation part
        matrix[0][3] = 0;
        matrix[1][3] = 0;
        matrix[2][3] = 0;

        // Push the normalized matrix onto the stack
        GzMatrix result;
        multiplyMatrix(Xnorm[matlevelNormal - 1], matrix, result);
        memcpy(Xnorm[matlevelNormal], result, sizeof(GzMatrix));
    }

    matlevelNormal++;

    return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
    /* HW 3.10
    - pop a matrix off the Ximage stack
    - check for stack underflow
    */
    if (matlevel <= 0) {
        return GZ_FAILURE; // Stack underflow
    }

    matlevel--;
    return GZ_SUCCESS;
}


int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
    // Check bounds and values
    if (i < 0 || i >= xres || j < 0 || j >= yres ||
        a < 0 || a > 1 ||
        r < 0 || g < 0 || b < 0 ||
        z < 0 || z > INT_MAX)
    {
        return GZ_FAILURE;
    }

    // Clamp values because values above 4095 are not valid
    r = r > 4095 ? 4095 : r;
    g = g > 4095 ? 4095 : g;
    b = b > 4095 ? 4095 : b;

    // Write pixel values into the buffer
    pixelbuffer[i + j * xres] = { r, g, b, a, z };

    return GZ_SUCCESS;
}

int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
    // Retrieve a pixel information from the pixel buffer
    *r = pixelbuffer[i + j * xres].red;
    *g = pixelbuffer[i + j * xres].green;
    *b = pixelbuffer[i + j * xres].blue;
    *a = pixelbuffer[i + j * xres].alpha;
    *z = pixelbuffer[i + j * xres].z;

    return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2File(FILE* outfile)
{
    // Write image to ppm file
    fprintf(outfile, "P6 %d %d 255\n", xres, yres);
    for (int i = 0; i < xres * yres; ++i) {
        char colors[3];
        colors[0] = pixelbuffer[i].red >> 4; // 12-bit to 8-bit
        colors[1] = pixelbuffer[i].green >> 4;
        colors[2] = pixelbuffer[i].blue >> 4;
        fwrite(colors, 1, 3, outfile);
    }
    return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
    // Write pixels to framebuffer: store as blue, green, red
    for (int i = 0; i < xres * yres; ++i) {
        int index = i * 3;
        framebuffer[index] = pixelbuffer[i].blue >> 4;
        framebuffer[index + 1] = pixelbuffer[i].green >> 4;
        framebuffer[index + 2] = pixelbuffer[i].red >> 4;
    }
    return GZ_SUCCESS;
}

int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList) {
    for (int i = 0; i < numAttributes; ++i) {
        switch (nameList[i]) {
        case GZ_RGB_COLOR: {
            // Set default flat shading color
            GzColor* color = static_cast<GzColor*>(valueList[i]);
            flatcolor[RED] = (*color)[RED];
            flatcolor[GREEN] = (*color)[GREEN];
            flatcolor[BLUE] = (*color)[BLUE];
            break;
        }
        case GZ_INTERPOLATE: {
            // Set interpolation mode (flat, Gouraud, Phong)
            interp_mode = *(int*)valueList[i];
            break;
        }
        case GZ_AMBIENT_LIGHT: {
            // Set ambient light parameters
            ambientlight = *(GzLight*)valueList[i];
            break;
        }
        case GZ_DIRECTIONAL_LIGHT: {
            // Set directional light (up to 10 lights)
            if (numlights < MAX_LIGHTS) {
                lights[numlights] = *(GzLight*)valueList[i];
                numlights++;
            }
            break;
        }
        case GZ_AMBIENT_COEFFICIENT: {
            // Set ambient reflection coefficient (Ka)
            GzColor* Ka = static_cast<GzColor*>(valueList[i]);
            this->Ka[RED] = (*Ka)[RED];
            this->Ka[GREEN] = (*Ka)[GREEN];
            this->Ka[BLUE] = (*Ka)[BLUE];
            break;
        }
        case GZ_DIFFUSE_COEFFICIENT: {
            // Set diffuse reflection coefficient (Kd)
            GzColor* Kd = static_cast<GzColor*>(valueList[i]);
            this->Kd[RED] = (*Kd)[RED];
            this->Kd[GREEN] = (*Kd)[GREEN];
            this->Kd[BLUE] = (*Kd)[BLUE];
            break;
        }
        case GZ_SPECULAR_COEFFICIENT: {
            // Set specular reflection coefficient (Ks)
            GzColor* Ks = static_cast<GzColor*>(valueList[i]);
            this->Ks[RED] = (*Ks)[RED];
            this->Ks[GREEN] = (*Ks)[GREEN];
            this->Ks[BLUE] = (*Ks)[BLUE];
            break;
        }
        case GZ_DISTRIBUTION_COEFFICIENT: {
            // Set specular power (shininess)
            spec = *(float*)valueList[i];
            break;
        }
		case GZ_TEXTURE_MAP: {
			// Set texture function pointer
            tex_fun = (GzTexture)valueList[i];
			break;

        }
    }
    return GZ_SUCCESS;
}


// Interpolate plane equation based on vertex attributes
void InterpolatePlane(GzCoord* vertices, GzColor* attributes, int index, float& A, float& B, float& C, float& D) {
    float X1 = vertices[1][0] - vertices[0][0];
    float Y1 = vertices[1][1] - vertices[0][1];
    float Z1 = attributes[1][index] - attributes[0][index];
    float X2 = vertices[2][0] - vertices[0][0];
    float Y2 = vertices[2][1] - vertices[0][1];
    float Z2 = attributes[2][index] - attributes[0][index];
    A = (Y1 * Z2) - (Z1 * Y2);
    B = -((X1 * Z2) - (Z1 * X2));
    C = (X1 * Y2) - (Y1 * X2);
    D = -1.0f * (A * vertices[0][0] + B * vertices[0][1] + C * attributes[0][index]);
}

// Interpolate color using plane equation
float InterpolateColor(int x, int y, float A, float B, float C, float D) {
    return -1.0f * (A * x + B * y + D) / C;
}
// Transform vertices and normals
void TransformVerticesAndNormals(const GzCoord* inputVertices, const GzCoord* inputNormals, GzCoord* outputVertices, GzCoord* outputNormals, GzMatrix* matrixStack, GzMatrix* normalMatrix, int matrixLevel) {
    float vertices4D[3][4], normals4D[3][4];
    float transformedVertices4D[3][4], transformedNormals4D[3][4];

    for (int v = 0; v < 3; ++v) {
        // Convert vertices and normals to 4D
        for (int i = 0; i < 3; ++i) {
            vertices4D[v][i] = inputVertices[v][i];
            normals4D[v][i] = inputNormals != nullptr ? inputNormals[v][i] : 0.0f;
        }
        vertices4D[v][3] = 1.0f;
        normals4D[v][3] = 1.0f;

        // Apply transformation matrix
        for (int j = 0; j < 4; ++j) {
            transformedVertices4D[v][j] = 0;
            transformedNormals4D[v][j] = 0;
            for (int i = 0; i < 4; ++i) {
                transformedVertices4D[v][j] += matrixStack[matrixLevel - 1][j][i] * vertices4D[v][i];
                transformedNormals4D[v][j] += normalMatrix[matrixLevel - 1][j][i] * normals4D[v][i];
            }
        }
    }

    // Convert back to 3D
    for (int v = 0; v < 3; ++v) {
        for (int i = 0; i < 3; ++i) {
            outputVertices[v][i] = transformedVertices4D[v][i] / transformedVertices4D[v][3];
            outputNormals[v][i] = transformedNormals4D[v][i] / transformedNormals4D[v][3];
        }
    }
}

// Calculate edge coefficients
void CalculateEdgeCoefficients(const GzCoord* inputVertices, float& coeffA12, float& coeffB12, float& coeffC12, float& coeffA23, float& coeffB23, float& coeffC23, float& coeffA31, float& coeffB31, float& coeffC31) {
    float deltaX12 = inputVertices[1][0] - inputVertices[0][0];
    float deltaY12 = inputVertices[1][1] - inputVertices[0][1];
    float deltaX23 = inputVertices[2][0] - inputVertices[1][0];
    float deltaY23 = inputVertices[2][1] - inputVertices[1][1];
    float deltaX31 = inputVertices[0][0] - inputVertices[2][0];
    float deltaY31 = inputVertices[0][1] - inputVertices[2][1];

    coeffA12 = deltaY12;
    coeffB12 = -deltaX12;
    coeffC12 = deltaX12 * inputVertices[0][1] - deltaY12 * inputVertices[0][0];

    coeffA23 = deltaY23;
    coeffB23 = -deltaX23;
    coeffC23 = deltaX23 * inputVertices[1][1] - deltaY23 * inputVertices[1][0];

    coeffA31 = deltaY31;
    coeffB31 = -deltaX31;
    coeffC31 = deltaX31 * inputVertices[2][1] - deltaY31 * inputVertices[2][0];
}

// Compute Z plane
void ComputeZPlane(const GzCoord* inputVertices, float& planeCoeffA, float& planeCoeffB, float& planeCoeffC, float& planeCoeffD) {
    float diffX1 = inputVertices[1][0] - inputVertices[0][0];
    float diffY1 = inputVertices[1][1] - inputVertices[0][1];
    float diffZ1 = inputVertices[1][2] - inputVertices[0][2];
    float diffX2 = inputVertices[2][0] - inputVertices[0][0];
    float diffY2 = inputVertices[2][1] - inputVertices[0][1];
    float diffZ2 = inputVertices[2][2] - inputVertices[0][2];

    planeCoeffA = (diffY1 * diffZ2) - (diffZ1 * diffY2);
    planeCoeffB = -((diffX1 * diffZ2) - (diffZ1 * diffX2));
    planeCoeffC = (diffX1 * diffY2) - (diffY1 * diffX2);
    planeCoeffD = -1.0f * (planeCoeffA * inputVertices[0][0] + planeCoeffB * inputVertices[0][1] + planeCoeffC * inputVertices[0][2]);
}

// Interpolating colors using Phong model based on interpolated normals and lighting
void InterpolateColorsAndNormals(const GzCoord* inputNormals, GzColor* interpolatedColor, const GzLight* lightSources, int numLights, const GzColor& Ks, const GzColor& Kd, const GzColor& Ka, const GzLight& ambientLight, float shininessFactor) {
    for (int v = 0; v < 3; ++v) {
        GzColor specColor = { 0 }, diffColor = { 0 }, ambientColor = { 0 };
        GzCoord viewDir = { 0, 0, -1 }; // Assuming camera is at the origin looking down -Z axis

        // Loop through all the light sources
        for (int i = 0; i < numLights; ++i) {
            GzCoord lightDir = { lightSources[i].direction[0], lightSources[i].direction[1], lightSources[i].direction[2] };
            float nDotL = dotProduct(inputNormals[v], lightDir);
            float nDotV = dotProduct(inputNormals[v], viewDir);

            if (nDotL * nDotV > 0) {
                // If both are positive or both are negative, proceed with shading
                GzCoord effectiveNormal;
                effectiveNormal[0] = inputNormals[v][0];
                effectiveNormal[1] = inputNormals[v][1];
                effectiveNormal[2] = inputNormals[v][2];

                if (nDotL < 0 && nDotV < 0) {
                    // Flip the normal if both dot products are negative
                    effectiveNormal[0] = -inputNormals[v][0];
                    effectiveNormal[1] = -inputNormals[v][1];
                    effectiveNormal[2] = -inputNormals[v][2];
                    nDotL = -nDotL;
                    nDotV = -nDotV;
                }

                // Calculate the reflection direction
                GzCoord reflectionDir;
                for (int j = 0; j < 3; ++j) {
                    reflectionDir[j] = 2 * nDotL * effectiveNormal[j] - lightDir[j];
                }
                normalize(reflectionDir);

                float reflectionDotView = max(0.0f, dotProduct(reflectionDir, viewDir));
                float specIntensity = pow(reflectionDotView, shininessFactor);

                // Calculate specular and diffuse contributions
                for (int j = 0; j < 3; ++j) {
                    specColor[j] += Ks[j] * specIntensity * lightSources[i].color[j];
                    diffColor[j] += Kd[j] * nDotL * lightSources[i].color[j];
                }
            }
        }

        // Add ambient color contribution
        for (int j = 0; j < 3; ++j) {
            ambientColor[j] += Ka[j] * ambientLight.color[j];
            interpolatedColor[v][j] = specColor[j] + diffColor[j] + ambientColor[j];
            interpolatedColor[v][j] = min(1.0f, max(0.0f, interpolatedColor[v][j])); // Clamp values between 0 and 1
        }
    }
}




// Calculate Phong shading based on interpolated normals and lighting
void ComputePhongShading(const GzCoord& normal, GzColor& resultColor, const GzLight* lightSources, int numLights, const GzColor& Ks, const GzColor& Kd, const GzColor& Ka, const GzLight& ambientLight, float shininess) {
    GzCoord viewDir = { 0, 0, -1 }; // Viewing direction

    GzColor specular = { 0 }, diffuse = { 0 }, ambient = { 0 };

    // Loop through all lights
    for (int i = 0; i < numLights; ++i) {
        GzCoord lightDir = { lightSources[i].direction[0], lightSources[i].direction[1], lightSources[i].direction[2] };
        float nDotL = normal[0] * lightDir[0] + normal[1] * lightDir[1] + normal[2] * lightDir[2];
        float nDotV = normal[0] * viewDir[0] + normal[1] * viewDir[1] + normal[2] * viewDir[2];

        // Ensure light and view directions are on the same side
        if (nDotL * nDotV > 0) {
            GzCoord reflectDir;
            for (int j = 0; j < 3; ++j) {
                reflectDir[j] = 2 * nDotL * normal[j] - lightDir[j];
            }
            normalize(reflectDir);

            float rDotV = max(0.0f, reflectDir[0] * viewDir[0] + reflectDir[1] * viewDir[1] + reflectDir[2] * viewDir[2]);
            float specularFactor = pow(rDotV, shininess);

            // Calculate specular and diffuse components
            for (int j = 0; j < 3; ++j) {
                specular[j] += Ks[j] * specularFactor * lightSources[i].color[j];
                if (nDotL > 0 && nDotV > 0) {
                    diffuse[j] += Kd[j] * nDotL * lightSources[i].color[j];
                }
                else {
                    diffuse[j] += Kd[j] * (-nDotL) * lightSources[i].color[j];
                }
            }
        }
    }

    // Combine the components with ambient lighting
    for (int j = 0; j < 3; ++j) {
        ambient[j] += Ka[j] * ambientLight.color[j];
        resultColor[j] = specular[j] + diffuse[j] + ambient[j];
        resultColor[j] = min(1.0f, max(0.0f, resultColor[j]));
    }
}

// Compute the surface normal of a triangle given three vertices
void ComputeSurfaceNormal(const GzCoord* vertices, GzCoord& surfaceNormal) {
    GzCoord vector1, vector2;

    // Compute two edge vectors of the triangle
    vector1[0] = vertices[1][0] - vertices[0][0];
    vector1[1] = vertices[1][1] - vertices[0][1];
    vector1[2] = vertices[1][2] - vertices[0][2];

    vector2[0] = vertices[2][0] - vertices[0][0];
    vector2[1] = vertices[2][1] - vertices[0][1];
    vector2[2] = vertices[2][2] - vertices[0][2];

    // Compute the cross product of the two edge vectors to get the surface normal
    surfaceNormal[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
    surfaceNormal[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
    surfaceNormal[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];

    // Normalize the resulting surface normal
    float length = sqrt(surfaceNormal[0] * surfaceNormal[0] +
        surfaceNormal[1] * surfaceNormal[1] +
        surfaceNormal[2] * surfaceNormal[2]);

    if (length != 0) {
        surfaceNormal[0] /= length;
        surfaceNormal[1] /= length;
        surfaceNormal[2] /= length;
    }
}


// 双线性插值函数
void BilinearInterpolation(float u, float v, GzColor* image, int xs, int ys, GzColor& color) {
    // 将 u 和 v 限制在 [0, 1] 范围内，避免超出边界
    u = min(max(u, 0.0f), 1.0f);
    v = min(max(v, 0.0f), 1.0f);

    // 计算纹理图像中的真实坐标
    float scaledU = u * (xs - 1);
    float scaledV = v * (ys - 1);

    // 获取四个周围像素的整数坐标
    int u0 = static_cast<int>(std::floor(scaledU));
    int u1 = min(u0 + 1, xs - 1);
    int v0 = static_cast<int>(std::floor(scaledV));
    int v1 = min(v0 + 1, ys - 1);

    // 获取双线性插值的四个颜色值
    GzColor color00 = { image[u0 + v0 * xs][0], image[u0 + v0 * xs][1], image[u0 + v0 * xs][2] };
    GzColor color01 = { image[u0 + v1 * xs][0], image[u0 + v1 * xs][1], image[u0 + v1 * xs][2] };
    GzColor color10 = { image[u1 + v0 * xs][0], image[u1 + v0 * xs][1], image[u1 + v0 * xs][2] };
    GzColor color11 = { image[u1 + v1 * xs][0], image[u1 + v1 * xs][1], image[u1 + v1 * xs][2] };

    // 计算插值权重
    float s = scaledU - u0;
    float t = scaledV - v0;

    // 进行双线性插值
    for (int i = 0; i < 3; ++i) { // 对 RGB 三个通道分别进行插值
        color[i] = (1 - s) * (1 - t) * color00[i] +
            (1 - s) * t * color01[i] +
            s * (1 - t) * color10[i] +
            s * t * color11[i];
    }
}



int GzRender::GzPutTriangle(int numParts, GzToken* nameList, GzPointer* valueList) {
    GzCoord* vertices = nullptr;
    GzCoord* normals = nullptr;
    GzTextureIndex* uvList = nullptr; // Add variable to store texture indices

    // Extract triangle vertices, normals, and texture indices
    for (int i = 0; i < numParts; ++i) {
        switch (nameList[i]) {
        case GZ_POSITION:
            vertices = static_cast<GzCoord*>(valueList[i]);
            break;
        case GZ_NORMAL:
            normals = static_cast<GzCoord*>(valueList[i]);
            break;
        case GZ_TEXTURE_INDEX: // Handle texture indices
            uvList = static_cast<GzTextureIndex*>(valueList[i]);
            break;
        case GZ_NULL_TOKEN:
            break;
        }
    }

    if (vertices == nullptr) return GZ_FAILURE;

    GzCoord transformedVertices[3];
    GzCoord transformedNormals[3];

    // Step 1: Transform vertices and normals to camera space using global function
    TransformVerticesAndNormals(vertices, normals, transformedVertices, transformedNormals, Ximage, Xnorm, matlevel);


    // Sort vertices by Y coordinate
    if (transformedVertices[1][1] < transformedVertices[0][1]) {
        std::swap(transformedVertices[1], transformedVertices[0]);
        std::swap(transformedNormals[1], transformedNormals[0]);
    }
    if (transformedVertices[2][1] < transformedVertices[1][1]) {
        std::swap(transformedVertices[2], transformedVertices[1]);
        std::swap(transformedNormals[2], transformedNormals[1]);
    }
    if (transformedVertices[1][1] < transformedVertices[0][1]) {
        std::swap(transformedVertices[1], transformedVertices[0]);
        std::swap(transformedNormals[1], transformedNormals[0]);
    }

    // Step 2: Check if all vertices are behind the camera (z < 0)
    if (transformedVertices[0][2] < 0 && transformedVertices[1][2] < 0 && transformedVertices[2][2] < 0) {
        return GZ_SUCCESS;  // Skip this triangle
    }

    // Step 3: Convert each vertex's (u, v) to perspective space (U, V)
    GzTextureIndex perspectiveUV[3];
    float Vz[3]; // Store Vz for each vertex

    for (int i = 0; i < 3; ++i) {
        // Calculate Vz for each vertex using its depth
        Vz[i] = vertices[i][2] / (INT_MAX - vertices[i][2]);

        // Convert (u, v) to perspective space (U, V)
        perspectiveUV[i][0] = uvList[i][0] / (vertices[i][2] + Vz[i]);
        perspectiveUV[i][1] = uvList[i][1] / (vertices[i][2] + Vz[i]);
    }






    // Step 4: Calculate edge coefficients using global function
    float coeffA12, coeffB12, coeffC12, coeffA23, coeffB23, coeffC23, coeffA31, coeffB31, coeffC31;
    CalculateEdgeCoefficients(transformedVertices, coeffA12, coeffB12, coeffC12, coeffA23, coeffB23, coeffC23, coeffA31, coeffB31, coeffC31);

    // Step 5: Compute the Z-plane equation using global function
    float planeA, planeB, planeC, planeD;
    ComputeZPlane(transformedVertices, planeA, planeB, planeC, planeD);

    // Step 6: Interpolate colors and normals using global function
    GzColor finalColor[3];
    InterpolateColorsAndNormals(transformedNormals, finalColor, lights, numlights, Ks, Kd, Ka, ambientlight, spec);

    // Define the interpolation planes for red, green, blue, and normals
    float redPlane[4], greenPlane[4], bluePlane[4];
    InterpolatePlane(transformedVertices, finalColor, 0, redPlane[0], redPlane[1], redPlane[2], redPlane[3]);
    InterpolatePlane(transformedVertices, finalColor, 1, greenPlane[0], greenPlane[1], greenPlane[2], greenPlane[3]);
    InterpolatePlane(transformedVertices, finalColor, 2, bluePlane[0], bluePlane[1], bluePlane[2], bluePlane[3]);

    float normalXPlane[4], normalYPlane[4], normalZPlane[4];
    InterpolatePlane(transformedVertices, transformedNormals, 0, normalXPlane[0], normalXPlane[1], normalXPlane[2], normalXPlane[3]);
    InterpolatePlane(transformedVertices, transformedNormals, 1, normalYPlane[0], normalYPlane[1], normalYPlane[2], normalYPlane[3]);
    InterpolatePlane(transformedVertices, transformedNormals, 2, normalZPlane[0], normalZPlane[1], normalZPlane[2], normalZPlane[3]);

    // Step 7: Scanline rasterization
    if (transformedVertices[1][1] < transformedVertices[0][1]) std::swap(transformedVertices[1], transformedVertices[0]);
    if (transformedVertices[2][1] < transformedVertices[1][1]) std::swap(transformedVertices[2], transformedVertices[1]);
    if (transformedVertices[1][1] < transformedVertices[0][1]) std::swap(transformedVertices[1], transformedVertices[0]);

    // Setup DDA for each edge with texture coordinates
    EdgeDDA dda12 = EdgeDDA::setupEdgeDDA(transformedVertices[0], transformedVertices[1], perspectiveUV[0], perspectiveUV[1]);
    EdgeDDA dda23 = EdgeDDA::setupEdgeDDA(transformedVertices[1], transformedVertices[2], perspectiveUV[1], perspectiveUV[2]);
    EdgeDDA dda13 = EdgeDDA::setupEdgeDDA(transformedVertices[0], transformedVertices[2], perspectiveUV[0], perspectiveUV[2]);


    EdgeDDA* leftDDA;
    EdgeDDA* rightDDA;

    // Determine left and right edges
    float slope12 = (transformedVertices[1][0] - transformedVertices[0][0]) / (transformedVertices[1][1] - transformedVertices[0][1]);
    float slope13 = (transformedVertices[2][0] - transformedVertices[0][0]) / (transformedVertices[2][1] - transformedVertices[0][1]);

    int reverse = (slope12 < slope13) ? 0 : 1;

    // Scanline rasterization
    for (int y = dda12.yStart; y < dda23.yEnd; ++y) {
        if (y < dda23.yStart) {
            leftDDA = &dda12;
            rightDDA = &dda13;
        }
        else {
            leftDDA = &dda23;
            rightDDA = &dda13;
        }

        if (reverse) {
            std::swap(leftDDA, rightDDA);
        }

        int xStart = std::ceil(leftDDA->x);
        int xEnd = std::ceil(rightDDA->x);

        // Interpolate Z, red, green, blue colors along the scanline
        float z = leftDDA->z;
        float dz = (rightDDA->z - leftDDA->z) / (xEnd - xStart);

        float u = leftDDA->u;
        float du = (rightDDA->u - leftDDA->u) / (xEnd - xStart);

        float v = leftDDA->v;
        float dv = (rightDDA->v - leftDDA->v) / (xEnd - xStart);

        interp_mode = GZ_NORMALS;

        for (int x = xStart; x < xEnd; ++x) {
            if (x >= 0 && x < xres && y >= 0 && y < yres) {
                // Use interpolated Z, U, and V for further processing
                float interpolatedZ = z;
                float interpolatedU = u;
                float interpolatedV = v;

                // Apply perspective correction to get the final (u, v)
                float Vz_prime = interpolatedZ / (INT_MAX - interpolatedZ);
                //float finalU = interpolatedU * (Vz_prime + 1.0f);
                //float finalV = interpolatedV * (Vz_prime + 1.0f);
                float finalU = interpolatedU;
                float finalV = interpolatedV;

                // Calculate texture color using final (u, v)
                GzColor textureColor;
                BilinearInterpolation(finalU, finalV, textureImage, textureWidth, textureHeight, textureColor);

                GzIntensity redIntensity, greenIntensity, blueIntensity;

                if (interp_mode == GZ_FLAT) {
                    // Flat Shading
                    redIntensity = ctoi(finalColor[0][0]);
                    greenIntensity = ctoi(finalColor[0][1]);
                    blueIntensity = ctoi(finalColor[0][2]);
                }
                else if (interp_mode == GZ_COLOR) {
                    // Gouraud Shading
                    GzColor intensity;
                    intensity[0] = InterpolateColor(x, y, redPlane[0], redPlane[1], redPlane[2], redPlane[3]);
                    intensity[1] = InterpolateColor(x, y, greenPlane[0], greenPlane[1], greenPlane[2], greenPlane[3]);
                    intensity[2] = InterpolateColor(x, y, bluePlane[0], bluePlane[1], bluePlane[2], bluePlane[3]);

                    redIntensity = ctoi(intensity[0]);
                    greenIntensity = ctoi(intensity[1]);
                    blueIntensity = ctoi(intensity[2]);
                }
                else if (interp_mode == GZ_NORMALS) {
                    // Phong Shading
                    GzCoord interpolatedNormal;
                    interpolatedNormal[0] = InterpolateColor(x, y, normalXPlane[0], normalXPlane[1], normalXPlane[2], normalXPlane[3]);
                    interpolatedNormal[1] = InterpolateColor(x, y, normalYPlane[0], normalYPlane[1], normalYPlane[2], normalYPlane[3]);
                    interpolatedNormal[2] = InterpolateColor(x, y, normalZPlane[0], normalZPlane[1], normalZPlane[2], normalZPlane[3]);

                    GzColor intensity;
                    ComputePhongShading(interpolatedNormal, intensity, lights, numlights, Ks, Kd, Ka, ambientlight, spec);

                    redIntensity = ctoi(intensity[0]);
                    greenIntensity = ctoi(intensity[1]);
                    blueIntensity = ctoi(intensity[2]);
                }

                if (interpolatedZ < pixelbuffer[ARRAY(x, y)].z) {
                    GzPut(x, y, redIntensity, greenIntensity, blueIntensity, 1, interpolatedZ);
                }
            }
            z += dz;
        }

        // Update DDA for next scanline
        leftDDA->x += leftDDA->dx;
        rightDDA->x += rightDDA->dx;
        leftDDA->z += leftDDA->dz;
        rightDDA->z += rightDDA->dz;
    }

    return GZ_SUCCESS;
}
