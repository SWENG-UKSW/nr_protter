#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <vector>
#include <string>
#include <utility>

#include "CLI/CLI.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "plotter.h"
#include "overlay.h"

/////////////////////////////////////////////////////////////////////////

PlotData plot_data;

/////////////////////////////////////////////////////////////////////////

double wyznacznik(const double *macierz, int rozmiar) {
    if (rozmiar == 1) {
        return macierz[0];
    }
    if (rozmiar == 2) {
        return macierz[0] * macierz[3] - macierz[1] * macierz[2]; // (0*3 - 1*2)
    }

    double wynik = 0.0;


    // rozwiniêcie Laplace'a
    for (int j = 0; j < rozmiar; j++) {
        double *minor = (double*)malloc(sizeof(double)*(rozmiar-1)*(rozmiar-1));
        int minorIdx = 0;

        // tworzymy minor przez pomijanie wiersza 0 i kolumny j
        for (int i = 0; i < rozmiar; i++) {
            if (i == 0) continue; // pomijamy pierwszy wiersz
            for (int k = 0; k < rozmiar; k++) {
                if (k == j) continue; // pomijamy kolumnê j
                minor[minorIdx++] = macierz[i * rozmiar + k];
            }
        }

        int sign = (j % 2 == 0) ? 1 : -1;
        wynik += sign * macierz[j] * wyznacznik(minor, rozmiar - 1);
        free(minor);
    }

    return wynik;
}

///////////////////////////////////////////////////////////////////////////////
// matrix row_width x row_width-1
void zamienWiersze(double *A, double *b, int n, int wiersz1, int wiersz2) {
    for (int j = 0; j < n; j++) {
        double temp = A[wiersz1 * n + j];
        A[wiersz1 * n + j] = A[wiersz2 * n + j];
        A[wiersz2 * n + j] = temp;
    }
    double temp_b = b[wiersz1];
    b[wiersz1] = b[wiersz2];
    b[wiersz2] = temp_b;
}

// Funkcja eliminacji Gaussa z pivotacj¹
void gaussEliminationWithPivoting(double *A, double *b, int n) {
    for (int k = 0; k < n; k++) {
        // ZnajdŸ max element w kolumnie k od wiersza k w dó³
        int maxRow = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i * n + k]) > fabs(A[maxRow * n + k])) {
                maxRow = i;
            }
        }
        // Zamieñ wiersze jeœli trzeba
        if (maxRow != k) {
            zamienWiersze(A, b, n, k, maxRow);
        }
        // Eliminacja
        for (int i = k + 1; i < n; i++) {
            if (A[k * n + k] == 0) {
                printf("Dzielenie przez zero!\n");
                exit(1);
            }
            double c = A[i * n + k] / A[k * n + k];
            for (int j = k; j < n; j++) {
                A[i * n + j] -= c * A[k * n + j];
            }
            b[i] -= c * b[k];
        }
    }
}

// Funkcja rozwi¹zuj¹ca uk³ad przez podstawianie wstecz
void backSubstitution(double *A, double *b, double *x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i * n + j] * x[j];
        }
        if (A[i * n + i] == 0) {
            printf("Brak rozwi¹zañ lub nieskoñczenie wiele!\n");
            exit(1);
        }
        x[i] /= A[i * n + i];
    }
}



///////////////////////////////////////////////////////////////////////////////

// Funkcja do rozwi¹zywania uk³adu równañ metod¹ Jacobiego
void jacobi(double* A, double* b, double* x, int n) 
{
    std::vector<double> x_new (n, 0.0);

    constexpr int MAX_ITER= 10;
    constexpr double TOLERANCE = 1e-10;

    int iter = 0;
    while (iter < MAX_ITER) {
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    s += A[i * n + j] * x[j];
                }
            }
            x_new[i] = (b[i] - s) / A[i * n + i];
        }

        // Sprawdzenie kryterium zbie¿noœci
        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = fabs(x_new[i] - x[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
            x[i] = x_new[i];
        }

        if (max_diff < TOLERANCE) {
            break;
        }
        iter++;
    }
}

namespace {

const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexCoord;
out vec2 TexCoord;
uniform mat4 projection;
void main()
{
    gl_Position = projection * vec4(aPos, 1.0);
    TexCoord = aTexCoord;
}
)";
const char* fragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;
in vec2 TexCoord;
uniform sampler2D ourTexture;
void main()
{
    FragColor = texture(ourTexture, TexCoord);
}
)";

const char* vertexOvlSource = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
uniform mat4 projection;
void main()
{
    gl_Position = projection * vec4(aPos, 1.0);
}
)";
const char* fragmentOvlSource = R"(
#version 330 core
out vec4 FragColor;
uniform vec3 objectColor;
void main()
{
    FragColor = vec4(objectColor, 1.0);
}
)";

int texture_width_ = 0;
int texture_height_ = 0;
unsigned int texture_id_;

std::unique_ptr<RenderObject> points_;
std::unique_ptr<RenderObject> lines_x_;
std::unique_ptr<RenderObject> lines_y_;

int key_pressed_ = 0;

const std::string plot_filename_ = "plot.png";

} // end of anonymous namespace

/////////////////////////////////////////////////////////////////////////

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void key_callback(GLFWwindow* window, int key, int scancode, int action,
                  int mods);
void mouse_button_callback(GLFWwindow* window, int button, int action,
                           int mods);
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
void reloadTexture(GLFWwindow* window, const std::string& filename);


/////////////////////////////////////////////////////////////////////////
// Example of creating plot from given function

void on_key_a_pressed(GLFWwindow* window)
{
    const double PI = 3.141592653589793;

    plot_data.plot_name = L"sinus";
    plot_data.line_type = L"solid";
    plot_data.rgb[0] = 0.0;
    plot_data.rgb[1] = 0.0;
    plot_data.rgb[2] = 1.0;

    GeneratePlotFromFunc(
        "plot.png", [](double x) { return sinf((float)x); }, 16, -PI, PI);
    reloadTexture(window, "plot.png");
}

/////////////////////////////////////////////////////////////////////////
// Example of creating hardcoded plot

void on_key_s_pressed(GLFWwindow* window)
{
    std::vector<double> xs = { -2, -1, 0, 1, 2 };
    std::vector<double> ys = { 2, -1, -2, -1, 2 };

    plot_data.plot_name = L"x^2-2";
    plot_data.line_type = L"dotted";
    plot_data.rgb[0] = 0.0;
    plot_data.rgb[1] = 1.0;
    plot_data.rgb[2] = 0.0;

    if (!GeneratePlotFromPoints(plot_filename_, xs, ys))
    {
        std::cerr << "Failed to generate initial plot image." << std::endl;
    }
    reloadTexture(window, "plot.png");
}

/////////////////////////////////////////////////////////////////////////
// Example of creating plot from click-points

void on_key_d_pressed(GLFWwindow* window)
{
    plot_data.plot_name = L"clicked";
    plot_data.line_type = L"dashed";
    plot_data.rgb[0] = 1.0;
    plot_data.rgb[1] = 0.0;
    plot_data.rgb[2] = 0.0;

    if (!GeneratePlotFromPoints(plot_filename_, plot_data.xs, plot_data.ys))
    {
        std::cerr << "Failed to generate initial plot image." << std::endl;
    }
    reloadTexture(window, "plot.png");

    // clear previous points
    points_->vertices.clear();
    plot_data.xs.clear();
    plot_data.ys.clear();

    updateRenderObject(points_.get());
}

/////////////////////////////////////////////////////////////////////////
// Empty slots for the students

void on_key_f_pressed(GLFWwindow* window) {}

void on_key_g_pressed(GLFWwindow* window) {}

// clear user interaction buffers
void on_key_clear(GLFWwindow* window)
{
    // clear previous points
    points_->vertices.clear();
    plot_data.xs.clear();
    plot_data.ys.clear();
    updateRenderObject(points_.get());

    lines_x_->vertices.clear();
    plot_data.user_range_x[0] = 0.0;
    plot_data.user_range_x[1] = 0.0;
    updateRenderObject(lines_x_.get());

    lines_y_->vertices.clear();
    plot_data.user_range_y[0] = 0.0;
    plot_data.user_range_y[1] = 0.0;
    updateRenderObject(lines_y_.get());
}

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Reloads the PNG from disk and updates the active OpenGL texture.
 * @param filename The path to the PNG image file.
 */
void reloadTexture(GLFWwindow* window, const std::string& filename)
{
    int new_width, new_height, nrChannels;
    // Force loading the image with 4 channels (RGBA) for consistency
    unsigned char* data =
        stbi_load(filename.c_str(), &new_width, &new_height, &nrChannels, 4);

    if (data)
    {
        texture_width_ = new_width;
        texture_height_ = new_height;

        glBindTexture(GL_TEXTURE_2D, texture_id_);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture_width_, texture_height_,
                     0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glfwSetWindowSize(window, texture_width_, texture_height_);

        stbi_image_free(data);
        std::cout << "Texture '" << filename << "' reloaded (" << texture_width_
                  << "x" << texture_height_ << ")." << std::endl;
    }
    else
    {
        std::cerr << "Failed to reload texture from '" << filename << "'."
                  << std::endl;
    }
}

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Handles all key press events for the application.
 */
void key_callback(GLFWwindow* window, int key, int scancode, int action,
                  int mods)
{
    // We only care about key press events, not releases
    if (action != GLFW_PRESS)
    {
        return;
    }

    key_pressed_ = key;

    switch (key)
    {
        case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, true); break;
        case GLFW_KEY_R: reloadTexture(window, "plot.png"); break;
        case GLFW_KEY_A: on_key_a_pressed(window); break;
        case GLFW_KEY_S: on_key_s_pressed(window); break;
        case GLFW_KEY_D: on_key_d_pressed(window); break;
        case GLFW_KEY_F: on_key_f_pressed(window); break;
        case GLFW_KEY_G: on_key_g_pressed(window); break;

        case GLFW_KEY_C: on_key_clear(window); break;
    }
}

/////////////////////////////////////////////////////////////////////////

void mouse_button_x_range(const double& xpos, const double& ypos)
{
    if (lines_x_->vertices.size() <= 2)
    {
        Vertex v0 = { (float)(xpos / plot_data.pix_x) * 2.f - 1.f, -1.f };
        Vertex v1 = { (float)(xpos / plot_data.pix_x) * 2.f - 1.f, 1.f };

        lines_x_->vertices.push_back(v0);
        lines_x_->vertices.push_back(v1);

        float plot_x = float(xpos - plot_data.pad_x)
            / (plot_data.pix_x - plot_data.pad_x * 2);
        plot_x = plot_data.range_x_min
            + plot_x * (plot_data.range_x_max - plot_data.range_x_min);
        plot_data.user_range_x[0] = plot_x;
    }
    else
    {
        float plot_x = float(xpos - plot_data.pad_x)
            / (plot_data.pix_x - plot_data.pad_x * 2);
        plot_x = plot_data.range_x_min
            + plot_x * (plot_data.range_x_max - plot_data.range_x_min);
        plot_data.user_range_x[1] = plot_x;

        if (plot_data.user_range_x[1] < plot_data.user_range_x[0])
            std::swap(plot_data.user_range_x[1], plot_data.user_range_x[0]);

        key_pressed_ = 0;
    }
    updateRenderObject(lines_x_.get());
}

/////////////////////////////////////////////////////////////////////////

void mouse_button_y_range(const double& xpos, const double& ypos)
{
    if (lines_y_->vertices.size() <= 2)
    {
        Vertex v0 = { -1.f, float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };
        Vertex v1 = { 1.f, float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };

        lines_y_->vertices.push_back(v0);
        lines_y_->vertices.push_back(v1);

        float plot_y = 1.f
            - (float(ypos - plot_data.pad_y)
               / (plot_data.pix_y - plot_data.pad_y * 2));
        plot_y = plot_data.range_y_min
            + plot_y * (plot_data.range_y_max - plot_data.range_y_min);
        plot_data.user_range_y[0] = plot_y;
    }
    else
    {
        float plot_y = 1.f
            - (float(ypos - plot_data.pad_y)
               / (plot_data.pix_y - plot_data.pad_y * 2));
        plot_y = plot_data.range_y_min
            + plot_y * (plot_data.range_y_max - plot_data.range_y_min);
        plot_data.user_range_y[1] = plot_y;

        if (plot_data.user_range_y[1] < plot_data.user_range_y[0])
            std::swap(plot_data.user_range_y[1], plot_data.user_range_y[0]);

        key_pressed_ = 0;
    }
    updateRenderObject(lines_y_.get());
}

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Handles mouse button events, storing left-click coordinates in a
 * buffer.
 */
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        if (key_pressed_ == GLFW_KEY_X)
        {
            mouse_button_x_range(xpos, ypos);
        }
        else if (key_pressed_ == GLFW_KEY_Y)
        {
            mouse_button_y_range(xpos, ypos);
        }
        else
        {
            // update plot resources
            {
                float plot_x = float(xpos - plot_data.pad_x)
                    / (plot_data.pix_x - plot_data.pad_x * 2);
                plot_x = plot_data.range_x_min
                    + plot_x * (plot_data.range_x_max - plot_data.range_x_min);

                float plot_y = 1.f
                    - (float(ypos - plot_data.pad_y)
                       / (plot_data.pix_y - plot_data.pad_y * 2));
                plot_y = plot_data.range_y_min
                    + plot_y * (plot_data.range_y_max - plot_data.range_y_min);

                plot_data.xs.push_back(plot_x);
                plot_data.ys.push_back(plot_y);
            }

            // update rendering resources
            {
                Vertex v = { float(xpos / plot_data.pix_x) * 2.f - 1.f,
                             float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };
                points_->vertices.push_back(v);

                updateRenderObject(points_.get());
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (key_pressed_ == GLFW_KEY_X)
    {
        Vertex v0 = { float(xpos / plot_data.pix_x) * 2.f - 1.f, -1.f };
        Vertex v1 = { float(xpos / plot_data.pix_x) * 2.f - 1.f, 1.f };
        uint8_t i0 = 0, i1 = 1;

        if (lines_x_->vertices.size() <= 2)
        {
            lines_x_->vertices.resize(2);
        }
        else
        {
            i0 = 2;
            i1 = 3;
        }
        lines_x_->vertices[i0] = v0;
        lines_x_->vertices[i1] = v1;
        updateRenderObject(lines_x_.get());
    }
    else if (key_pressed_ == GLFW_KEY_Y)
    {
        Vertex v0 = { -1.f, float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };
        Vertex v1 = { 1.f, float(1.f - ypos / plot_data.pix_y) * 2.f - 1.f };
        uint8_t i0 = 0, i1 = 1;

        if (lines_y_->vertices.size() <= 2)
        {
            lines_y_->vertices.resize(2);
        }
        else
        {
            i0 = 2;
            i1 = 3;
        }
        lines_y_->vertices[i0] = v0;
        lines_y_->vertices[i1] = v1;
        updateRenderObject(lines_y_.get());
    }
}

/////////////////////////////////////////////////////////////////////////
/**
 * @brief Callback for window resize events. Kept for completeness.
 */
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // The viewport is handled dynamically in the main render loop to maintain
    // aspect ratio. This function is still required by GLFW but can be left
    // empty for our purpose.
}

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    ////////////////////////////////
    // zadanie 1
    double matA[4][4] = {
        {1, 0, 2, -1},
        {3, 0, 0, 5},
        {2, 1, 4, -3},
        {1, 0, 5, 0}
    }; // det = 30

    double res = wyznacznik((double*)matA, 4);

    ////////////////////////////////
    // zadanie 2
    {
#if 0
        constexpr int N=3;
        double mat[N][N] = {
            {2, 1, -1},
            {-3, -1, 2},
            {-2, 1, 2}
        };

        double wektor_b[N] = {8, -11, -3};
#else

        constexpr int N=4;
        double mat[4][4] = {
            {2, 1, 0, 0},
            {1, 2, 1, 0},
            {0, 1, 2, 1},
            {0, 0, 1, 2}
        };
        double wektor_b[4] = {1, 2, 3, 4};
#endif

        double x[N];

        gaussEliminationWithPivoting((double*)mat, (double*)wektor_b, N);
        backSubstitution((double*)mat, (double*)wektor_b, x, N);

        printf("Rozwi¹zanie:\n");
        for (int i=0; i<N; i++)
            printf("x[%d] = %f\n", i, x[i]);

        double det=mat[0][0];
        for(int i=1; i<N; i++)
            det*=mat[i][i];

        printf("wyznacznik = %f\n", det);

    }


    ////////////////////////////////////////////////////////
    // zadanie 3
    {
        constexpr int N=4;
        double A[N][N] = {
            {10, -1, 2, 0},
            {-1, 11, -1, 3},
            {2, -1, 10, -1},
            {0, 3, -1, 8}
        };
        double b[4] = {6, 25, -11, 15};  // Rozwi¹zanie: [1, 2, -1, 1]

        double x[4] = {0, 0, 0, 0};

        jacobi((double*)A, b, x, N);

        // Wyœwietlenie wyników
        printf("Rozwi¹zanie:\n");
        for (int i = 0; i < N; i++) {
            printf("x[%d] = %f\n", i, x[i]);
        }
    }


    // --- initialize from command line ---
    plot_data.pad_x = plot_data.pix_x / 20;
    plot_data.pad_y = plot_data.pix_y / 20;

    // --- parse command line options ---
    CLI::App app{ "Numerical Recipes Plotter" };

    app.add_option("-x, --plot-width", plot_data.pix_x, "Plot width in pixels")
        ->default_val(640);
    app.add_option("-y, --plot-height", plot_data.pix_y,
        "Plot height in pixels")
        ->default_val(480);

    app.add_option("-a,--pad-width", plot_data.pad_x,
        "Plot padding width in pixels")
        ->default_val(32);
    app.add_option("-s,--pad-height", plot_data.pad_y,
        "Plot padding height in pixels")
        ->default_val(24);

    app.add_option("-z,--range-minx", plot_data.range_x_min, "Plot min X")
        ->default_val(-10.f);
    app.add_option("-c,--range-maxx", plot_data.range_x_max, "Plot max X")
        ->default_val(10.f);

    app.add_option("-t,--range-miny", plot_data.range_y_min, "Plot min Y")
        ->default_val(-10.f);
    app.add_option("-u,--range-maxy", plot_data.range_y_max, "Plot max Y")
        ->default_val(10.f);

    // 2. Parse the command line.
    // This macro includes a try/catch block and will exit cleanly on --help.
    CLI11_PARSE(app, argc, argv);


    if (!GenerateEmptyPlot(plot_filename_))
    {
        std::cerr << "Failed to generate initial plot image." << std::endl;
        return -1;
    }

    // --- GLFW initialization ---
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Load initial image dimensions to create a well-sized window
    int initial_width, initial_height, nrChannels;
    unsigned char* initial_data =
        stbi_load(plot_filename_.c_str(), &initial_width, &initial_height,
            &nrChannels, 0);
    if (!initial_data)
    {
        std::cerr << "Failed to load initial texture to determine size."
            << std::endl;
        glfwTerminate();
        return -1;
    }
    stbi_image_free(initial_data);

    // --- GLFW window creation ---
    GLFWwindow* window =
        glfwCreateWindow(initial_width, initial_height,
            "PNG Viewer | Press 'R' to refresh", NULL, NULL);
    if (window == NULL)
    {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // --- Set Callbacks ---
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);

    // --- GLAD ---
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    {
        float pcol[] = { 1.f, 0.f, 0.f };
        points_.reset(createRenderObject(GL_POINTS, pcol, 4.f));
    }

    {
        float lcol[] = { 0.f, 1.f, 0.f };
        lines_x_.reset(createRenderObject(GL_LINES, lcol, 2.f));
    }

    {
        float lcol[] = { 0.f, 0.f, 1.f };
        lines_y_.reset(createRenderObject(GL_LINES, lcol, 2.f));
    }

    // --- Shaders (compile and link) ---
    unsigned int shaderProgram = 0, overlayProgram = 0;
    {
        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glCompileShader(vertexShader);
        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(fragmentShader);
        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
    }

    {
        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexOvlSource, NULL);
        glCompileShader(vertexShader);
        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentOvlSource, NULL);
        glCompileShader(fragmentShader);
        overlayProgram = glCreateProgram();
        glAttachShader(overlayProgram, vertexShader);
        glAttachShader(overlayProgram, fragmentShader);
        glLinkProgram(overlayProgram);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
    }


    // --- Vertex Data for a Plot Quad ---
    float vertices[] = { // positions      // texture coords
        1.0f, 1.0f,  0.0f, 1.0f,  0.0f,  1.0f, -1.0f,
        0.0f, 1.0f,  1.0f, -1.0f, -1.0f, 0.0f, 0.0f,
        1.0f, -1.0f, 1.0f, 0.0f,  0.0f,  0.0f
    };
    unsigned int indices[] = { 0, 1, 3, 1, 2, 3 };
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices,
        GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
        (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
        (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // --- Texture Generation ---
    glGenTextures(1, &texture_id_); // Use the global texture ID
    glBindTexture(GL_TEXTURE_2D, texture_id_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Load the initial texture data using our new function
    reloadTexture(window, plot_filename_);

    // --- Render loop ---
    while (!glfwWindowShouldClose(window))
    {
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // --- Dynamic Viewport for Aspect Ratio Correction ---
        int win_width, win_height;
        glfwGetFramebufferSize(window, &win_width, &win_height);
        float texture_aspect = (texture_width_ > 0 && texture_height_ > 0)
            ? (float)texture_width_ / (float)texture_height_
            : 1.0f;
        int view_width = win_width;
        int view_height = (int)(view_width / texture_aspect);
        if (view_height > win_height)
        {
            view_height = win_height;
            view_width = (int)(view_height * texture_aspect);
        }
        int view_x = (win_width - view_width) / 2;
        int view_y = (win_height - view_height) / 2;
        glViewport(view_x, view_y, view_width, view_height);

        // --- Draw Quad ---
        glUseProgram(shaderProgram);
        float projection[4][4] = { { 1.0f, 0.0f, 0.0f, 0.0f },
            { 0.0f, 1.0f, 0.0f, 0.0f },
            { 0.0f, 0.0f, 1.0f, 0.0f },
            { 0.0f, 0.0f, 0.0f, 1.0f } };
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1,
            GL_FALSE, &projection[0][0]);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture_id_);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // --- Draw Overlays ---
        glUseProgram(overlayProgram);
        glUniformMatrix4fv(glGetUniformLocation(overlayProgram, "projection"),
            1, GL_FALSE, &projection[0][0]);
        drawRenderObject(points_.get(), overlayProgram);
        drawRenderObject(lines_x_.get(), overlayProgram);
        drawRenderObject(lines_y_.get(), overlayProgram);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // --- Cleanup ---
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glDeleteProgram(shaderProgram);
    glDeleteTextures(1, &texture_id_);

    glfwTerminate();
    return 0;
}
