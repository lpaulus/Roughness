#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <limits>

// See https://www.kennethmoreland.com/color-advice/
/*
#define NCOLORS 3
double colormap[3][3] = {
    { 85.0/255,  72.0/255, 193.0/255},
    {221.0/255, 221.0/255, 221.0/255},
    {177.0/255,   1.0/255,  39.0/255}
};
*/
#define NCOLORS 9
double colormap[9][3] = {
    {  0.0/255,   0.0/255,   0.0/255},
    { 44.0/255,   5.0/255, 103.0/255},
    {  3.0/255,  67.0/255,  67.0/255},
    {  5.0/255, 103.0/255,  13.0/255},
    {117.0/255, 124.0/255,   6.0/255},
    {246.0/255, 104.0/255,  74.0/255},
    {250.0/255, 149.0/255, 241.0/255},
    {232.0/255, 212.0/255, 253.0/255},
    {255.0/255, 255.0/255, 255.0/255}
};

void color_map (double *color, double x, double min, double max, bool apply_sqrt) {
    //printf("%lf %lf %lf\n", min, x, max);
    if (x < min) {
        memcpy(color, colormap[0], 3);
    } else if (x > max) {
        memcpy(color, colormap[NCOLORS-1], 3);
    } else {
        double scaling = (x - min) / (max - min);
        if (apply_sqrt) {
            //scaling = sqrt(scaling);
            scaling = pow(scaling, 0.8);
        }
        int i = 0;
        for (i = 0; i < NCOLORS - 2 && scaling > (i + 1.0) / (NCOLORS - 1.0); i++);
        //printf("%lf %lf %d\n", scaling, (i + 1.0) / (NCOLORS - 1.0), i);
        scaling *= NCOLORS - 1;
        scaling -= i;
        for (int j = 0; j < 3; j++) {
            color[j] = colormap[i][j] + scaling * (colormap[i+1][j] - colormap[i][j]);
        }
    }
}

#define LINE_LENGTH 256
char line[LINE_LENGTH];

// Similar to `fgets(line, LINE_LENGTH, f)` when `eol=true`.
void read(FILE *f, bool eol) {
    int i = 0;
    for (i = 0; i < LINE_LENGTH - 1; i++) {
	char c = fgetc(f);
	if (c == '\n' || c == EOF) {
	    break;
        }
        line[i] = c;
    }
    if (i == 0) {
        printf("Empty line\n");
	exit(EXIT_FAILURE);
    }
    if (i == LINE_LENGTH - eol) {
        printf("Line too long, increase `LINE_LENGTH`\n");
	exit(EXIT_FAILURE);
    }
    if (eol) {
        line[i] = '\n';
        i++;
    }
    line[i] = '\0';
}

bool endsWith(const char *str1, const char *str2)
{
    int l1 = strlen(str1);
    int l2 = strlen(str2);
    return l1 >= l2 && !strcmp(str1+l1-l2, str2);
}

int main(int argc, char *argv[])
{
    if (argc <= 2 || argc > 5) {
        printf("Color model.off attribute.txt [output.off] [-l]\n");
	return EXIT_FAILURE;
    }
    if (!endsWith(argv[1], ".off")) {
        printf("First argument should have .off extension, got %s\n", argv[1]);
	return EXIT_FAILURE;
    }
    char *input_mesh_filename = argv[1];
    if (!endsWith(argv[2], ".txt")) {
        printf("First argument should have .txt extension, got %s\n", argv[2]);
	return EXIT_FAILURE;
    }
    char *input_attr_filename = argv[2];
    char *output_mesh_filename = new char[16];
    strcpy(output_mesh_filename, "output.off");
    bool apply_sqrt = false;
    if (argc > 3) {
        if (!strcmp(argv[3], "-s")) {
            apply_sqrt = true;
        } else if (endsWith(argv[3], ".off")) {
    	delete[] output_mesh_filename;
            output_mesh_filename = argv[3];
        } else {
            printf("Invalid argument %s, should be -s or output.off\n", argv[3]);
            return EXIT_FAILURE;
        }
    }
    if (argc > 4) {
        if (!strcmp(argv[4], "-s")) {
            apply_sqrt = true;
        }
    }
    // See https://people.sc.fsu.edu/~jburkardt/data/off/off.html
    FILE *input_mesh = fopen(input_mesh_filename, "r");
    if (input_mesh == NULL) {
        fprintf(stderr, "Invalid file name %s\n", input_mesh_filename);
        return EXIT_FAILURE;
    }
    FILE *input_attr = fopen(input_attr_filename, "r");
    if (input_attr == NULL) {
        fprintf(stderr, "Invalid file name %s\n", input_attr_filename);
        return EXIT_FAILURE;
    }
    FILE *output_mesh = fopen(output_mesh_filename, "w");
    if (output_mesh == NULL) {
        fprintf(stderr, "Invalid file name %s\n", output_mesh_filename);
        return EXIT_FAILURE;
    }
    read(input_mesh, false);
    if (strcmp(line, "OFF") != 0) {
        fprintf(stderr, "Invalid %s: Does not start with OFF but starts with: %s\n", input_mesh_filename, line);
        return EXIT_FAILURE;
    }
    fputs("COFF\n", output_mesh);
    read(input_mesh, true);
    while (line[0] == '#') {
        read(input_mesh, true);
    }
    int vertex_count = 0, face_count = 0, edge_count = 0;
    int num = sscanf(line, "%d %d %d", &vertex_count, &face_count, &edge_count);
    if (num != 3) {
        fprintf(stderr, "Invalid %s: Does not provide 3 integers but provides: %s\n", input_mesh_filename, line);
        return EXIT_FAILURE;
    }
    fprintf(output_mesh, "%d %d %d\n", vertex_count, face_count, edge_count);
    double attr = 0.0;
    double min_attr = std::numeric_limits<double>::max();
    double max_attr = std::numeric_limits<double>::min();
    for (int i = 1; i <= vertex_count; i++) {
        read(input_attr, false);
        num = sscanf(line, "%lf", &attr);
        if (num != 1) {
	    fprintf(stderr, "Invalid %s\n", input_attr_filename);
	    return EXIT_FAILURE;
        }
        if (attr < min_attr)
            min_attr = attr;
        if (attr > max_attr)
            max_attr = attr;
    }
    printf("min = %lf\n", min_attr);
    printf("max = %lf\n", max_attr);
    if (apply_sqrt) {
        printf("Applying sqrt\n");
    } else {
        printf("Not applying sqrt\n");
    }
    rewind(input_attr);
    double color[3];
    for (int i = 1; i <= vertex_count; i++) {
	read(input_mesh, false);
        fputs(line, output_mesh);
        read(input_attr, false);
        num == sscanf(line, "%lf", &attr);
        if (num != 1) {
	    fprintf(stderr, "Invalid %s\n", input_attr_filename);
	    return EXIT_FAILURE;
        }
        color_map(color, attr, min_attr, max_attr, apply_sqrt);
        fprintf(output_mesh, " %lf %lf %lf 1.0\n",  color[0], color[1], color[2]);
    }
    for (int i = 1; i <= face_count; i++) {
	read(input_mesh, true);
        fputs(line, output_mesh);
    }
    fclose(input_mesh);
    fclose(input_attr);
    fclose(output_mesh);

    return EXIT_SUCCESS;
}
