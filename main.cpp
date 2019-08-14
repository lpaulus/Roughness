#include <iostream>
#include "Roughness.h"

using namespace std;

bool endsWith(const char *str1, const char *str2)
{
    int l1 = strlen(str1);
    int l2 = strlen(str2);
    return l1 >= l2 && !strcmp(str1+l1-l2, str2);
}

int main(int argc, char *argv[])
{
    if(argc <= 2 || (!endsWith(argv[1], ".off") && !endsWith(argv[1], ".obj")))
        printf("Roughness model.[off|obj] SmoothRadius [CurvatureRadius] [AverageRadius]\n");
    else
    {
        PolyhedronPtr poly = PolyhedronPtr(new Polyhedron());
        printf("\e[1;1m\e[38;5;087m┌ Loading\n");
        if (endsWith(argv[1], ".obj")) {
            poly->load_mesh_obj(argv[1]);
        } else {
            poly->load_mesh(argv[1]);
        }
        printf("\e[1;1m\e[38;5;087m└\e[0m\n");

        printf("\e[1;1m\e[38;5;087m┌ Normalise mesh to have unit bounding box\n");
        // The mesh is recentered and every vertices position are divided by `scaling`.
        // We will need to adapt the roughness accordingly.
        double scaling = poly->Normalise();
        poly->write_off("normalised.off", false, false);
        printf("\e[1;1m\e[38;5;087m└\e[0m written to normalised.off\n");
        poly->compute_bounding_box();
        poly->compute_normals();
        poly->compute_type();
        poly->calc_nb_components();
        poly->calc_nb_boundaries();

        double SmoothRadius;
        int err = sscanf(argv[2], "%lf", &SmoothRadius);
        if (err != 1) {
            fprintf(stderr, "[ Error: Invalid SmoothRadius: %s\n", argv[2]);
        }

        double CurvatureRadius = 0.005;
        if(argc > 3) {
            err = sscanf(argv[3], "%lf", &CurvatureRadius);
            if (err != 1) {
                fprintf(stderr, "[ Error: Invalid CurvatureRadius: %s\n", argv[3]);
            }
        }

        double AverageRadius = 2 * SmoothRadius;
        if(argc > 4) {
            err = sscanf(argv[4], "%lf", &AverageRadius);
            if (err != 1) {
                fprintf(stderr, "[ Error: Invalid AverageRadius: %s\n", argv[4]);
            }
        }

        CRoughness<Polyhedron> roughness(poly.get());
        roughness.compute_Roughness(AverageRadius, SmoothRadius, CurvatureRadius, scaling);
    }

    return EXIT_SUCCESS;
}
