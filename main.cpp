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
        printf("Roughness model.[off|obj] epsilon [CurvatureRadius] [roughness.txt]\n");
    else
    {
        PolyhedronPtr poly = PolyhedronPtr(new Polyhedron());
        printf("Loading\n");
        if (endsWith(argv[1], ".obj")) {
            poly->load_mesh_obj(argv[1]);
        } else {
            poly->load_mesh(argv[1]);
        }

        printf("Normalise\n");
        poly->Normalise();
        poly->write_off("normalised.off", false, false);
        printf("Bounding box\n");
        poly->compute_bounding_box();
        printf("Normals\n");
        poly->compute_normals();
        printf("Type\n");
        poly->compute_type();
        printf("# components\n");
        poly->calc_nb_components();
        printf("# boundaries\n");
        poly->calc_nb_boundaries();

        double epsilon;
        sscanf(argv[2], "%lf", &epsilon);

        double CurvatureRadius = epsilon / 2;
        std::string outputFilename = "roughness.txt";
        if(argc > 3) {
            int err = sscanf(argv[3], "%lf", &CurvatureRadius);
            if (err != 1) {
                // It wasn't a double, might be the filename
                outputFilename = argv[3];
            }
        }
        if(argc > 4)
            outputFilename = argv[3];

        CRoughness<Polyhedron> roughness(poly.get());
        printf("Compute roughness\n");
        roughness.compute_Roughness(2*epsilon, epsilon, CurvatureRadius);

        printf("Write roughness\n");
        FILE *fichier = fopen(outputFilename.c_str(), "w");
        for (Vertex_iterator pVertex = poly->vertices_begin(); pVertex != poly->vertices_end(); pVertex++)
        {
            fprintf(fichier, "%lf\n", pVertex->Roughness());
        }
        fclose(fichier);
    }

    return EXIT_SUCCESS;
}
