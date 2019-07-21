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
    if(argc <= 2 || !endsWith(argv[1], ".obj"))
        printf("Roughness model.obj epsilon (roughness.txt)\n");
    else
    {
        PolyhedronPtr poly = PolyhedronPtr(new Polyhedron());
        poly->load_mesh_obj(argv[1]);

        poly->Normalise();
        poly->compute_bounding_box();
        poly->compute_normals();
        poly->compute_type();
        poly->calc_nb_components();
        poly->calc_nb_boundaries();

        double epsilon;
        sscanf(argv[2], "%lf", &epsilon);

        std::string outputFilename = "roughness.txt";
        if(argc > 3)
            outputFilename = argv[3];

        CRoughness<Polyhedron> roughness(poly.get());
        roughness.compute_Roughness(2*epsilon, epsilon);

        FILE *fichier = fopen(outputFilename.c_str(), "w");
        for (Vertex_iterator pVertex = poly->vertices_begin(); pVertex != poly->vertices_end(); pVertex++)
        {
            fprintf(fichier, "%lf\n", pVertex->Roughness());
        }
        fclose(fichier);
        printf("success\n");
    }

    return 0;
}
