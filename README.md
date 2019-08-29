# Roughness analysis of 3D meshes

This tool was written by Léa Paulus as part of her master thesis supervised by [Prof. Greet Kerckhofs](https://uclouvain.be/fr/repertoires/greet.kerckhofs).
This started as a fork of [the code provided by Guillaume Lavoué](https://perso.liris.cnrs.fr/guillaume.lavoue/rech/soft.html) but I fixed several bugs and added new features as highlighted in my [thesis report](https://dial.uclouvain.be/downloader/downloader_thesis.php?pid=thesis:22284&datastream=PDF_01&key=2eeb063b64f1a3defa14fd0b7d6151aa).

To use the tool, first install it as follows:
```
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```
Note that CGAL and BOOST need to be installed for the installation to work.

It should produce the two executables `build/Roughness` and `build/Color`.
There are scripts available to automate the use of these programs provided that the meshes are stored in `data` and the output is produced in `results`.
For instance, if your data is available at `/media/user/MyExternalDrive/` you can do a symbolic link from `data` to `/media/username/MyExternalDrive/`
to avoid having to copy it to your disk while still having it available in the `data` folder:
```
$ ln -s /media/username/MyExternalDrive data
```
Then given a mesh `meshname.obj` in `data`, you can compute its roughness using
```
$ bash process.sh meshname 1e-2 5e-3 2e-2
```
The results will be stored in `results/meshname/1e-2/5e-3/2e-2`.
The roughness and curvature maps are exported in [`.off` files](https://en.wikipedia.org/wiki/OFF_(file_format)) and can for instance be visualized with [MeshLab](http://www.meshlab.net/).

Note that the smooth mesh is stored in `results/meshname/1e-2/smooth.off` so with the following command:
```
$ bash process.sh meshname 1e-2 6e-3 2e-2
```
as the `SmoothRadius` is still `1e-2`, the smoothing is not computed again and the smooth mesh is simply read from `results/meshname/1e-2/smooth.off`.
Moreover, the curvature is stored in `results/meshname/1e-2/5e-3/smooth_kmax.txt` and `results/meshname/1e-2/5e-3/original_kmax.txt` so with the following command:
```
$ bash process.sh meshname 1e-2 5e-3 1e-2
```
the curvature will not be computed again and will simply be read from the files just mentioned.
This allows a significant time saving when trying the tool for the same mesh with different radii.

To analyse the roughness parameters of the obtained results, some Octave/Matlab scripts are available as well:
`roughness_stats.m` plots statistics of the roughness parameters with confidence intervals and
`roughness_table.m` prints LaTeX tables with the values of all roughness parameters.
