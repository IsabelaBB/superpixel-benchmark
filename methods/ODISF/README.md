## Object-based Dynamic and Iterative Spanning Forest (ODISF)

This is the implementation of the superpixel segmentation method _Object-based Dynamic and Iterative Spanning Forest_(ODISF) as proposed in
- **F.Belém, B.Perret, J.Cousty, S.Guimarães, A.Falcão.** [_Towards a Simple and Efficient Object-based Superpixel Delineation Framework_](https://ieeexplore.ieee.org/document/9643123). In 34th International Conference on Graphics, Patterns and Images (SIBGRAPI), pg 346-353. 2021.

This software includes a program for running the ODISF method, and another one for assisting the visualization of the segmentation by overlaying the superpixel borders. Please cite the aforementioned paper if you use any of this software in your own project.

### Hardware, Setup and Requirements

The project was developed in **C** under a **Linux-based** operational system; therefore, it is **NOT GUARANTEED** to work properly in other systems (_e.g._ Windows and macOS). Moreover, the same applies for **non-GCC** compilers, such as Clang and MinGW.

All code within this project were developed, compiled and tested using the following programs:
- **[GNU GCC](https://gcc.gnu.org/)**: version 7.5.0
- **[GNU Make](https://www.gnu.org/software/make/)**: version 4.1

This code was implemented and evaluated in a computer with the following specifications:
- **Model:** Acer X555LB
- **Operational System:** Linux Mint v20.2 x86_64 kernel version 5.4.0-86-generic
- **Order:** Little-Endian
- **CPU:** 4x Dual-core Intel(R) Core(TM) i5-5200 @ 2.20 GHz
- **Memory:** 8GB RAM ; 480 SSD

The library has in-built support for handling **PNM** images. For enabling external library support, please refer to the [README](externals/README.md) file within the **externals** folder.

### Compiling

If your computer meets the aforementioned requirements, you may run the commands below for compiling the library and all the demonstration programs.
```bash
make lib
make demo
```
Or simply run one of the following commands for compiling both latter at once.
```bash
make
make all
```

For removing the files generated from compilation, one may run the following rule.
```bash
make clean
```

### Running

After compiling, one may run ODISF for segmenting an image through the following command
```bash
./bin/RunODISF --img path/to/image.ppm --objsm path/to/objsm.pgm --out path/to/segm.pgm
```
Briefly, `--img`,`--objsm`, and `--out` indicate the paths to the image to be segmented and its object saliency map, and to the resulting segmentation, respectively. For other arguments, one may run
```bash
./bin/RunODISF --help
```
for more information.

Given such segmentation, it is possible to draw the superpixel borders over the original image for a better visualization by the following command
```bash
./bin/RunOvlayBorders --img path/to/image.ppm --labels path/to/segm.pgm --out path/to/overlay.ppm
```
In this case, `--img`,`--labels`, and `--out` indicate the paths to the original image and segmentation, and to the resulting overlay image, respectively. Likewise, for other arguments, one may run
```bash
./bin/RunOvlayBorders --help
```
for more information.

### License

All codes within this project are under the **MIT License**. See the [LICENSE](LICENSE) file for more details.

### Acknowledgements

This work was financially supported by the following brazilian research funding agencies: 
- Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq)
- Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES)
- Fundação de Amparo à Pesquisa do Estado de Minas Gerais (FAPEMIG)
- Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP)

### Contact

If you have any questions or faced an unexpected behavior (_e.g._ bugs), please feel free to contact the authors through the following email addresses:
- **Felipe C. Belém**:  [felipe.belem@ic.unicamp.br](mailto:felipe.belem@ic.unicamp.br)
- **Benjamin Perret**:  [benjamin.perret@esiee.fr](mailto:benjamin.perret@esiee.fr)
- **Jean Cousty**:  [jean.cousty@esiee.fr](mailto:jean.cousty@esiee.fr)
- **Silvio Jamil F. Guimarães**:  [sjamil@pucminas.br](mailto:sjamil@pucminas.br)
- **Alexandre X. Falcão**: [afalcao@ic.unicamp.br](mailto:afalcao@ic.unicamp.br)
