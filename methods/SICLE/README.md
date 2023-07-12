## Superpixels through Iterative CLEarcutting (SICLE)

This is the implementation of the superpixel segmentation method _Superpixels through Iterative CLEarcutting_(SICLE) as proposed in:
- **F.Belém, B.Perret, J.Cousty, S.Guimarães, A.Falcão.** [_Efficient Multiscale Object-based Superpixel Framework_](https://doi.org/10.48550/arXiv.2204.03533). ArXiv preprint.
- **F.Belém, I.Borlido, L.João, B.Perret, J.Cousty, S.Guimarães, A.Falcão.** [_Fast and Effective Superpixel Segmentation Using Accurate Saliency Estimation_](https://doi.org/10.1007/978-3-031-19897-7_21). Discrete Geometry and Mathematical Morphology, pg. 261-273, 2022.

This software includes four programs:
- One for running the SICLE method;
- One for assisting the visualization of the segmentation by overlaying the superpixel borders; 
- One for creating an object's minimum bounding box image from a ground-truth;
- One for assessing a superpixel segmentation (Boundary Recall and Under-segmentation Error);

Please cite the aforementioned papers if you use any of this software in your own project.

### Hardware, Setup and Requirements

The project was developed in **C** under a **Linux-based** operational system; therefore, it is **NOT GUARANTEED** to work properly in other systems (_e.g._ Windows and macOS). Moreover, the same applies for **non-GCC** compilers, such as Clang and MinGW.

All code within this project were developed, compiled and tested using the following programs:
- **[GNU GCC](https://gcc.gnu.org/)**: version 7.5.0
- **[GNU Make](https://www.gnu.org/software/make/)**: version 4.1

This code was implemented and evaluated in a computer with the following specifications:
- **CPU:** 64-bit Intel(R) Core(TM) i7-4790 @ 3.20 GHz
- **Memory:** 8GB RAM

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
./bin/RunSICLE --img path/to/image.ppm --objsm path/to/objsm.pgm --out path/to/segm.pgm
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

To create an object saliency map by drawing the object's minimum bounding box, one may run the following command
```bash
./bin/RunBBFromGT --gt path/to/groundtruth.pgm --out path/to/minbbimage.pgm
```
Clearly, `--gt` refers to the object's groundtruth, whereas `--out`, its respective minimum bounding box image.

Finally, you may run the following command for calculating the Boundary Recall and Under-segmentation metrics for a given superpixel segmentation.
```bash
./bin/RunImgMetrics --labels path/to/image.pgm --gt path/to/groundtruth.pgm
```
The `--labels` parameter holds the path to the superpixel segmentation to be evaluated, which will be compared against the ground-truth indicated by `--gt`. For other arguments, one may run
```bash
./bin/RunImgMetrics --help
```

### License

All codes within this project are under the **MIT License**. See the [LICENSE](LICENSE) file for more details.

### Acknowledgements

This work was financially supported by the following funding agencies: 
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
