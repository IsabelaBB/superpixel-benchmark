# Evaluation measures

This code contais the following superpixel evalution measures: Similarity between Image and Reconstruction from Superpixels (SIRS)[1], Boundary Recall (BR)[2], Explained Variation (EV)[3], Undersegmentation Error (UE)[4], and Compactness (CO)[5]. 

### Requirements
The project was developed in **C/C++** under a **Linux-based** operational system; therefore, it is **NOT GUARANTEED** to work properly in other systems (_e.g._ Windows and macOS). It's also required to install OpenCV 4.
        
### Compiling and cleaning
- To compile all files: `make`
- For removing all generated files from source: `make clean`

### Running
Usage: `./bin/main [OPTIONS]`

Options:
```
--img 		:	Original image (eval 1,2, and 7) or ground-truth (eval 3 and 4) file/path
--eval          :       Superpixel evaluation option. Type: int. {1:SIRS, 2:EV, 3:BR, 4:UE, 5:CO, 6:Enforce connectivity, 7:Enforce superpixels' number}
--label 	: 	Segmented image file/path (pgm/png images)
--ext 		: 	Extension of segmented image (defaut: pgm)

--buckets 	: 	Number of color subsets in SIRS evaluation (eval 1) (default:16)
--alpha 	: 	Number of subsets used to represent a superpixel in SIRS evaluation (eval 1) (default:4)
--gaussVar 	: 	Variance of gaussian in SIRS evaluation (eval 1) (default: 0.01)
--k             :       Desired number of superpixels. Used in eval 7. Type: int
--imgScores 	: 	File/Path of the colored result of color homogeneity in SIRS/EV evaluation (eval 1 or 2) (optional)
--drawScores 	: 	Boolean option {0,1} to write scores in the colored image result (imgScores option). Used in SIRS/EV evaluation (eval 1 or 2) (optional)
--log   	: 	txt log file with the mean evaluation results of a measure for a directory (optional)
--dlog 		: 	txt log file with the evaluation results of a measure for all images (optional)
--recon 	: 	File/Path of image reconstruction. Can be used in SIRS/EV (eval 1 or 2) (optional)
--save          :       Save image superpixels after enforce coonectivity/minimum number of superpixels. Can be used when enforce connectivity or enforce superpixels' number (eval 6 or 7) (optional)
```

**Examples:**
- Simple example: `./bin/main --img ./image.jpg --label ./label_500.pgm --imgScores ./result.png`
- Example with image scores: `./bin/main --img ./image.jpg --label ./label_100.pgm --imgScores ./result.png --drawScores 1`

## Cite
If this work was useful for your research, please cite our paper:
```
@InProceedings{barcelos2023review,
  title={A comprehensive review and new taxonomy on superpixel segmentation},
  author={Barcelos, Isabela Borlido and Bel{\'e}m, Felipe and Melo, Leonardo, JR, Zenilton K. G. do Patroc{\'i}nio, and Falc{\~a}o, Alexandre Xavier and Guimar{\~a}es, Silvio Jamil F},
  booktitle={ACM Computing Surveys},
  pages={},
  year={2023},
  organization={ACM}
}
```
## References
[1] Isabela B Barcelos, Felipe De C Belém, Leonardo De M João, Alexandre X Falcão, and Guimarães Silvio JF. 2022.
Improving color homogeneity measure in superpixel segmentation assessment. In 2022 35th SIBGRAPI Conference on
Graphics, Patterns and Images (SIBGRAPI), Vol. 1. 79–84. https://doi.org/10.1109/SIBGRAPI55357.2022.9991772.

[2] David R Martin, Charless C Fowlkes, and Jitendra Malik. 2004. Learning to detect natural image boundaries using
local brightness, color, and texture cues. IEEE transactions on pattern analysis and machine intelligence 26, 5 (2004),
530–549. https://doi.org/10.1109/TPAMI.2004.1273918

[3] Alastair P Moore, Simon JD Prince, Jonathan Warrell, Umar Mohammed, and Graham Jones. 2008. Superpixel lattices.
In 2008 IEEE conference on computer vision and pattern recognition. IEEE, 1–8.

[4] Peer Neubert and Peter Protzel. 2012. Superpixel benchmark and comparison. In Proc. Forum Bildverarbeitung, Vol. 6.
KIT Scientific Publishing, 1–12.

[5] Alexander Schick, Mika Fischer, and Rainer Stiefelhagen. 2012. Measuring and evaluating the compactness of
superpixels. In Proceedings of the 21st international conference on pattern recognition (ICPR2012). IEEE, 930–934.

