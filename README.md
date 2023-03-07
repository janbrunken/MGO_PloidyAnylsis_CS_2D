# PloidyAnalysis_2D
Analyze ploidy from 2D fluorescecne microscopy image data. 
This repository contains the PloidyAnalysis_2D python package as well as Juypter Notebooks to perform ploidy analysis on 2D immunofluorescence images. Image segmentation is performed using pre-trained models of [StarDist](https://github.com/stardist/stardist/) (i.e. '2D_versatile_fluo') and [Cellpose 2.0](https://github.com/mouseland/cellpose) ('cyto'). The whole analysis pipeline consists of 4 major steps:

  1. Nuclei segmentation
  2. Marker segmentation
  3. Image preprocessing (matching labels, removing unmatched labels, reindexing labels)
  4. Ploidy analysis (filtering out cycling cells, ploidy determination, cell type filtering)

## Image requirements
Cells to be analyzed should be stained against the following markers:
  - Nuclear labeling reflecting DNA content (e.g. DAPI, HOECHST)
  - Cell cycle marker (e.g. Geminin)
  - Cell type marker of interest (e.g. Cytokeratin 8 for luminal mammary gland cells)
To perform 2D segmentation, maximum intensity projections are recommended. Ploidy measurements should be performed on sum intensity projections.

