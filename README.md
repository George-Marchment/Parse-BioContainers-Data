# Parse-BioContainers-Data

This repository contains : 
* A notebook which downloads the [BioContainers](https://biocontainers.pro/) data from this adress : https://github.com/BioContainers/containers. Parsers the files, and extarcts for each tool its corresponding tags and summaries. It then saves the structured data into a *.json* file, then deletes the original downloaded BioContainers data.
* It also containes multiple notebooks which performs some small studies on the BioContainers data:
  * The first does a little study on the capacity of the annotations to cluster similar tools together (see [here](studies/clustering_study.ipynb))
  * The second performs a diversity study, on a [corpus](https://github.com/George-Marchment/Process-Comparison-Dataset/tree/main/data) of 30 (rna-seq) workflows (see [here](studies/diversity_study.ipynb))