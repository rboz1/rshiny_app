<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
  <h3 align="center">Gene Expression Analysis App</h3>

  <p align="center">
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Script</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
This Shiny app provides an interactive interface to analyze gene expression data through various modules, including data summary, counts matrix processing, differential expression analysis, and trajectory analysis.

### Features
Summary Tab
- Summary Table: Displays general descriptive statistics (e.g., age, gender distribution).
- Data Table: View and explore the uploaded metadata in a sortable table.
- Plots: Visualize sample-level data using violin plots.

Counts Matrix Tab
- Diagnostic scatter plots for variance and counts.
- Heatmaps and PCA plots of filtered data.

Differential Expression Tab
- Volcano Plot: Customize and visualize differential expression results interactively.
- Differential Expression Results Table: Explore DESeq2 results interactively.

Trajectory Analysis Tab
- Upload Seurat object: Accepts an .rds file containing a processed Seurat object.
- Filtering Options: Filter data based on conditions, clusters, and selected genes.
- Plots: Generate trajectory and pseudotime plots with customizable color palettes.
- Filtered Data Table: View data filtered by the selected criteria.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Built With
- R
- DESeq2
- Seurat
- Monocle3
- Data from: https://www.nature.com/articles/s41586-024-07493-y

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

### Prerequisites
* R

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/rboz1/rshiny_app.git

2. Install packages
   ```
   install.packages(c("tidyverse", "dplyr", "DESeq2", "Seurat", "SeuratObject", 
                   "ggplot2", "gplots", "RColorBrewer", "shiny", 
                   "shinycssloaders", "colourpicker", "DT", "biomaRt"))
   
3. Run the script 
   ```
   shiny::runApp('rshiny_app')
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage
<img width="1439" alt="Screenshot 2024-12-21 at 1 38 46 PM" src="https://github.com/user-attachments/assets/2c7faa01-e7ed-46cf-9034-acaa3da633bc" />
<img width="1439" alt="trajectory" src="https://github.com/user-attachments/assets/29aeda22-6bda-4dbe-9a8d-bcbaaadf1a2a" />

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Rachel - rbozadjian@gmail.com

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: www.linkedin.com/in/rachel-bozadjian-203999109
