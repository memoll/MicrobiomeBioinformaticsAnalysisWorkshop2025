<h1 id="title">Bioinformatics Analysis of Microbial Communities</h1>

<p id="GI">
<h2 id="general">General Information</h2>
</p>

<p id="by">
  <strong>Led by:</strong>
  <a href="https://www.arrietalab.com/home">ArrietaLab</a>
</p>

<p id="loc">
  <strong>Location:</strong>
  IUCBC Instituto Universitario de Ciencias Biomédicas de Córdoba - Argentina (Online)
</p>

<p id="date">
  <strong>Dates:</strong>
  July 21-23, 2025
</p>

<p id="obj">
  <strong>Objective:</strong>
  Analyzing microbial DNA sequences
</p>

<p id="requirements">
  <strong>Requirements:</strong> 
</p>

  - Laptop with a Mac, Linux, or Windows operating system 
    (not a tablet, Chromebook, etc.) with administrative privileges 
  - Access to Wifi 
  - R and Rstudio installed (instructions <a href="#setup">below</a>)
  - Excel or any text editor installed (e.g. TextWrangler, Notepad, BBEdit, etc.)
  

<p id="contact">
  <strong>Contact:</strong>
  <a href="mailto:{{mona.parizadeh@ucalgary.ca}}">mona.parizadeh@ucalgary.ca</a> 
</p>

<h2 id="content">Course Content</h2>  
This course covers the following material:

<h3 id="day1">Day 1</h3>
<h4 id="day1h">(9 a.m. to 12 p.m. & 1 p.m. to 4 p.m.)</h4>

 - Introduction to metabarcoding and amplicon sequencing
 - Introduction to R
 - <a href="https://benjjneb.github.io/dada2/tutorial.html">DADA2 Tutorial</a>

<h3 id="day2">Day 2</h3>
<h4 id="day2h">(9 a.m. to 12 p.m. & 1 p.m. to 4 p.m.)</h4>

 - Introduction to phyloseq package in R
 - Data exploration
 - Statistical analyses of 16S rRNA gene sequences:
   - Taxonomic composition
   - Alpha diversity
   - Beta diversity (ordination)
   - PERmutational Multivariate ANalysis Of VAriance (PERMANOVA)
   - Differential analysis

<h3 id="day2">Day 3</h3>
<h4 id="day2h">(9 a.m. to 12 p.m.)</h4>

- Discussion
  - Wrap up of the results of 16S rRNA gene sequence analysis
  - Differences between DADA2 pipeline workflows in analyzing
  <a href="https://benjjneb.github.io/dada2/tutorial.html">16S</a> vs 
  <a href="https://benjjneb.github.io/dada2/ITS_workflow.html">ITS</a> sequences  

     
<h2 id="setup">Setup</h2> 
To participate in this Workshop, please install the following <a href="#r">software</a> 
and <a href="#r-pkg">packages</a>, complete the <a href="#r-course">R course for beginners</a>, and let us know if you need any help before attending.

<div id="r">
  <h3>Install R and RStudio</h3>
  <p>
    <a href="http://www.r-project.org">R</a> is a free and open-source programming 
    language that is particularly powerful for data exploration, visualization, and 
    statistical analysis. We use <a href="https://posit.co/downloads/">RStudio</a> 
    to interact with R.
  </p>
 
 <div class="row">
   <div class="col-md-4">
     <h4 id="r-windows">Windows</h4>
    <p>
     Please download R for Windows
        from <a href="http://cran.r-project.org/index.html">CRAN</a> to install R, and 
        also install the <a href="http://www.rstudio.com/ide/download/desktop">RStudio IDE</a>.
        If you have separate user and admin accounts, please run the installers as an 
        administrator by right-clicking on the .exe file and selecting "Run as administrator" 
        instead of double-clicking. Otherwise, problems may arise later when installing R packages.
    </p>
     <a href="https://www.youtube.com/watch?v=q0PjTAylwoU">Video Tutorial</a>
 </div> 
   
 <div class="col-md-4">
   <p> 
   <h4 id="r-macosx">Mac OS X</h4>
   </p>
   <p>
    Please download R for macOS
       from <a href="http://cran.r-project.org/index.html">CRAN</a> to install R, and also install 
       the <a href="http://www.rstudio.com/ide/download/desktop">RStudio IDE</a>.
   </p>
    <a href="https://www.youtube.com/watch?v=5-ly3kyxwEg">Video Tutorial</a>
  </div> 
   
  <div class="col-md-4">
    <h4 id="r-linux">Linux (Debian, Fedora/Redhat, Ubuntu)</h4>
   <p>
    Please download the binary files for your distribution from
    <a href="http://cran.r-project.org/index.html">CRAN</a> to install R, or use a package manager 
     (e.g. run <code>sudo apt-get install r-base</code> for Debian/Ubuntu and run
        <code>sudo yum install R</code> for Fedora/Redhat). Additionally, please install the
        <a href="http://www.rstudio.com/ide/download/desktop">RStudio IDE</a>.
   </p>
  </div> 
 </div>
</div>
 
 <h3 id="r-course">R for beginners</h3>
To follow the workshop, you must have a basic understanding of R.
Before attending the workshop, please go through the following course:

  - <a href="http://swcarpentry.github.io/r-novice-inflammation/">Programming with R</a>


<h3 id="r-pkg">Install R packages</h3>
Please install the following <a href="https://github.com/memoll/MicrobiomeBioinformaticsAnalysisWorkshop2025/blob/main/package_installation">R packages</a> in RStudio: 


 - dada2
 - phyloseq
 - vegan
 - ggplot2
 - tidyverse
 - dplyr
 - BiocGenerics
 - SummarizedExperiment
 - rstatix
 - DESeq2
 - ggpubr
 - RColorBrewer
 - gridExtra
 
<h3 id="db">Download database for taxonomic assignment with DADA2</h3>

 - <a href="https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1">silva_nr99_v138.2_toGenus_trainset.fa.gz</a>
 - <a href="https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1">silva_nr99_v138.2_toSpecies_trainset.fa.gz</a> 

<h3 id="blast"><a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch">BLAST</a></h3>
<a href="https://github.com/memoll/MicrobiomeBioinformaticsAnalysisWorkshop2025/blob/main/data/sequenceExample">Sequence Example</a>
