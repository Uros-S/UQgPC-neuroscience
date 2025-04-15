# UQgPC_neuroscience
Code used to generate the results in the conference paper "Efficient gPC-based quantification of probabilistic robustness for systems in neuroscience" by U. Sutulovic, D. Proverbio, R. Katz, and G. Giordano, presented at the European Control Conference (ECC) 2025 in Thessaloniki (Greece).

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".

# Toolbox required

The results of this paper rely on the functions from the toolbox PoCET. Details on installation and usage can be found here:
https://www.tu-chemnitz.de/etit/control/research/PoCET/index.php.en

# Usage

- Download the data from the following Zenodo repository and put the folders in the same working directory of 'main.m' and the folder 'util':
  https://zenodo.org/records/15226040

- In 'main.m' select the model (1 --> Hindmarsh-Rose, 2 --> Jansen-Rit, 3 --> Epileptor) and whether to recompute everything from scratch (recompute = 1) or to use the pre-computed data generated for the paper's figures 
  (recompute = 0). Note that the numerical values will slightly differ from those of the paper since sampling is involved in some parts of the algorithm.
