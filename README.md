## Nonlinear Beamforming
A Nonlinear Beamforming for Enhanced Spatiotemporal Sensitivity in High Frame Rate Ultrasound Flow Imaging

Repository to share scripts and functions for Nonlinear Beamforming for Ultrasound Flow Imaging. All functions are usable with agreement from their owner. 

###### DATE 05-05-2022 : VERSION 1.0


#### AUTHORS: 
A. N. Madhavanunni and Mahesh Raveendranatha Panicker <br />
Center for Computational Imaging, Indian Institute of Technology Palakkad, Kerala, India.

#### Instructions for the execution of codes 

Please follow the below instructions to run the code.
1. Download the source_code and data folder from the repository
2. Please note that the source code has the following dependencies:<br />
      (a) MUST Toolbox (ver.2.0): used for robust smoothing and vector flow visualisation <br />
      (b) Time-Frequency Toolbox (tftb-ver.0.2): used to obtain localised time-frequency plots of the beamformed signal <br />
      (c) The _show_Figure6.m_ uses _tiledlayout_ and requires MATLAB 2019b or higher for its execution
3. Download the [MUST toolbox](https://www.biomecardio.com/MUST/index.html),  [Time-Frequency Toolbox](https://tftb.nongnu.org/) and place it in _source_code/lib/_  
4. Set appropriate values for _bfParams.beamApod_ and _bfParams.DMAS_ in _mainCode_expPWIdisc.m_ and run the same. Please ensure that the paths provided in _mainCode_expPWIdisc.m_ is valid.
5. The _show_Figure3.m_ should be executed only after executing _mainCode_expPWIdisc.m_ with _saveEnable_ set to 1.


##### Code Available under Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) (see https://creativecommons.org/licenses/by-nc-nd/4.0/)

##### Please see the Academic references and Acknowledgements that are to be cited for any usage of the code and/or data.

#### ACADEMIC REFERENCES TO BE CITED:
1. Details of nonlinear beamforming is available in the article by A. N. Madhavanunni and Mahesh Raveendranatha Panicker <br />
[A Nonlinear Beamforming for Enhanced Spatiotemporal Sensitivity in High Frame Rate Ultrasound Flow Imaging](https://doi.org/10.48550/arXiv.2108.02688)<br />
Video results are available at https://www.youtube.com/playlist?list=PLiuuVhVNWBZSYikqhd20FsVr8NTKRlZ4F

2. More details of dual apodization applied to vector flow imaging is available in the following article:<br />
A. N. Madhavanunni and Mahesh Raveendranatha Panicker, [Directional Beam Focusing Based Dual Apodization Approach for Improved Vector Flow Imaging](https://doi.org/10.1109/ISBI45749.2020.9098494) _IEEE 17th International Symposium on Biomedical Imaging (ISBI)_, Iowa, USA, 2020. 

3. Description of the triangulation based single transmit dual angle multi-receive scheme for velocity estimation is available in the following article:<br />
A. N. Madhavanunni and Mahesh Raveendranatha Panicker, [Triangulation based vector flow imaging with non-steered plane waves for transverse flows](https://doi.org/10.1117/12.2549253), Medical Imaging 2020: Ultrasonic Imaging and Tomography. Vol. 11319. International Society for Optics and Photonics, Houston, USA, 2020

#### Acknowledgements:
##### [MUST Toolbox](https://www.biomecardio.com/MUST/index.html) for the _in-vitro_ rotating disk dataset, robust smoothing (smoothn function) and vector flow visualisation (vplot function): 
1. Damien Garcia. [Make the most of MUST, an open-source Matlab UltraSound Toolbox](https://doi.org/10.1109/IUS52206.2021.9593605). IEEE International Ultrasonics Symposium, IUS, 2021.
2. Craig Madiena, Julia Faurie, Jonathan Poree, and Damien Garcia., [Color and Vector Flow Imaging in Parallel Ultrasound with Sub-Nyquist Sampling](https://doi.org/10.1109/TUFFC.2018.2817885). IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 65(5):795–802, 2018.
3. Garcia D. [Robust smoothing of gridded data in one and higher dimensions with missing values](https://www.biomecardio.com/publis/csda10.pdf). Computational Statistics & Data Analysis 2010; 54:1167-1178.
4. Garcia D. [A fast all-in-one method for automated post-processing of PIV data](https://www.biomecardio.com/publis/expfluids11.pdf). Exp Fluids 2011; 50:1247-1259.


##### [Time-Frequency Toolbox](https://tftb.nongnu.org/) for the localized time-frequency plots using reassigned pseudo Wigner-Ville distribution (RPWVD) method (tfrrpwv function in the toolbox)
1. Franqois Auger and Patrick Flandrin. [Improving the Readability of Time-Frequency and Time-Scale Representations by the Reassignment Method](https://doi.org/10.1109/78.382394). IEEE Transactions on Signal Processing, 43(5), 1995.
