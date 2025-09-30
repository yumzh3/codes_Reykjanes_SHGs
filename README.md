# Overview
Please read the README.md and LICENSE before using and downloading codes.<br> 
Please follow the CITATION.cff to for citation.<br> 
Python codes for crystallization calculation shown in Figure S3, cross-correlation, and spectral analysis are provided here to replicate results shown in the paper.<br> 
Codes and demo data for crystallization modeling are listed in the folder 'xtal_Yu2025'.<br> 
This is a part of Supplementary Materials for paper "Sensitivity of Reykjanes Ridge volcanism to glacial cycles over the last 130,000 years" written by Mingzhen Yu (corresponding author, myu@g.harvard.edu) et al. being submitted to Science Advances. Codes are corresponding to Mingzhen Yu (Department of Earth and Planetary Sciences, Harvard University, Cambridge, MA 02138, USA).<br> 
Codes are written in Python.<br> 

# File Introduction
## SHG_cross_correlation.py
This python code can be used to replicate the cross-correlation results shown in the paper with data loaded.<br>
The time series data of SHGs are provided in the Supplementary DataS3.<br> 
## SHG_spectral_analysis.py
This python code can be used to replicate the spectral analysis results shown in the paper with data loaded.<br>
The time series data of SHGs are provided in the Supplementary DataS3.<br> 
## xtal_Yu2025 folder
In the folder 'xtal_Yu2025', there are five '.py' files and one '.csv' file.<br> 
### data file
The '.csv' file named 'parental_magma_data.csv' provides users with starting compositions used in fractional crystallization modeling shown in figure S3.<br> 
### code files
#### wl1989stoich_2023.py, wl1989kdcalc_2023.py, wlState_2023.py, wl1989models_2023.py
These codes defined fucntions used in calculating melt and mineral (olivine, plagioclase, clinopyroxene) compositions for two types of crystallization: fractional crystallization and equilibrium crystallization. Compositions calculated include SiO2, TiO2, Al2O3, FeO(t), MgO, CaO, Na2O, K2O, P2O5, MnO, NiO for both melt and minerals, and Rb, Sr, Nb, Ba, La, Ce, Pr, Nd, Sm, Eu, Gd, Dy, Er, Yb, Lu, Ta, Pb, Th, U, and the totally incompatible component with bulk partition coefficients equal to 0 for the melt. Note that iron is treated as ferrous iron in the calculation and we assume a ratio of FeO/FeO(t) as 0.9 for mid-ocean ridge basalts in our modeling.<br> 
Fundamental algorithms are give by Weaver, J.S. and Langmuir, C.H., 1990. Calculation of phase equilibrium in mineral-melt systems. Computers & Geosciences, 16(1), pp.1-19. Partitioning coefficients for major and trace components are provided in wl1989kdcalc_2023.py and wl1989models_2023.py, respectively.<br> 
A brief introduction is commented at the beginning of each code. These codes will be called by 'xtal_main_2025.py'.<br> 
#### xtal_main_2025.py
This code calls all the functions defined for crystallization calculations. With given starting compositions and pressure, users will get liquid and mineral compositions for fractional and equilibrium crystallization.<br> 
Output is in a data frame format which can be saved as .csv or excel file.<br> 
See comments in codes for detailed explanation of each variable.<br> 
### References cited in codes
Bédard, J.H., 2023. Trace element partitioning coefficients between terrestrial silicate melts and plagioclase feldspar: Improved and simplified parameters. Geochimica et Cosmochimica Acta, 350, pp.69-86.<br> 
Davis, F.A., Humayun, M., Hirschmann, M.M. and Cooper, R.S., 2013. Experimentally determined mineral/melt partitioning of first-row transition elements (FRTE) during partial melting of peridotite at 3 GPa. Geochimica et Cosmochimica Acta, 104, pp.232-260.<br> 
Hill, E., Wood, B.J. and Blundy, J.D., 2000. The effect of Ca-Tschermaks component on trace element partitioning between clinopyroxene and silicate melt. Lithos, 53(3-4), pp.203-215.<br> 
Langmuir, C.H., Klein, E.M. and Plank, T., 1992. Petrological systematics of mid-ocean ridge basalts: constraints on melt generation beneath ocean ridges. Geophysical monograph series, 71, pp.183-280.<br> 
LaTourrette, T.Z. and Burnett, D.S., 1992. Experimental determination of U and Th partitioning between clinopyroxene and natural and synthetic basaltic liquid. Earth and Planetary Science Letters, 110(1-4), pp.227-244.<br> 
Le Roux, V., Dasgupta, R. and Lee, C.T., 2011. Mineralogical heterogeneities in the Earth's mantle: Constraints from Mn, Co, Ni and Zn partitioning during partial melting. Earth and Planetary Science Letters, 307(3-4), pp.395-408.<br> 
Salters, V.J. and Stracke, A., 2004. Composition of the depleted mantle. Geochemistry, Geophysics, Geosystems, 5(5).<br> 
Sobolev, A.V., Hofmann, A.W., Sobolev, S.V. and Nikogosian, I.K., 2005. An olivine-free mantle source of Hawaiian shield basalts. Nature, 434(7033), pp.590-597.<br> 
St C. O’Neill, H. and Jenner, F.E., 2012. The global pattern of trace-element distributions in ocean floor basalts. Nature, 491(7426), pp.698-704.<br> 
Walter, M.J., 1998. Melting of garnet peridotite and the origin of komatiite and depleted lithosphere. Journal of petrology, 39(1), pp.29-60.<br> 
Weaver, J.S. and Langmuir, C.H., 1990. Calculation of phase equilibrium in mineral-melt systems. Computers & Geosciences, 16(1), pp.1-19.<br> 
Toplis, M.J., 2005. The thermodynamics of iron and magnesium partitioning between olivine and liquid: criteria for assessing and predicting equilibrium in natural and experimental systems. Contributions to Mineralogy and Petrology, 149(1), pp.22-39.<br> 
Yu, M. and Langmuir, C.H., 2023. The origin of Ni and Mn variations in Hawaiian and MORB olivines and associated basalts. Chemical Geology, 640, p.121745.<br> 
# Updates and Cite Policy
The updates of the code will be posted timely. Comments and suggestions are welcome and can be sent to Mingzhen Yu (myu@g.harvard.edu).<br> 
Any publications using this code need to follow the CITATION.cff to cite the codes.<br> 













