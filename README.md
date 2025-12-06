# STPC_Lite
A Lite version of Spatio-Temporal Principal Component (STPC) filter for basin-scale hydrological, cryospheric and oceanographic applications of GRACE/GRACE-FO. By reconciling two eigenvalue problems in spherical Slepian method and Multichannel Singular Spectrum Analysis (MSSA) based on different noise significance levels, The STPC filter establishes a statistical framework for basin-scale analyses via automatically determined, region-specific physical parameters governing spatio-temporal localization, leakage-in/out effects differentiation, uncertainty minimization across data centers, and mitigation of north-south striping.

This version incorporates additional functionality for ocean applications (e.g., recovering GAD and correcting Inverted barometer in [Slepian_ocean_add-on package](https://github.com/KMartin0013/Slepian_ocean_add-on/tree/master)). 

## Notice before use
1. The Slepian procedure of this code is based on the [foundational Software from the Simons Laboratories](https://geoweb.princeton.edu/people/simons/software.html), and [its add-on version](https://github.com/KMartin0013/Slepian_ocean_add-on) for ocean application (e.g., GAD, GIA, IB correction).
2. It is strongly recommended to familiarize yourself with the foundational work by <a href="https://polarice.geo.arizona.edu/">C. Harig</a> &amp; <a href="http://www.frederik.net">F. J. Simons</a> ([Ref 1](https://doi.org/10.1073/pnas.1206785109)), and our previous work on ocean ([Ref 2](https://doi.org/10.1016/j.jag.2024.104065)) before using this add-on code.
3. **As the code has not yet been fully optimised, it is advisable to reserve 5–10GB of space before running each case.**


**Required software:**<br>
[slepian_alpha](https://github.com/csdms-contrib/slepian_alpha)  
[slepian_bravo](https://github.com/csdms-contrib/slepian_bravo)  
[slepian_delta](https://github.com/csdms-contrib/slepian_delta)    
[slepian_zero](https://github.com/csdms-contrib/slepian_zero) (actually only guyotphysics.m)  

**Required software for mapping:**<br>
[m_map](https://www.eoas.ubc.ca/~rich/map.html)  
[tight_subplot](https://ch.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)  
[hatchfill2](https://ch.mathworks.com/matlabcentral/fileexchange/53593-hatchfill2)  
[special heatmap](https://ch.mathworks.com/matlabcentral/fileexchange/125520-special-heatmap/)

## Usage of this package
1. Download the required software (or codes) and this codes as mentioned above.
2. Download the Spherical Harmonic coefficients and the associated correction files (e.g., degree 1, degree 20 correction) from different data centers. 
  2.1 Notice that for now this version supports the Release 05 (RL05) and Release 06 (RL06) data from the Center of Space Research (CSR) at the University of Texas, the Jet Propulsion Laboratory (JPL), GeoForschungsZentrum (GFZ), and the Institute of Geodesy at the Graz University of Technology (ITSG). If you use data centers other than these four, please revise the function 'grace2plmt_m.m' and 'gracedeg1_m' accordingly.  
3. Run 'Lite_Case1_Greenland.m' for a case study in the Greenland Ice Sheet. Three main procedures were included: Lite_Slepian_V1.m, Lite_MSSA_V1.m, Lite_STPC_Greenland_V1.m.   
  3.1 The fist procedure 'Lite_Slepian_V1.m' calculates spherical Slepian functions and corresponding spherical Slepian coefficients, primarily following the guideline in [slepian_delta](https://github.com/csdms-contrib/slepian_delta). Notice that the buffer zone for the study of interest (parameter R) should be first defined (normally by pre-simulating). You may need to run this procedure several times if you want to remove the uncertainty of GRACE/GRACE-FO across different data centers.
  3.2 The second procedure 'Lite_MSSA_V1.m' decomposes the corresponding spherical Slepian coefficients by MSSA and conducts a noise detection process at given significance level p. Notice that the sliding window size (parameter M) is optimally selected and a iterative MSSA procedure for gap filling proposed by <a href="https://www.linkedin.com/in/louis-marie-gauer-8a99a7176/?originalSubdomain=fr/">L. Gauer</a> et al. ([Ref 3](https://doi.org/10.1029/2022JB025561)) is automatically conducted in this procedure .
   3.3 The third procedure 'Lite_STPC_Greenland_V1.m' determines the optimal truncation number S for spherical Slepian method and cutoff number N for MSSA. Then, the area-weighted time series and spatial maps of equivalent water height (EWH), or non-steric sea level anomalies (NSLA), or mass loss at different significance level p was provided.
4. Case study for South China Sea is also provided in 'Lite_Case1_SCS.m'.  
   
## Citation Information:
Please cite our work and the foundational work by C. Harig, F. J. Simons, and L. Gauer as appropriate if you find this package useful.  

[Ref 1] Harig, Christopher and Frederik J. Simons. 
Mapping Greenland's mass loss in space and time.
<i>Proc. Natl. Acad. Sc.</i>, 109(49), 19934-19937, 2012.
doi:https://doi.org/10.1073/pnas.1206785109

[Ref 2] Zhongtian, Ma, Hok Sum Fok, Robert Tenzer, Jianli Chen. 
A Novel Slepian Approach for Determining Mass-term Sea Level from GRACE over the South China Sea. <i>Int. J. Appl. Earth. Obs. Geoinf.</i>, 132, 104065, 2024.
doi:https://doi.org/10.1016/j.jag.2024.104065

[Ref 3] L. Gauer, K. Chanard, L. Fleitout. 
Data‐driven gap filling and spatio‐temporal filtering of the GRACE and GRACE‐FO records. <i>J. Geophys. Res. Solid Earth</i>, 128(5), e2022JB025561, 2023. doi:https://doi.org/10.1029/2022JB025561
