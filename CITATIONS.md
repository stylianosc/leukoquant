# Citations

LeukoQuant builds on several open-source tools and published methods. If you use LeukoQuant in your research, please also cite the relevant underlying tools listed below.

---

## FreeSurfer

Used for cortical reconstruction, subcortical segmentation, and brain parcellation (`process-recon`).

- Fischl B (2012). FreeSurfer. *NeuroImage*, 62(2):774–781. https://doi.org/10.1016/j.neuroimage.2012.01.021

---

## SynthStrip

Used for skull-stripping (`--skull-strip` flag in `process-dti` and `process-noddi`).

- Hoopes A, Mora JS, Dalca AV, Fischl B, Hoffmann M (2022). SynthStrip: skull-stripping for any brain image. *NeuroImage*, 260:119474. https://doi.org/10.1016/j.neuroimage.2022.119474

---

## TRACULA

Used for probabilistic tractography (`process-tracula`).

- Yendiki A, Panneck P, Srinivasan P, Stevens A, Zöllei L, Augustinack J, Wang R, Salat D, Ehrlich S, Behrens T, Jbabdi S, Gollub R, Fischl B (2011). Automated probabilistic reconstruction of white-matter pathways in health and disease using an atlas of the underlying anatomy. *Frontiers in Neuroinformatics*, 5:23. https://doi.org/10.3389/fninf.2011.00023

---

## FSL

Used for DWI pre-processing, registration, and brain extraction.

- Smith SM, Jenkinson M, Woolrich MW, Beckmann CF, Behrens TEJ, Johansen-Berg H, Bannister PR, De Luca M, Drobnjak I, Flitney DE, Niazy R, Saunders J, Vickers J, Zhang Y, De Stefano N, Brady JM, Matthews PM (2004). Advances in functional and structural MR image analysis and implementation as FSL. *NeuroImage*, 23(Suppl 1):S208–S219. https://doi.org/10.1016/j.neuroimage.2004.07.051
- Jenkinson M, Beckmann CF, Behrens TEJ, Woolrich MW, Smith SM (2012). FSL. *NeuroImage*, 62:782–790. https://doi.org/10.1016/j.neuroimage.2011.09.015

---

## GIF (Geodesic Information Flows)

Used for brain parcellation and segmentation (`process-gif`).

- Cardoso MJ, Modat M, Wolz R, Melbourne A, Cash D, Rueckert D, Ourselin S (2015). Geodesic information flows: spatially-variant graphs and their application to segmentation and fusion. *IEEE Transactions on Medical Imaging*, 34(9):1976–1988. https://doi.org/10.1109/TMI.2015.2418298

---

## BaMoS

Used for white matter hyperintensity segmentation (`process-bamos`).

- Sudre CH, Cardoso MJ, Bouvy WH, Biessels GJ, Barnes J, Ourselin S (2015). Bayesian model selection for pathological neuroimaging data applied to white matter lesion segmentation. *IEEE Transactions on Medical Imaging*, 34(10):2079–2102. https://doi.org/10.1109/TMI.2015.2419072

---

## NiftyReg

Used for image registration within the pipeline.

- Modat M, Cash DM, Daga P, Winston GP, Duncan JS, Ourselin S (2014). Global image registration using a symmetric block-matching approach. *Journal of Medical Imaging*, 1(2):024003. https://doi.org/10.1117/1.JMI.1.2.024003
- Modat M, Ridgway GR, Taylor ZA, Lehmann M, Barnes J, Hawkes DJ, Fox NC, Ourselin S (2010). Fast free-form deformation using graphics processing units. *Computer Methods and Programs in Biomedicine*, 98(3):278–284. https://doi.org/10.1016/j.cmpb.2009.09.002

---

## NODDI

Used for NODDI model fitting (`process-noddi`).

- Zhang H, Schneider T, Wheeler-Kingshott CA, Alexander DC (2012). NODDI: Practical in vivo neurite orientation dispersion and density imaging of the human brain. *NeuroImage*, 61(4):1000–1016. https://doi.org/10.1016/j.neuroimage.2012.03.072
- Daducci A, Canales-Rodríguez EJ, Zhang H, Dyrby TB, Alexander DC, Thiran JP (2015). Accelerated Microstructure Imaging via Convex Optimization (AMICO) from diffusion MRI data. *NeuroImage*, 105:32–44. https://doi.org/10.1016/j.neuroimage.2014.10.026

---

## NiftySeg

Used for label fusion and segmentation within BaMoS.

- Cardoso MJ, Leung K, Modat M, Keihaninejad S, Cash D, Barnes J, Fox NC, Ourselin S; ADNI (2013). STEPS: Similarity and truth estimation for propagated segmentations and its application to hippocampal segmentation and brain parcellation. *Medical Image Analysis*, 17(6):671–684. https://doi.org/10.1016/j.media.2013.02.006

---

## Snakemake

Used as the workflow engine.

- Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, Forster J, Lee S, Twardziok SO, Kanitz A, Wilm A, Holtgrewe M, Rahmann S, Nahnsen S, Köster J (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10:33. https://doi.org/10.12688/f1000research.29032.2

---

## Apptainer / Singularity

Used for containerised execution of all compute steps.

- Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. *PLoS ONE*, 12(5):e0177459. https://doi.org/10.1371/journal.pone.0177459

---

## dcm2niix

Used for DICOM to NIfTI conversion.

- Li X, Morgan PS, Ashburner J, Smith J, Rorden C (2016). The first step for neuroimaging data analysis: DICOM to NIfTI conversion. *Journal of Neuroscience Methods*, 264:47–56. https://doi.org/10.1016/j.jneumeth.2016.03.001

---

## Python Libraries

Core scientific Python libraries used throughout the pipeline.

- **NumPy** - Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, Wieser E, Taylor J, Berg S, Smith NJ, Kern R, Picus M, Hoyer S, van Kerkwijk MH, Brett M, Haldane A, del Río JF, Wiebe M, Peterson P, Gérard-Marchant P, Sheppard K, Reddy T, Weckesser W, Abbasi H, Gohlke C, Oliphant TE (2020). Array programming with NumPy. *Nature*, 585:357–362. https://doi.org/10.1038/s41586-020-2649-2
- **SciPy** - Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, Burovski E, Peterson P, Weckesser W, Bright J, van der Walt SJ, Brett M, Wilson J, Millman KJ, Mayorov N, Nelson ARJ, Jones E, Kern R, Larson E, Carey CJ, Polat İ, Feng Y, Moore EW, VanderPlas J, Laxalde D, Perktold J, Cimrman R, Henriksen I, Quintero EA, Harris CR, Archibald AM, Ribeiro AH, Pedregosa F, van Mulbregt P, SciPy 1.0 Contributors (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17:261–272. https://doi.org/10.1038/s41592-019-0686-2
- **pandas** - McKinney W (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61. https://doi.org/10.25080/Majora-92bf1922-00a
- **Matplotlib** - Hunter JD (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3):90–95. https://doi.org/10.1109/MCSE.2007.55
- **scikit-image** - van der Walt S, Schönberger JL, Nunez-Iglesias J, Boulogne F, Warner JD, Yager N, Gouillart E, Yu T, scikit-image contributors (2014). scikit-image: image processing in Python. *PeerJ*, 2:e453. https://doi.org/10.7717/peerj.453
- **scikit-learn** - Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, Blondel M, Prettenhofer P, Weiss R, Dubourg V, Vanderplas J, Passos A, Cournapeau D, Brucher M, Perrot M, Duchesnay É (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12:2825–2830. http://jmlr.org/papers/v12/pedregosa11a.html
- **nilearn** - Abraham A, Pedregosa F, Eickenberg M, Gervais P, Mueller A, Kossaifi J, Gramfort A, Thirion B, Varoquaux G (2014). Machine learning for neuroimaging with scikit-learn. *Frontiers in Neuroinformatics*, 8:14. https://doi.org/10.3389/fninf.2014.00014
- **dipy** - Garyfallidis E, Brett M, Amirbekian B, Rokem A, van der Walt S, Descoteaux M, Nimmo-Smith I (2014). Dipy, a library for the analysis of diffusion MRI data. *Frontiers in Neuroinformatics*, 8:8. https://doi.org/10.3389/fninf.2014.00008
