This repository holds the code that I used for reproducing the main findings of the single cell analysis from the following paper:

https://stm.sciencemag.org/content/early/2020/11/25/scitranslmed.abe4282

The analysis comprises of analysing single cells from a Covid19 patient pre and post transplant, a healthy lung donor and a second Covid patient.
This analysis showed that when postmortem lung tissue from COVID19 patients were batch corrected and integrated with same tissue type of fibrosis patients and controls, idiopathic pulmonary fibrosis patients showed similarity of epithelial, immune and stromal cells with COVID19 patients and while control lungs from the two studies shared similarity in the latent state after integration of the two datasets. This becomes clear when looking at the Habermann_Bharat_all_tissue_type_integrated_obs_UMAP_By_Tissue_type.pdf plot within the 'Integration with Habermann' folder. 

The Bharat et. al. dataset was analysed individually and then subsequently that dataset was used for latent space analysis of Bharat et. al. and Habermann et. al. datasets.

Paper data used for integration with Bharat data: https://www.science.org/doi/10.1126/sciadv.aba1972?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed

The code used in the 'src' folder.
