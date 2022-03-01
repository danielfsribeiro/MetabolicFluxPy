# MetabolicFluxPy

## Description
A package to perform flux balance analysis (FBA) on RNA-seq data.  

## Dependencies
Packages:  
pandas, numpy, cobrapy, gurobipy, scipy.

## Rationale
### Data preparation
The RNA-seq and empiral flux tables are imported. The genome-scale metabolic model (GEM) written in SMBL format is loaded.  
Genes names form the RNA-seq table are matched with existing genes on the metabolic model.  
[cobrapy][1] (1) is used to load, save, modify and query metabolic models before and after FBA.  
Data preparation, FBA analysis and possible follow up analysis, are based on [Granata, 2019][2] (2) paper. Particularly, [Lee-12][3] (3) algorithm is used here.  

### FBA analysis
[smop][4] (4) was used to automtic translate Matlab-to-Python code of [call_lee12][5] (5) and [Lee-12][3] algorithm (3). The generated code was verified and further and adapted to the this context. A modified `smop` run-time library `libsomp` is used in this package and is necessary to run these code.
`call_lee12`begins the analysis. `easyLP` and  `solveCobraLP` do the FBA heavy lifting.  
The [Cobra Toolbox][6] (6) was also written in Matlab. The minimum essential functions from the Cobra Toolbox, particularly those required to execute `solveCobraLP`, were hand converted from Matlab to Python. In `solveCobraLP`, only the Gurobi solver part is implemented, therefore this package requires [Gurobi][7] (7).  
The package, ran with the human GEM Recon3D (quite big), takes around 1.5 days to run on a laptop. 

### Output
A table is generated with the model's **Reaction Name**, **Reaction ID** and calculated **Flux** values.
In addition, a `json` and `xml` version of the obtained metabolic model is saved.
The output folder may be configured.

## Usage
```console
~$ python metabolic_flux.py -f FILE -g GENE_COL -e EXPRESSION_COL -s SD_COL -F MEASURED_FLUXES -m MODEL
                            [-o OUTPUT_PATH] [-n NAME] [-h]
```

Arguments:  
**-h, --help**            show this help message and exit.  
**-f *FILE*, --file *FILE***  path/to/rnaseq/file. Requires a Gene per line (gene column), a Mean expression (norm counts) column and a std dev column.  
**-g *GENE_COL*, --gene_col *GENE_COL*** Gene column name form *FILE*. Build with **gene symbols**.  
**-e *EXPRESSION_COL*, --expression_col *EXPRESSION_COL*** Expression column name.  
**-s *SD_COL*, --sd_col *SD_COL*** Std dev column name.  
**-F *MEASURED_FLUXES*, --measured_fluxes *MEASURED_FLUXES***  path/to/measured/fluxes/file. Format: First column are reaction names in the chosen metabolic
                      									 model. Second column is the measured flux. By convention, exchange reactions are written as
                      									 export reactions (e.g. ‘glc[e] <==>’), so import of a metabolite is a negative flux and
                      									 export of a metabolite is a positive flux.  
**-m MODEL, --model MODEL** path/to/metabolic/model/file, in XML format. Tested with Recon3D.  
**-o OUTPUT_PATH, --output_path OUTPUT_PATH** path/to/output/folder. Default is current path.
**-n NAME, --name NAME** Optional name for the analysis.

## Licence:  
[Granata, 2019][1], [Lee et al. 2012][2], [Machado et al. 2014][3], : CC BY 4.0 licence  
[OpenCobra Toolbox v.3.0][4]: GPL-3.0 licence  
[COBRApy 0.21][5]: GPL-2.0+ license  
[smop][6]: MIT licence  
[gurobi 9.5.1][7]: Individual academic licence  
  

## References
[1]: <https://github.com/opencobra/cobrapy> "COBRApy"  
[1]: COBRApy. (https://github.com/opencobra/cobrapy;https://doi.org/10.1186/1752-0509-7-74).  

[2]:  <https://doi.org/10.1186/s12859-019-2685-9> "Granata et al. 2019, BMC Bioinformatics"  
[2]:  Granata et al. 2019, BMC Bioinformatics. (https://doi.org/10.1186/s12859-019-2685-9).  

[3]: <https://doi.org/10.1186/1752-0509-6-73> "Lee et al. 2012, BMC Syst. Biol."  
[3]: Lee et al. 2012, BMC Syst. Biol. (https://doi.org/10.1186/1752-0509-6-73).  

[4]: <https://github.com/victorlei/smop> "smop - Small Matlab to Python compiler."  
[4]: smop - Small Matlab to Python compiler. (https://github.com/victorlei/smop)  

[5]: <https://doi.org/10.1371/journal.pcbi.1003580> "Machado et al. 2014, PLoS Comput. Biol."  
[5]: Machado et al. 2014, PLoS Comput. Biol. (https://doi.org/10.1371/journal.pcbi.1003580).  

[6]: <https://doi.org/10.1038/s41596-018-0098-2> "Heirendt et al. 2019, Nat Protoc."  
[6]: Heirendt et al. 2019, Nat Protoc. (https://doi.org/10.1038/s41596-018-0098-2). (OpenCobra Toolbox v.3.0)  

[7]: <https://www.gurobi.com> "Gurobi solver."  
[7]: Gurobi solver. (https://www.gurobi.com)  

