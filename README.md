# sia-upsurging


## Description
This is a python implementation of the UPstream SURface Gradient cappING (UPSURGING) scheme which inhibits mass conservation violation (negative ice thicknesses) in shallow ice approximation (SIA) ice flow models based on finite differences. The key assets of this method are its simple implementation and low computational cost. For more details see Appendix C of the dissertation of Michael Andreas Imhof, "Combined climate-ice flow modelling of the Alpine Ice Field during the Last Glacial Maximum", 2021, ETH Zurich, dissertation no. 27416, DOI:???


## Requirements
The code is compatible with Python 2 and 3. The following packages are required:
- numpy
- time
- copy
- os.path
- matplotlib


## Script list
- **'launcher.py'** contains an example of how to start the model
- **'sia-upsurging.py'** contains the code for an SIA model with the option to use the UPSURGING scheme or the muscl SUPERBEE flux limiter of [Jarosch et al., (2013)](https://doi.org/10.5194/tc-7-229-2013). 


## sia-upsurging.py
The model loads the bedrock topography of Rhone Valley at a resolution of 1x1 km (1km-rhone_valley.txt) and uses an elevation dependent mass balance model. A directory with a time stamp and a name is created into which ice thickness data is written at a regular base. A figure from the end state is also saved in there as well as performance data. 
Three methods to calculate the ice flow by the SIA are available (-model_choice):
- **m2** This is the SIA default version which is typically used. It is prone to problems with mass conservation. 
- **upsurging** The surface gradient capping method introduced in Appendix C of _link to my dissertation_. It resolves the mass conservation problem. 
- **muscl** Flux limiter of [Jarosch et al., (2013)](https://doi.org/10.5194/tc-7-229-2013). This is an other method to resolve the mass conservation problem. 


### Eismint
The model has also the option to perform two [Eismint benchmark tests](https://www.cambridge.org/core/journals/annals-of-glaciology/article/eismint-benchmarks-for-testing-icesheet-models/F8563050E59F7161FAD3EA55329E70E6). To do so, set the mass balance model (-smbm) to eismint1fm (fixed margin) or eismint1mm (moving margin) and provide no input bedrock topography (omit -ib). 


## How to run
To start the model type "./launcher.py" in the terminal. 


## TODO
- list with input parsers
- describe output data in more detail
- add reference to bedrock topography
