# Slow extraction with crystals

The aim is to simulate the Slow  Extraction from the SPS towards the Noth Area. 

In particular this study aims to implement a new scheme for slow extraction which is non-resonant. This technique uses TECA crystal to give a kick to the particles that will be extracted. The beam is slowly pushed towards the crystal via reference momentum change.

To run these codes is mandatory to have installed Xsuite.
You can install it in the following link:

https://xsuite.readthedocs.io/en/latest/installation.html

It is highly reccomended to create an appropriate environment, 
in my case I did it with the miniforge that you can find in the link:

https://xsuite.readthedocs.io/en/latest/installation.html#install-miniforge

or simly use the requirements.txt file to replicate the environment. 

The best example of how the EverestCrystal is used is in the directory Non-Resonant Extraction where there is the file RecordingInteraction.ipynb. In this file you can see how the internal logger is used to save the interactions of the particles with the crystal.