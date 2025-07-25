# Slow extraction with crystals

The aim is to simulate the Slow Resonant Extraction from the SPS into the LSS2 towards the Noth Area. 
In particular we need to make this simulation with crystals to see if the efficiency is higher than the ZS septa.

To run these codes is mandatory to have installed Xsuite.
You can install it in the following link:

https://xsuite.readthedocs.io/en/latest/installation.html

It is highly reccomended to create an appropriate environment, 
in my case I did it with the miniforge that you can find in the link:

https://xsuite.readthedocs.io/en/latest/installation.html#install-miniforge

The best example of how the EverestCrystal is used is in the directory Non-Resonant Extraction where there is the file RecordingInteraction.ipynb. In this file you can see how the internal logger is used to save the interactions of the particles with the crystal.