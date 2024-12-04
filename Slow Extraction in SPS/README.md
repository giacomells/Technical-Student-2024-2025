To run the model (ModelFinalDoubleLine.ipynb) you have to download the dictionaty (septa.py) and the original sequence (sps_with_extraction_sliced_quads.json).

In this model a line is created for the circulating beam and there are some particles that get extracted if they hit an aperture(septum = xt.LimitRect(...)).
Once these particles are extracted a copy of them is created and a new bunch of particles is created with the same properties of the extracted particles.

A new line then is created (LSS2) where the septum is installed with the proper functions.
Every element of the septum is sliced in many smaller parts and every part is created as a BeamInteraction object which has an
interaction_process = SeptumInteraction (....)
Every particle that is outside the aperture of this element gets lost.

Once the new line LSS2 is created with every element installed, the extracted particles are tracked trhough the LSS2 in only one turn.