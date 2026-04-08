## Animations

This directory contains scripts for generating SPS extraction animations and for saving the SPS sequence used by the animation workflow.

### Files

- `save_sequence_SPS.py`: builds and saves `sps_for_sx.json`.
- `phaseSpaceAnimation.py`: generates a phase-space animation from a local SPS JSON sequence.


### Usage notes

- `phaseSpaceAnimation.py` expects `sps_for_sx.json` to be present in this same directory.
- Run `save_sequence_SPS.py` first if the JSON sequence has not been generated yet.
- `save_sequence_SPS.py` depends on CERN model files fetched over the network.

To run this code you need to run the file save_sequence_SPS.py to store the SPS sequence into your directory. In this case the sequence is saved in sps_for_sx.json 
Then you can run the phaseSpaceAnimation.py file properly.