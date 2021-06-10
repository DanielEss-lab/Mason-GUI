# Purpose
This GUI connects with [Mason](https://github.com/DanielEss-lab/Mason) which connects to [Taylor](https://github.com/DanielEss-lab/Taylor).
# Installation and Setup
Do the following in miniconda:
```bash
conda create -n gui
```
```bash
conda activate gui
```
```bash
conda install -c conda-forge openbabel
```
```bash
pip install PyQt5 pyqtgraph PyOpenGL
```

Download this repository to your computer. In the miniconda command prompt, cd to where you downloaded it.

After everything is set up, run
```bash
python newgui.py
```
This should open up the gui in a separate window.