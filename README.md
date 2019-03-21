# prf-explainer

Interactive (matlab) code illustrating the population receptive field mapping technique (Dumoulin &amp; Wandell)

## To get going

1. ``git clone`` this repository (assuming macOS). Use Windows magic to achieve the same:
```bash
cd ~
git clone https://github.com/schluppeck/prf-explainer.git
```

2. In `matlab`, open the live script
```MATLAB
cd ~/prf-explainer
open('pRFfittingExample')
```

3. Hit the play button and walk through the different cells:
<center>
<img src="./screenshot.png" width="70%">
</center>

## Notes and references

Tested with MATLAB Version: 9.4.0.813654 (R2018a) on macOS. If you find errors, typos, better ways of doing things, submit and issue on the github repo.

- Serge Dumoulin's homepage (original paper, software) – http://www.spinozacentre.nl/dumoulin/
- Brian Wandell's homepage – https://web.stanford.edu/group/vista/cgi-bin/wandell/
- Justin Gardner's homepage – http://gru.stanford.edu/doku.php/shared/home
http://gru.stanford.edu/doku.php/mrTools/tutorialsprf – mrTools implementation (and tutorial for how to use)
- my *Data Analysis for Neuroimaging* lab / teaching materials are hosted here: https://schluppeck.github.io/dafni/
