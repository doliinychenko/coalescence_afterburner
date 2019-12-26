# Coalescence afterburner
Constructing light nuclei from nucleons out of transport code

Currently the intention is the following:
1. Read extended SMASH output (see https://smash-transport.github.io/)
2. Loop over nucleon pairs, triplets, quadruplets: boost the to their CM frame and try to combine into light nuclei
3. Decay unstable light nuclei
4. Print out the light nuclei and their spectra

The final goal so far is to test
1. How well one can describe deuteron? Can I do it as well as Sombun et al. (http://inspirehep.net/record/1675249)?
2. How well one can describe tritons, helium-3, helium-4, and hypertriton following the same ideas?
